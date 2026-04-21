%% ========= FIG3E_ABS_ZSCORE_EARLY_VS_LATE__GLOBALWARP__FIG2TRIALS_EXACT.m =========
% NOW (OPTION A: NO WARPING) + SESSION-LEVEL ANALYSIS:
% 1) Uses ONLY CLASSIFIED units (cell_type in {1,2,3}) (empties and 0=Unclassified skipped)
% 2) Loads GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat ROBUSTLY (auto-searches common paths)
% 3) Uses the saved GLOBAL window endpoints (cue/press/lick) as LAG WINDOWS RELATIVE TO EACH EVENT
%    on EACH TRIAL (NO global-warp timeline):
%       Epoch 1 = Cue window   (relative to cue=0 on each trial)
%       Epoch 2 = Press window (relative to press time on each trial)
%       Epoch 3 = Lick window  (relative to lick time on each trial)
%
% IMPORTANT:
% This script does NOT re-fit the GLM; it only uses GLOBAL_WINDOWS to pick time intervals
% for computing mean firing rate. Windows are applied per-trial around each event time (Option A).
%
% SESSION-LEVEL CHANGE (to avoid pseudoreplication):
% - Compute per-unit metrics as before (trial -> unit mean).
% - Then average across units WITHIN EACH SESSION to yield one value per epoch per session.
% - Use those session-level means as independent observations for Early vs Late comparisons.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/3_GENERAL_CHANGE_EARLY_LATE';
outPng  = fullfile(baseOut, '3_GENERAL_CHANGE_EARLY_LATE_FIRING_RATE.png'); % (name kept as requested)

% >>> ADDED: also save SVG and MAT <<<
[outDir, outBase, ~] = fileparts(outPng);
outSvg = fullfile(outDir, [outBase '.svg']);
outMat = fullfile(outDir, [outBase '.mat']);

% >>> REQUIRED: load classification output <<<
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skipped)

% Analysis window used to BUILD SDFs (kept, but we will use an expanded fixed window below)
preCueMs   = 500;
postLickMs = 3000;

% FIG2-like trial selection window (relative to cue) — MUST match your code 2/3C
fixedWin = [-1000, 6000];

% Trial filters (must match your code 2/3C)
MIN_RT_MS    = 100;
requireValid = true;

% PER-UNIT inclusion (must match your code 2/3C)
minTrialsPerUnit = 10;

% SDF settings (paper-style)
dt_ms        = 10;
gaussSigmaMs = 25;

% Groups
earlyRange = [1 8];
lateRange  = [29 36];

% Plot style
figPos  = [120, 120, 1050, 520];

% >>> DOUBLED (as requested) <<<
fsTitle = 40;    % was 20
fsAx    = 32;    % was 16
lwAx    = 5.0;   % was 2.5

% Colors (only for points/boxes)
colEarly = [0.10 0.60 0.15];
colLate  = [0.90 0.70 0.10];

rng(0);

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

%% ---- PRECOMPUTE GAUSS KERNEL (unit area; Hz after dividing by dt) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- LOAD GLOBAL WINDOWS (ROBUST PATH SEARCH) ----
candidates = { ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(fileparts(outPng),'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
    };

winFile = '';
for i = 1:numel(candidates)
    if exist(candidates{i},'file') == 2
        winFile = candidates{i};
        break;
    end
end

assert(~isempty(winFile), ...
    'Missing GLOBAL windows file. Looked in:\n  %s', strjoin(candidates, sprintf('\n  ')));

W = load(winFile);
assert(isfield(W,'GLOBAL_WINDOWS') && isstruct(W.GLOBAL_WINDOWS), ...
    'GLOBAL_WINDOWS missing in: %s', winFile);
GW = W.GLOBAL_WINDOWS;

fprintf('\nLoaded GLOBAL windows from:\n  %s\n', winFile);

% Windows are stored as [start end] in ms RELATIVE TO EACH EVENT (lag axis).
wCue_rel   = GW.cue_ms;    % [ms, ms] relative to cue
wPress_rel = GW.press_ms;  % [ms, ms] relative to press
wLick_rel  = GW.lick_ms;   % [ms, ms] relative to lick

fprintf('\n=== WINDOWS (RELATIVE TO EACH EVENT) ===\n');
fprintf('Cue window (rel cue):     [%d, %d] ms\n', round(wCue_rel(1)),   round(wCue_rel(2)));
fprintf('Press window (rel press): [%d, %d] ms\n', round(wPress_rel(1)), round(wPress_rel(2)));
fprintf('Lick window (rel lick):   [%d, %d] ms\n', round(wLick_rel(1)),  round(wLick_rel(2)));

%% ---- BUILD A SINGLE CUE-ALIGNED TIME AXIS THAT GUARANTEES WINDOW COVERAGE ----
minLag = min([wCue_rel(1), wPress_rel(1), wLick_rel(1)]);
maxLag = max([wCue_rel(2), wPress_rel(2), wLick_rel(2)]);

winLeft  = fixedWin(1) + minLag;
winRight = fixedWin(2) + maxLag;

winLeft  = min(winLeft, -preCueMs);
winRight = max(winRight, fixedWin(2));

tgrid_ms = winLeft:dt_ms:winRight; %#ok<NASGU>
nT = numel(winLeft:dt_ms:winRight);

fprintf('\n=== CUE-ALIGNED SDF AXIS USED FOR ALL TRIALS (Option A) ===\n');
fprintf('SDF axis: [%d, %d] ms (dt=%d ms, nT=%d)\n', round(winLeft), round(winRight), dt_ms, nT);

%% ---- COMPUTE SESSION-LEVEL METRICS FOR EACH GROUP (NO WARP) ----
M_early = computeAbsZMetricsForSessions_(earlyRange, beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, ...
    winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel);

M_late  = computeAbsZMetricsForSessions_(lateRange, beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, ...
    winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel);

fprintf('\nEarly sessions: %d | Late sessions: %d\n', size(M_early,1), size(M_late,1));

%% ---- PLOT (Fig3E-like: Early vs Late for each window) ----
fig = figure('Color','w','Position',figPos);

epochNames = { ...
    sprintf('Cue',   round(wCue_rel(1)),   round(wCue_rel(2))), ...
    sprintf('Lever press', round(wPress_rel(1)), round(wPress_rel(2))), ...
    sprintf('Reward lick',  round(wLick_rel(1)),  round(wLick_rel(2)))};

fprintf('\n=== Rank-sum p-values (Early vs Late; session-level) ===\n');

for e = 1:3
    ax = subplot(1,3,e); hold(ax,'on');

    yE = M_early(:,e);
    yL = M_late(:,e);

    boxchart(ax, ones(size(yE)), yE, 'BoxFaceColor', colEarly, 'MarkerStyle','none');
    boxchart(ax, 2*ones(size(yL)), yL, 'BoxFaceColor', colLate,  'MarkerStyle','none');

    jitter = 0.18;
    scatter(ax, 1 + (rand(size(yE))-0.5)*2*jitter, yE, 28, 'filled', ...
        'MarkerFaceColor', colEarly, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 2 + (rand(size(yL))-0.5)*2*jitter, yL, 28, 'filled', ...
        'MarkerFaceColor', colLate,  'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);

    p = ranksum(yE, yL);

    fprintf('%s: p = %.6g\n', epochNames{e}, p);

    stars = '';
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    end

    if ~isempty(stars)
        yMax = max([yE(:); yL(:)]);
        yMin = min([yE(:); yL(:)]);
        yRange = yMax - yMin;
        if ~isfinite(yRange) || yRange == 0
            yRange = 1;
        end

        yBar = yMax + 0.08*yRange;
        hBar = 0.03*yRange;

        % >>> DOUBLED (as requested): LineWidth 2 -> 4 <<<
        plot(ax, [1 1 2 2], [yBar yBar+hBar yBar+hBar yBar], 'k-', 'LineWidth', 4);

        text(ax, 1.5, yBar+hBar+0.01*yRange, stars, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', fsTitle, 'FontWeight','bold');
    end

    title(ax, sprintf('%s', epochNames{e}), 'FontSize', fsTitle, 'FontWeight','bold');
    set(ax, 'XLim',[0.4 2.6], 'XTick',[1 2], 'XTickLabel',{'Early','Late'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    ylabel(ax, 'Mean FR (Hz)');
    box(ax,'off');
end

% >>> DOUBLED offset too: fsTitle+2 -> fsTitle+4 <<<
sgtitle(sprintf('Event-window firing rates: Early vs Late', ...
    earlyRange(1), earlyRange(2), lateRange(1), lateRange(2)), ...
    'FontSize', fsTitle+4, 'FontWeight','bold');

saveas(fig, outPng);
fprintf('Saved: %s\n', outPng);

saveas(fig, outSvg);
fprintf('Saved: %s\n', outSvg);

save(outMat, 'fig', 'M_early', 'M_late', 'earlyRange', 'lateRange', ...
    'wCue_rel', 'wPress_rel', 'wLick_rel', 'fixedWin', 'dt_ms', 'gaussSigmaMs', ...
    'winLeft', 'winRight');
fprintf('Saved: %s\n', outMat);

%% ================= LOCAL HELPERS =================

function M = computeAbsZMetricsForSessions_(sessRange, beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, ...
    winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel)

    nSessions = numel(beh);

    s0 = max(1, min(nSessions, sessRange(1)));
    s1 = max(1, min(nSessions, sessRange(2)));
    if s1 < s0
        M = nan(0,3);
        return;
    end

    M = nan(0,3); % nSessionsUsed x 3

    for sIdx = s0:s1
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        keepTrial = false(numel(trials),1);
        cueAbs   = nan(numel(trials),1);
        pressLat = nan(numel(trials),1);
        lickLat  = nan(numel(trials),1);

        for k = 1:numel(trials)
            tr = trials(k);

            if requireValid
                if ~isfield(tr,'valid') || ~tr.valid, continue; end
            end
            if ~all(isfield(tr, {'cue','press','lick'})), continue; end

            cue   = double(tr.cue);
            press = double(tr.press);
            lick  = double(tr.lick);
            if any(~isfinite([cue press lick])), continue; end

            rt = press - cue;
            if rt < MIN_RT_MS, continue; end

            pr = press - cue;
            lr = lick  - cue;

            if ~(pr >= fixedWin(1) && pr <= fixedWin(2)), continue; end
            if ~(lr >= fixedWin(1) && lr <= fixedWin(2)), continue; end

            keepTrial(k) = true;
            cueAbs(k)    = cue;
            pressLat(k)  = rt;
            lickLat(k)   = lr;
        end

        idxKeep = find(keepTrial);
        if isempty(idxKeep), continue; end

        cueAbsK   = cueAbs(idxKeep);
        pressLatK = pressLat(idxKeep);
        lickLatK  = lickLat(idxKeep);

        M_units = nan(0,3);

        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end

            for u = 1:numel(uc)

                ct = safe_get_celltype_(cell_type, sIdx, ch, u);
                if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                    continue;
                end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                activeTrials = false(numel(idxKeep),1);
                spk = spk_abs;
                if ~issorted(spk), spk = sort(spk); end
                j = 1;
                nSpk = numel(spk);

                for iTr = 1:numel(idxKeep)
                    cue0 = cueAbsK(iTr);
                    w0 = cue0 + fixedWin(1);
                    w1 = cue0 + fixedWin(2);

                    while j <= nSpk && spk(j) < w0
                        j = j + 1;
                    end
                    if j <= nSpk && spk(j) <= w1
                        activeTrials(iTr) = true;
                    end
                end

                if nnz(activeTrials) < minTrialsPerUnit
                    continue;
                end

                cueU = cueAbsK(activeTrials);
                rtU  = pressLatK(activeTrials);
                lkU  = lickLatK(activeTrials);
                nTrU = numel(cueU);

                m1_tr = nan(nTrU,1);
                m2_tr = nan(nTrU,1);
                m3_tr = nan(nTrU,1);

                for iTr = 1:nTrU
                    cue0   = cueU(iTr);
                    tPress = rtU(iTr);
                    tLick  = lkU(iTr);

                    spk_rel = spk_abs - cue0;
                    spk_rel = spk_rel(spk_rel >= winLeft & spk_rel <= winRight);

                    counts = zeros(1, nT);
                    if ~isempty(spk_rel)
                        jj = round((spk_rel - winLeft)/dt_ms) + 1;
                        jj = jj(jj >= 1 & jj <= nT);
                        counts = accumarray(jj(:), 1, [nT 1], @sum, 0).';
                    end

                    sm = conv(counts, g, 'same');
                    y  = sm / (dt_ms/1000);

                    idxCue   = msWindowToIdx_( [0      + wCue_rel(1),   0      + wCue_rel(2)],   winLeft, dt_ms, nT );
                    idxPress = msWindowToIdx_( [tPress + wPress_rel(1), tPress + wPress_rel(2)], winLeft, dt_ms, nT );
                    idxLick  = msWindowToIdx_( [tLick  + wLick_rel(1),  tLick  + wLick_rel(2)],  winLeft, dt_ms, nT );

                    m1_tr(iTr) = mean(y(idxCue(1):idxCue(2)),     'omitnan');
                    m2_tr(iTr) = mean(y(idxPress(1):idxPress(2)), 'omitnan');
                    m3_tr(iTr) = mean(y(idxLick(1):idxLick(2)),   'omitnan');
                end

                m1 = mean(m1_tr, 'omitnan');
                m2 = mean(m2_tr, 'omitnan');
                m3 = mean(m3_tr, 'omitnan');

                M_units(end+1,:) = [m1 m2 m3]; %#ok<AGROW>
            end
        end

        if isempty(M_units)
            continue;
        end

        M_sess = mean(M_units, 1, 'omitnan');
        M(end+1,:) = M_sess; %#ok<AGROW>
    end
end

function idx12 = msWindowToIdx_(w_ms, winLeft_ms, dt_ms, nT)
    if any(~isfinite(w_ms))
        idx12 = [1 1];
        return;
    end
    i1 = round((w_ms(1) - winLeft_ms)/dt_ms) + 1;
    i2 = round((w_ms(2) - winLeft_ms)/dt_ms) + 1;
    i1 = max(1, min(nT, i1));
    i2 = max(1, min(nT, i2));
    if i2 < i1, tmp = i1; i1 = i2; i2 = tmp; end
    idx12 = [i1 i2];
end

function ct = safe_get_celltype_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
       ch   <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
       u    <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
    halfWidth = ceil(5*sigmaMs/dtMs);
    x = (-halfWidth:halfWidth) * dtMs;
    g = exp(-0.5*(x./sigmaMs).^2);
    g = g / sum(g);
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k=1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k})), beh = S.(cands{k}); return; end
    end
    f = fieldnames(S);
    for i=1:numel(f)
        v = S.(f{i});
        if isstruct(v) && numel(v)>1, beh = v; return; end
    end
    for i=1:numel(f)
        v = S.(f{i});
        if isstruct(v), beh = v; return; end
    end
end

function tf = isGoodTrialsStruct_(session)
    tf = false;
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    if ~isstruct(session.trials), return; end
    req = {'valid','cue','press','lick'};
    tf = all(isfield(session.trials, req));
end

function spikes = getSpikesForSession_(session, S, sessionIdx)
    spikes = [];
    if isfield(session,'spikes') && ~isempty(session.spikes)
        spikes = session.spikes; return
    end
    if isfield(S,'spikes_session') && ~isempty(S.spikes_session) && ...
            numel(S.spikes_session) >= sessionIdx && ~isempty(S.spikes_session{sessionIdx})
        spikes = S.spikes_session{sessionIdx}; return
    end
    if isfield(S,'spikes_persession') && ~isempty(S.spikes_persession) && ...
            numel(S.spikes_persession) >= sessionIdx && ~isempty(S.spikes_persession{sessionIdx})
        spikes = S.spikes_persession{sessionIdx}; return
    end
    if isfield(S,'spikes') && iscell(S.spikes) && ~isempty(S.spikes)
        spikes = S.spikes; % last resort
    end
end