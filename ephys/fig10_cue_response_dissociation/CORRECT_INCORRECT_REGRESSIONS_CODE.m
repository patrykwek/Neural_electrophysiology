%% ========= FIG_NEURAL_BEHAVIOR_COUPLING__MI_VS_PERFORMANCE.m =========
% Neural–behavior coupling: per-session median(|MI|) vs %correct
% Uses the SAME GLOBAL windows and MI definition as your per-session script:
% - SDF built on common cue-aligned axis [winLeft..winRight] (no warping)
% - Epoch windows applied per-trial relative to cue / press / lick (GLOBAL_WINDOWS)
% - Baseline window fixed relative to cue: [-500 -100] ms
% - MI_trial = (FR_epoch - FR_base) / (FR_epoch + FR_base + eps)
% - Per-unit summary: median(|MI_trial|) across trials
% - Per-session summary: median across units
%
% Output: 3 separate panels (one per epoch) with regression line.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE';
outPng  = fullfile(baseOut, 'NEURAL_BEHAVIOR_COUPLING__MI_VS_PERFORMANCE.png');

% >>> REQUIRED: load classification output <<<
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

% Analysis window used to BUILD SDFs
preCueMs   = 500;
postLickMs = 3000;

% Trial selection window (relative to cue)
fixedWin = [-1000, 6000];

% Trial filters
MIN_RT_MS    = 100;
requireValid = true;

% PER-UNIT inclusion
minTrialsPerUnit = 10;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Baseline window for MI (relative to cue)
baselineWin_relCue = [-500, -100];  % ms

% Plot style
figPos  = [120, 120, 1200, 420];
fsTitle = 16;
fsAx    = 13;
lwAx    = 2.0;
lwLine  = 2.5;

% Colors
colCue   = [0.15 0.50 0.80];
colPress = [0.85 0.35 0.15];
colLick  = [0.20 0.65 0.25];

rng(0);

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

%% ---- PRECOMPUTE GAUSS KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- LOAD GLOBAL WINDOWS ----
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

wCue_rel   = GW.cue_ms;
wPress_rel = GW.press_ms;
wLick_rel  = GW.lick_ms;

%% ---- BUILD COMMON CUE-ALIGNED AXIS ----
minLag = min([wCue_rel(1), wPress_rel(1), wLick_rel(1), baselineWin_relCue(1)]);
maxLag = max([wCue_rel(2), wPress_rel(2), wLick_rel(2), baselineWin_relCue(2)]);

winLeft  = fixedWin(1) + minLag;
winRight = fixedWin(2) + maxLag;

winLeft  = min(winLeft, -preCueMs);
winRight = max(winRight, fixedWin(2));

nT = numel(winLeft:dt_ms:winRight);

%% ---- COMPUTE PER-SESSION medAbsMI + perfPct ----
[medAbsMI, perfPct] = computeSessionMedAbsMI_and_Perf_(beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue);

%% ---- PLOT: 3 PANELS SCATTER + REGRESSION ----
fig = figure('Color','w','Position',figPos);

epochNames = {'Cue','Lever press','Reward lick'};
epochCols  = {colCue, colPress, colLick};

for e = 1:3
    ax = subplot(1,3,e); hold(ax,'on');

    x = perfPct(:);
    y = medAbsMI(:,e);

    ok = isfinite(x) & isfinite(y);
    x = x(ok);
    y = y(ok);

    scatter(ax, x, y, 38, 'filled', 'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0);

    % Regression line (simple linear fit)
    if numel(x) >= 2
        p = polyfit(x, y, 1);
        xx = linspace(min(x), max(x), 100);
        yy = polyval(p, xx);
        plot(ax, xx, yy, '-', 'LineWidth', lwLine, 'Color', epochCols{e});

        % Correlation (Spearman, robust)
        rho = corr(x, y, 'Type','Spearman', 'Rows','complete');
        title(ax, sprintf('%s\n\\rho=%.2f', epochNames{e}, rho), ...
            'FontSize', fsTitle, 'FontWeight','bold');
    else
        title(ax, epochNames{e}, 'FontSize', fsTitle, 'FontWeight','bold');
    end

    xlabel(ax, '% correct');
    ylabel(ax, 'Median |MI|');

    set(ax, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    box(ax,'off');
end

saveas(fig, outPng);
fprintf('Saved: %s\n', outPng);

%% ================= LOCAL HELPERS =================

function [medAbsMI, perfPct] = computeSessionMedAbsMI_and_Perf_(beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue)

    nSessions = numel(beh);
    medAbsMI = nan(nSessions,3);
    perfPct  = nan(nSessions,1);

    for sIdx = 1:nSessions
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        % keepTrial: valid/cue/press/lick/RT/fixedWin
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

        % %correct (best-effort)
        perfPct(sIdx) = tryComputePerfPct_(trials(idxKeep));

        unitMedAbsMI = nan(0,3);

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

                % Active trials: >=1 spike within [cue+fixedWin]
                activeTrials = false(numel(idxKeep),1);
                spk = spk_abs;
                if ~issorted(spk), spk = sort(spk); end
                j = 1; nSpk = numel(spk);

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

                MI_cue_tr   = nan(nTrU,1);
                MI_press_tr = nan(nTrU,1);
                MI_lick_tr  = nan(nTrU,1);

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

                    idxBase  = msWindowToIdx_(baselineWin_relCue, winLeft, dt_ms, nT);
                    idxCue   = msWindowToIdx_([0      + wCue_rel(1),   0      + wCue_rel(2)],    winLeft, dt_ms, nT);
                    idxPress = msWindowToIdx_([tPress + wPress_rel(1), tPress + wPress_rel(2)], winLeft, dt_ms, nT);
                    idxLick  = msWindowToIdx_([tLick  + wLick_rel(1),  tLick  + wLick_rel(2)],   winLeft, dt_ms, nT);

                    frBase  = mean(y(idxBase(1):idxBase(2)),   'omitnan');
                    frCue   = mean(y(idxCue(1):idxCue(2)),     'omitnan');
                    frPress = mean(y(idxPress(1):idxPress(2)), 'omitnan');
                    frLick  = mean(y(idxLick(1):idxLick(2)),   'omitnan');

                    MI_cue_tr(iTr)   = (frCue   - frBase) ./ (frCue   + frBase + eps);
                    MI_press_tr(iTr) = (frPress - frBase) ./ (frPress + frBase + eps);
                    MI_lick_tr(iTr)  = (frLick  - frBase) ./ (frLick  + frBase + eps);
                end

                mCue   = median(abs(MI_cue_tr),   'omitnan');
                mPress = median(abs(MI_press_tr), 'omitnan');
                mLick  = median(abs(MI_lick_tr),  'omitnan');

                unitMedAbsMI(end+1,:) = [mCue mPress mLick]; %#ok<AGROW>
            end
        end

        if isempty(unitMedAbsMI)
            continue;
        end

        medAbsMI(sIdx,:) = median(unitMedAbsMI, 1, 'omitnan');
    end
end

function pct = tryComputePerfPct_(trials)
    pct = nan;
    if isempty(trials) || ~isstruct(trials), return; end
    candidates = {'correct','hit','success','rewarded','choice_correct','isCorrect'};
    for i = 1:numel(candidates)
        f = candidates{i};
        if isfield(trials, f)
            v = [trials.(f)];
            if isnumeric(v) || islogical(v)
                v = double(v(:));
                v = v(isfinite(v));
                if ~isempty(v)
                    pct = 100 * mean(v > 0);
                    return;
                end
            end
        end
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