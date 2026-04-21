%% ========= FIG3E_ABSFR_EARLY_AND_LATE_WITHIN_GROUP__SQUARE__EVENTCOLORS__TTEST_STARS__DOUBLED_STYLE.m =========
% PURPOSE:
%   Make TWO square Fig3E-style boxplot+points figures using the SAME DATA
%   computed by:
%     FIG3E_ABS_ZSCORE_EARLY_VS_LATE__GLOBALWARP__FIG2TRIALS_EXACT.m
%   (i.e., session-level mean firing rate in Cue/Press/Lick windows; classified units only)
%
% PLOTS (like your last within-group event plot):
%   1) EARLY sessions only: x-axis = event type (Cue/Press/Lick)
%   2) LATE  sessions only: x-axis = event type (Cue/Press/Lick)
%
% STATS (same style as last code):
%   - Within each plot/group, do paired t-tests (paired by session):
%       Cue vs Press, Cue vs Lick, Press vs Lick
%   - Draw significance bars + stars
%   - Print t-stat + p to command window
%
% STYLE:
%   - Square figure
%   - Event colors for boxes/points (Cue red, Press blue, Lick green)
%   - Font sizes and line widths DOUBLED (matching your last request)

clear; clc;

%% ---- USER SETTINGS (from your FR script) ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/3_GENERAL_CHANGE_EARLY_LATE';

% Required: load classification output
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skipped)

% Analysis window used to BUILD SDFs
preCueMs   = 500;
postLickMs = 3000;

% FIG2-like trial selection window (relative to cue)
fixedWin = [-1000, 6000];

% Trial filters
MIN_RT_MS    = 100;
requireValid = true;

% PER-UNIT inclusion
minTrialsPerUnit = 10;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Groups
earlyRange = [1 8];
lateRange  = [29 36];

%% ---- OUTPUT FILES (EARLY and LATE; FR in event windows) ----
outPng_E = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__EARLY.png');
outSvg_E = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__EARLY.svg');
outFig_E = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__EARLY.fig');

outPng_L = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__LATE.png');
outSvg_L = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__LATE.svg');
outFig_L = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__LATE.fig');

outMat   = fullfile(baseOut, 'FIG3E_WITHIN_GROUP_EVENTWINDOW_FR__EARLY_AND_LATE.mat');

%% ---- STYLE (DOUBLED like your last code) ----
figPos = [120, 120, 720, 720];  % square

fsTitle = 60;
fsAx    = 60;
lwAx    = 8.0;

ptSize   = 120;
boxLineW = 6.0;
statLineW = 6.0;

% EVENT COLORS
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

rng(0);

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

%% ---- PRECOMPUTE GAUSS KERNEL (unit area; Hz after dividing by dt) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- LOAD GLOBAL WINDOWS (ROBUST PATH SEARCH; EXACTLY like your FR script) ----
candidates = { ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut,'GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
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

% Windows stored as [start end] ms RELATIVE TO EACH EVENT (lag axis).
wCue_rel   = GW.cue_ms;
wPress_rel = GW.press_ms;
wLick_rel  = GW.lick_ms;

%% ---- BUILD A SINGLE CUE-ALIGNED TIME AXIS THAT GUARANTEES WINDOW COVERAGE ----
minLag = min([wCue_rel(1), wPress_rel(1), wLick_rel(1)]);
maxLag = max([wCue_rel(2), wPress_rel(2), wLick_rel(2)]);

winLeft  = fixedWin(1) + minLag;
winRight = fixedWin(2) + maxLag;

winLeft  = min(winLeft, -preCueMs);
winRight = max(winRight, fixedWin(2));

nT = numel(winLeft:dt_ms:winRight);

%% ---- COMPUTE SESSION-LEVEL METRICS (NO WARP; CLASSIFIED UNITS ONLY) ----
M_early = computeFRMetricsForSessions_(earlyRange, beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, ...
    winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel);

M_late  = computeFRMetricsForSessions_(lateRange, beh, S, cell_type, ...
    fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
    dt_ms, g, ...
    winLeft, winRight, nT, ...
    wCue_rel, wPress_rel, wLick_rel);

fprintf('\nEarly sessions used: %d | Late sessions used: %d\n', size(M_early,1), size(M_late,1));

epochNames = { ...
    sprintf('Cue  [%d,%d] ms',   round(wCue_rel(1)),   round(wCue_rel(2))), ...
    sprintf('Press  [%d,%d] ms', round(wPress_rel(1)), round(wPress_rel(2))), ...
    sprintf('Lick  [%d,%d] ms',  round(wLick_rel(1)),  round(wLick_rel(2)))};

%% ========================= PLOT: EARLY (x=event) =========================
plot_event_boxes_onegroup_with_tstats_( ...
    M_early(:,1), M_early(:,2), M_early(:,3), epochNames, ...
    sprintf('Event-window firing rate (Hz): Early (%d-%d)', earlyRange(1), earlyRange(2)), ...
    'Mean FR (Hz)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPng_E, outSvg_E, outFig_E);

%% ========================= PLOT: LATE (x=event) =========================
plot_event_boxes_onegroup_with_tstats_( ...
    M_late(:,1), M_late(:,2), M_late(:,3), epochNames, ...
    sprintf('Event-window firing rate (Hz): Late (%d-%d)', lateRange(1), lateRange(2)), ...
    'Mean FR (Hz)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPng_L, outSvg_L, outFig_L);

%% ---- SAVE MAT ----
save(outMat, ...
    'M_early','M_late','earlyRange','lateRange', ...
    'wCue_rel','wPress_rel','wLick_rel','fixedWin','dt_ms','gaussSigmaMs', ...
    'winLeft','winRight','minTrialsPerUnit','MIN_RT_MS','requireValid', ...
    'matFile','cellTypeFile','winFile');

fprintf('\nSaved MAT: %s\n', outMat);
fprintf('\nDONE.\n');

%% =============================== HELPERS ===============================

function plot_event_boxes_onegroup_with_tstats_(vCue, vPress, vLick, epochNames, ttl, ylab, ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, colCue, colPress, colLick, outPng, outSvg, outFig)

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    % For plotting (finite only)
    y1 = vCue;   y1 = y1(isfinite(y1));
    y2 = vPress; y2 = y2(isfinite(y2));
    y3 = vLick;  y3 = y3(isfinite(y3));

    % Boxes
    bc1 = boxchart(ax, 1*ones(size(y1)), y1, 'BoxFaceColor', colCue,   'MarkerStyle','none');
    bc2 = boxchart(ax, 2*ones(size(y2)), y2, 'BoxFaceColor', colPress, 'MarkerStyle','none');
    bc3 = boxchart(ax, 3*ones(size(y3)), y3, 'BoxFaceColor', colLick,  'MarkerStyle','none');
    bc1.LineWidth = boxLineW;
    bc2.LineWidth = boxLineW;
    bc3.LineWidth = boxLineW;

    % Points
    jitter = 0.18;
    scatter(ax, 1 + (rand(size(y1))-0.5)*2*jitter, y1, ptSize, 'filled', ...
        'MarkerFaceColor', colCue,   'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 2 + (rand(size(y2))-0.5)*2*jitter, y2, ptSize, 'filled', ...
        'MarkerFaceColor', colPress, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 3 + (rand(size(y3))-0.5)*2*jitter, y3, ptSize, 'filled', ...
        'MarkerFaceColor', colLick,  'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);

    % Axis style
    set(ax, 'XLim',[0.4 3.6], 'XTick',[1 2 3], 'XTickLabel',{'Cue','Press','Lick'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out', 'Layer','top');
    xlabel(ax, 'Event', 'FontSize', fsAx);
    ylabel(ax, ylab,    'FontSize', fsAx);
    title(ax, ttl, 'FontSize', fsTitle, 'FontWeight','bold');
    box(ax,'off');
    ax.Position = [0.14 0.16 0.82 0.76];

    % ---------- PAIRED T-TESTS + STARS (paired by session) ----------
    fprintf('\n=== Paired t-tests within group: %s ===\n', ttl);

    yAll = [y1(:); y2(:); y3(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        yMin = 0; yMax = 1;
    else
        yMin = min(yAll); yMax = max(yAll);
    end
    yRange = yMax - yMin;
    if ~isfinite(yRange) || yRange == 0, yRange = 1; end

    baseBarY = yMax + 0.08*yRange;
    hBar     = 0.03*yRange;
    stepY    = 0.09*yRange;

    comps = [1 2; 1 3; 2 3];
    for ci = 1:size(comps,1)
        a = comps(ci,1);
        b = comps(ci,2);

        switch a
            case 1, va = vCue;
            case 2, va = vPress;
            case 3, va = vLick;
        end
        switch b
            case 1, vb = vCue;
            case 2, vb = vPress;
            case 3, vb = vLick;
        end

        keep = isfinite(va) & isfinite(vb);
        va2 = va(keep);
        vb2 = vb(keep);

        if numel(va2) < 2
            fprintf('Comp %d-%d: not enough paired samples (n=%d)\n', a, b, numel(va2));
            continue;
        end

        [~, p, ~, stats] = ttest(va2, vb2);
        tstat = stats.tstat;

        stars = '';
        if p < 0.001
            stars = '***';
        elseif p < 0.01
            stars = '**';
        elseif p < 0.05
            stars = '*';
        end

        fprintf('%s vs %s: t = %.4f, p = %.6g, n = %d\n', ...
            epochNames{a}, epochNames{b}, tstat, p, numel(va2));

        if ~isempty(stars)
            yBar = baseBarY + (ci-1)*stepY;
            plot(ax, [a a b b], [yBar yBar+hBar yBar+hBar yBar], 'k-', 'LineWidth', statLineW);
            text(ax, mean([a b]), yBar+hBar+0.01*yRange, stars, ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', fsTitle, 'FontWeight','bold');
        end
    end
    % ---------------------------------------------------------------

    saveas(fig, outPng); fprintf('Saved: %s\n', outPng);
    saveas(fig, outSvg); fprintf('Saved: %s\n', outSvg);
    savefig(fig, outFig); fprintf('Saved: %s\n', outFig);
    close(fig);
end

function M = computeFRMetricsForSessions_(sessRange, beh, S, cell_type, ...
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

    M = nan(0,3); % one row per session

    for sIdx = s0:s1
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        % ---- Behavior-only keepTrial: EXACT (valid/cue/press/lick/RT/fixedWin) ----
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
            pressLat(k)  = rt;   % ms from cue
            lickLat(k)   = lr;   % ms from cue
        end

        idxKeep = find(keepTrial);
        if isempty(idxKeep), continue; end

        cueAbsK   = cueAbs(idxKeep);
        pressLatK = pressLat(idxKeep);
        lickLatK  = lickLat(idxKeep);

        % Per-UNIT metrics -> average within session
        M_units = nan(0,3);

        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end

            for u = 1:numel(uc)

                % ONLY CLASSIFIED UNITS (1/2/3)
                ct = safe_get_celltype_(cell_type, sIdx, ch, u);
                if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                    continue;
                end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                % PER-UNIT activeTrials: >=1 spike within [cue+fixedWin]
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

                    % SDF (Hz) on common cue-aligned axis [winLeft..winRight]
                    spk_rel = spk_abs - cue0;
                    spk_rel = spk_rel(spk_rel >= winLeft & spk_rel <= winRight);

                    counts = zeros(1, nT);
                    if ~isempty(spk_rel)
                        jj = round((spk_rel - winLeft)/dt_ms) + 1;
                        jj = jj(jj >= 1 & jj <= nT);
                        counts = accumarray(jj(:), 1, [nT 1], @sum, 0).';
                    end

                    sm = conv(counts, g, 'same');
                    y  = sm / (dt_ms/1000); % Hz

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
        spikes = S.spikes;
    end
end