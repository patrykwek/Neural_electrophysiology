%% ========= DECODING_STAGE__NEURONVECTOR_AND_WARPED_CONTINUOUS__LOSO_WITHIN_GROUP__PERMTESTS.m =========
% PURPOSE
%   KEEP your original NEURON-VECTOR decoding (Cue/Press/Lick classification) as before,
%   PLUS add the WARPED CONTINUOUS decoding on the whole warped trial (-0.5s pre-cue to +3s after lick).
%
%   (A) NEURON-VECTOR (paper-style FR-vector decoding):
%       - Labels: 1=cue, 2=press, 3=lick (event identity)
%       - Features: per-unit dFR in global event windows (same as before)
%       - Within-group LOSO: train/test within Early only, and within Late only
%       - Outputs: confusion heatmaps + accuracy bar + per-class recall bars
%       - NEW: add permutation test to EVERY barplot (stars drawn from permutation p)
%
%   (B) WARPED CONTINUOUS stage decoding (continuous prediction over warped time):
%       - Warped template: [-0.5s, Tlick_global+3s] in ms, dt=10ms
%       - Labels by warped-time segment:
%           class 1 = cue->press, class 2 = press->lick, class 3 = post-lick
%       - Features: per-unit dFR in sliding windows over warped time
%       - Within-group LOSO, Early vs Late curves
%       - Outputs: accuracy curve + mean-acc bar + peak-acc bar
%       - Bars use permutation tests too
%
% ONE CONTROL (kept):
%   - Drop EARLY units to match mean LATE unit count (VECTOR ONLY), as in your prior code.
%
% NOTE
%   - This script assumes your beh struct has trials with cue/press/lick (and optionally valid/correct/t_start/t_end)
%   - Warping targets are computed GLOBAL across ALL sessions together (median press latency and lick latency)

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);

% global windows file (same robust search pattern you used)
baseOut_forWins = '/Volumes/WD_BLACK/A_THESIS_FINAL';
candidates = { ...
    fullfile(baseOut_forWins,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'POPULATION_CHANGE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'DMS_SIMILARITY_HEATMAPS','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
    };
winFile = '';
for i = 1:numel(candidates)
    if exist(candidates{i},'file') == 2
        winFile = candidates{i};
        break;
    end
end

% output dir
outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/DECODING_STAGE_VECTOR_PLUS_WARPED_CONTINUOUS_WITHIN_GROUP';
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- OUTPUTS: NEURON-VECTOR (event identity decoding) ----
outConfPng_V = fullfile(outDir, 'DECODING_EVENT__CONFUSION__NEURONVECTOR.png');
outConfSvg_V = fullfile(outDir, 'DECODING_EVENT__CONFUSION__NEURONVECTOR.svg');
outConfFig_V = fullfile(outDir, 'DECODING_EVENT__CONFUSION__NEURONVECTOR.fig');

outBarPng_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__NEURONVECTOR.png');
outBarSvg_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__NEURONVECTOR.svg');
outBarFig_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__NEURONVECTOR.fig');

outCueBarPng_V   = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__CueRecall__NEURONVECTOR.png');
outCueBarSvg_V   = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__CueRecall__NEURONVECTOR.svg');
outCueBarFig_V   = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__CueRecall__NEURONVECTOR.fig');

outPressBarPng_V = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__PressRecall__NEURONVECTOR.png');
outPressBarSvg_V = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__PressRecall__NEURONVECTOR.svg');
outPressBarFig_V = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__PressRecall__NEURONVECTOR.fig');

outLickBarPng_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__LickRecall__NEURONVECTOR.png');
outLickBarSvg_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__LickRecall__NEURONVECTOR.svg');
outLickBarFig_V  = fullfile(outDir, 'DECODING_EVENT__EARLY_LATE_BAR__LickRecall__NEURONVECTOR.fig');

%% ---- OUTPUTS: WARPED CONTINUOUS ----
outCurvePng = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__ACCURACY_CURVE.png');
outCurveSvg = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__ACCURACY_CURVE.svg');
outCurveFig = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__ACCURACY_CURVE.fig');

outBarMeanPng = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__MEANACC.png');
outBarMeanSvg = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__MEANACC.svg');
outBarMeanFig = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__MEANACC.fig');

outBarPeakPng = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__PEAKACC.png');
outBarPeakSvg = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__PEAKACC.svg');
outBarPeakFig = fullfile(outDir, 'DECODING_WARPED_CONTINUOUS__BAR__PEAKACC.fig');

%% ---- GROUPS ----
earlySess = 1:8;
lateSess  = 29:36;

%% ---- RELAXED FILTERS (kept consistent) ----
requireValid   = false;
requireCorrect = false;
MIN_RT_MS      = 0;

% baseline relative to cue (ms) for dFR
BASELINE_WIN_MS = [-500 0];

% minimum sample counts to include session in event-identity decoding (vector)
minSamplesPerClass = 5;
minTotalSamples    = 30;
requireAll3InTest   = false;

rng(0);

%% ---- GENERAL SETTINGS ----
DT_MS = 10;
DO_ZSCORE = true;

%% ---- WARPED TEMPLATE SETTINGS ----
preCueMs   = 500;
postLickMs = 3000;
winLeft_ms = -preCueMs;

SW_WIN_MS  = 200;
SW_STEP_MS = 100;

%% ---- PLOT STYLE (match your decoding scripts) ----
figPos = [120, 120, 1350, 520];

fsTitle = 22;
fsAx    = 16;
lwAx    = 2.5;

fsConfLab  = 22;
fsConfTick = 20;
fsConfNum  = 18;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;
axesTickLW    = 4.0;

scatterMS     = 8;
scatterEdgeLW = 1.0;
jit = 0.10;
greyDot = [0.6 0.6 0.6];

colEarlyShade = [1 1 0]; % yellow
colLateShade  = [0 1 0]; % green

% Event colors
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

% Curve colors (continuous)
colEarlyLine = [0.00 0.45 0.74];
colLateLine  = [0.85 0.33 0.10];
curveLW = 5;

% Colormap for confusion plots
cmap = plasma(256);

% Permutation tests for ALL barplots
Nperm_bar = 10;

%% ---------------- LOAD BEH ----------------
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in MAT.');
nSessions = numel(beh);

earlySess = earlySess(earlySess>=1 & earlySess<=nSessions);
lateSess  = lateSess(lateSess>=1  & lateSess<=nSessions);

%% ---------------- PREPASS: DEFINE FIXED NEURON-VECTOR DIMENSION ----------------
[maxCh, maxUperCh, unitCountPerSession] = inferMaxChanUnitsAndCounts_(beh, S);
cumU = [0; cumsum(maxUperCh(:))];
vecDim = sum(maxUperCh);

fprintf('\n=== NEURON-VECTOR DIMENSION ===\n');
fprintf('maxCh=%d, vecDim=%d\n', maxCh, vecDim);

%% ---------------- ONE CHANGE: DROP EARLY UNITS TO MATCH LATE COUNTS ----------------
lateCounts = unitCountPerSession(lateSess);
lateCounts = lateCounts(isfinite(lateCounts) & lateCounts>0);

if isempty(lateCounts)
    targetLateUnits = round(nanmean(unitCountPerSession(unitCountPerSession>0)));
else
    targetLateUnits = round(mean(lateCounts));
end
fprintf('\n=== UNIT-DROPPING CONTROL (EARLY -> LATE MATCH) ===\n');
fprintf('Target #units for EARLY sessions (mean late) = %d\n', targetLateUnits);

keepMaskVec = cell(nSessions,1);
for sIdx = 1:nSessions
    session = beh(sIdx);
    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        keepMaskVec{sIdx} = [];
        continue;
    end

    presentIdx = [];
    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end
        for u = 1:numel(uc)
            sabs = uc{u};
            if isempty(sabs) || ~isnumeric(sabs), continue; end
            idx = unitIndex_(ch, u, maxUperCh, cumU);
            if idx > 0
                presentIdx(end+1) = idx; %#ok<AGROW>
            end
        end
    end
    presentIdx = unique(presentIdx);

    if isempty(presentIdx)
        keepMaskVec{sIdx} = false(1, vecDim);
        continue;
    end

    if ismember(sIdx, earlySess)
        nKeep = min(targetLateUnits, numel(presentIdx));
        keepIdx = presentIdx(randperm(numel(presentIdx), nKeep));
    else
        keepIdx = presentIdx;
    end

    m = false(1, vecDim);
    m(keepIdx) = true;
    keepMaskVec{sIdx} = m;
end

%% ---------------- LOAD WINDOWS (for event-identity decoding) ----------------
if ~isempty(winFile)
    W = load(winFile);
    assert(isfield(W,'GLOBAL_WINDOWS') && isstruct(W.GLOBAL_WINDOWS), 'GLOBAL_WINDOWS missing in: %s', winFile);
    GW = W.GLOBAL_WINDOWS;
    WIN_CUE_MS   = GW.cue_ms;
    WIN_PRESS_MS = GW.press_ms;
    WIN_LICK_MS  = GW.lick_ms;
    fprintf('\nLoaded GLOBAL windows from:\n  %s\n', winFile);
else
    WIN_CUE_MS   = [0 300];
    WIN_PRESS_MS = [-100 200];
    WIN_LICK_MS  = [-100 200];
    fprintf('\nGLOBAL windows file not found. Using defaults.\n');
end

fprintf('\n=== WINDOWS (RELATIVE TO EVENT) ===\n');
fprintf('Cue   [%d, %d] ms\n',   round(WIN_CUE_MS(1)),   round(WIN_CUE_MS(2)));
fprintf('Press [%d, %d] ms\n', round(WIN_PRESS_MS(1)), round(WIN_PRESS_MS(2)));
fprintf('Lick  [%d, %d] ms\n',  round(WIN_LICK_MS(1)),  round(WIN_LICK_MS(2)));

%% ============================================================
%  PART A: EVENT-IDENTITY DECODING (NEURON-VECTOR)  [KEEP AS BEFORE]
%% ============================================================
Xv = cell(nSessions,1);
Yv = cell(nSessions,1);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials) || ~isstruct(session.trials)
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        continue;
    end

    keepMask = keepMaskVec{sIdx};
    if isempty(keepMask), keepMask = false(1, vecDim); end

    xsv = []; ysv = [];

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end
        if requireCorrect
            if ~isfield(tr,'correct') || ~isfinite(double(tr.correct)) || double(tr.correct) ~= 1
                continue;
            end
        end

        req = {'cue','press','lick'};
        if ~all(isfield(tr, req)), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end

        % optional bounds (for clipping only)
        tStart = -inf; tEnd = inf;
        if isfield(tr,'t_start') && isfinite(double(tr.t_start)), tStart = double(tr.t_start); end
        if isfield(tr,'t_end')   && isfinite(double(tr.t_end)),   tEnd   = double(tr.t_end);   end

        % baseline (CLIPPED)
        b0 = cue + BASELINE_WIN_MS(1);
        b1 = cue + BASELINE_WIN_MS(2);
        [b0c, b1c] = clipWin_(b0, b1, tStart, tEnd);
        if ~isfinite(b0c) || ~isfinite(b1c) || b1c <= b0c
            continue;
        end

        % baseline per unit (vector)
        baseHz_vec = zeros(1, vecDim);
        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end
            for u = 1:numel(uc)
                idx = unitIndex_(ch, u, maxUperCh, cumU);
                if idx < 1 || ~keepMask(idx), continue; end
                sabs = uc{u};
                if isempty(sabs) || ~isnumeric(sabs), sabs = []; end
                baseHz_vec(idx) = rateHz_(double(sabs(:)), b0c, b1c);
            end
        end

        % cue window
        [w0c, w1c2] = clipWin_(cue + WIN_CUE_MS(1), cue + WIN_CUE_MS(2), tStart, tEnd);
        if isfinite(w0c) && isfinite(w1c2) && (w1c2 > w0c)
            cueHz_vec = zeros(1, vecDim);
            for ch = 1:numel(spikes_by_ch)
                uc = spikes_by_ch{ch};
                if ~iscell(uc) || isempty(uc), continue; end
                for u = 1:numel(uc)
                    idx = unitIndex_(ch, u, maxUperCh, cumU);
                    if idx < 1 || ~keepMask(idx), continue; end
                    sabs = uc{u};
                    if isempty(sabs) || ~isnumeric(sabs), sabs = []; end
                    cueHz_vec(idx) = rateHz_(double(sabs(:)), w0c, w1c2);
                end
            end
            dCue_vec = cueHz_vec - baseHz_vec;
            if all(isfinite(dCue_vec))
                xsv(end+1,:) = dCue_vec; %#ok<AGROW>
                ysv(end+1,1) = 1; %#ok<AGROW>
            end
        end

        % press window
        [w0p, w1p] = clipWin_(press + WIN_PRESS_MS(1), press + WIN_PRESS_MS(2), tStart, tEnd);
        if isfinite(w0p) && isfinite(w1p) && (w1p > w0p)
            pressHz_vec = zeros(1, vecDim);
            for ch = 1:numel(spikes_by_ch)
                uc = spikes_by_ch{ch};
                if ~iscell(uc) || isempty(uc), continue; end
                for u = 1:numel(uc)
                    idx = unitIndex_(ch, u, maxUperCh, cumU);
                    if idx < 1 || ~keepMask(idx), continue; end
                    sabs = uc{u};
                    if isempty(sabs) || ~isnumeric(sabs), sabs = []; end
                    pressHz_vec(idx) = rateHz_(double(sabs(:)), w0p, w1p);
                end
            end
            dPress_vec = pressHz_vec - baseHz_vec;
            if all(isfinite(dPress_vec))
                xsv(end+1,:) = dPress_vec; %#ok<AGROW>
                ysv(end+1,1) = 2; %#ok<AGROW>
            end
        end

        % lick window
        [w0l, w1l] = clipWin_(lick + WIN_LICK_MS(1), lick + WIN_LICK_MS(2), tStart, tEnd);
        if isfinite(w0l) && isfinite(w1l) && (w1l > w0l)
            lickHz_vec = zeros(1, vecDim);
            for ch = 1:numel(spikes_by_ch)
                uc = spikes_by_ch{ch};
                if ~iscell(uc) || isempty(uc), continue; end
                for u = 1:numel(uc)
                    idx = unitIndex_(ch, u, maxUperCh, cumU);
                    if idx < 1 || ~keepMask(idx), continue; end
                    sabs = uc{u};
                    if isempty(sabs) || ~isnumeric(sabs), sabs = []; end
                    lickHz_vec(idx) = rateHz_(double(sabs(:)), w0l, w1l);
                end
            end
            dLick_vec = lickHz_vec - baseHz_vec;
            if all(isfinite(dLick_vec))
                xsv(end+1,:) = dLick_vec; %#ok<AGROW>
                ysv(end+1,1) = 3; %#ok<AGROW>
            end
        end
    end

    % relaxed inclusion checks (same logic as your original)
    if isempty(xsv) || isempty(ysv), continue; end
    n1 = nnz(ysv==1); n2 = nnz(ysv==2); n3 = nnz(ysv==3);

    if size(xsv,1) < minTotalSamples
        continue;
    end

    if requireAll3InTest
        if any([n1 n2 n3] < minSamplesPerClass), continue; end
    else
        if max([n1 n2 n3]) < minSamplesPerClass, continue; end
    end

    Xv{sIdx} = xsv;
    Yv{sIdx} = ysv;

    fprintf('Session %02d: VECTOR n=%d feat=%d [Cue=%d Press=%d Lick=%d]\n', ...
        sIdx, size(xsv,1), size(xsv,2), n1, n2, n3);
end

% usable sessions for vector decoding
earlyUseV = earlySess(~cellfun(@isempty, Xv(earlySess)));
lateUseV  = lateSess( ~cellfun(@isempty, Xv(lateSess)));

fprintf('\nUsable EARLY sessions (VECTOR): %s\n', mat2str(earlyUseV(:)'));
fprintf('Usable LATE  sessions (VECTOR): %s\n', mat2str(lateUseV(:)'));

%% ---- LOSO within-group (VECTOR) ----
sessAcc_V = nan(nSessions,1);
sessRecall_V = nan(nSessions,3);

fprintf('\n=== [NEURON-VECTOR EVENT] EARLY-ONLY LOSO ===\n');
for sHold = earlyUseV(:)'
    [acc, ~, rec] = losoOneHoldout_(Xv, Yv, earlyUseV, sHold, DO_ZSCORE);
    sessAcc_V(sHold) = acc;
    sessRecall_V(sHold,:) = rec;
    fprintf('Holdout EARLY %02d: acc=%.3f | recall=[%.2f %.2f %.2f]\n', sHold, acc, rec);
end

fprintf('\n=== [NEURON-VECTOR EVENT] LATE-ONLY LOSO ===\n');
for sHold = lateUseV(:)'
    [acc, ~, rec] = losoOneHoldout_(Xv, Yv, lateUseV, sHold, DO_ZSCORE);
    sessAcc_V(sHold) = acc;
    sessRecall_V(sHold,:) = rec;
    fprintf('Holdout LATE  %02d: acc=%.3f | recall=[%.2f %.2f %.2f]\n', sHold, acc, rec);
end

earlyAcc_V = sessAcc_V(earlySess); lateAcc_V = sessAcc_V(lateSess);
earlyAcc_V = earlyAcc_V(isfinite(earlyAcc_V));
lateAcc_V  = lateAcc_V(isfinite(lateAcc_V));

confEarly_V = withinGroupConfusion_(Xv, Yv, earlyUseV, [1 2 3], DO_ZSCORE);
confLate_V  = withinGroupConfusion_(Xv, Yv, lateUseV,  [1 2 3], DO_ZSCORE);

% confusion plot (vector)
figConfV = figure('Color','w','Position',figPos);
subplot(1,2,1);
plotConf_(confEarly_V, sprintf('Early (%d-%d)', earlySess(1), earlySess(end)), fsTitle, fsAx, lwAx, cmap, fsConfLab, fsConfTick, fsConfNum);
subplot(1,2,2);
plotConf_(confLate_V,  sprintf('Late (%d-%d)',  lateSess(1),  lateSess(end)),  fsTitle, fsAx, lwAx, cmap, fsConfLab, fsConfTick, fsConfNum);
sgtitle('Event decoding confusion (neuron vector; early unit-dropped)', 'FontSize', fsTitle, 'FontWeight','bold');

saveas(figConfV, outConfPng_V); savefig(figConfV, outConfFig_V);
try, print(figConfV, outConfSvg_V, '-dsvg'); catch ME, warning('SVG save failed (vector confusion): %s', ME.message); end
close(figConfV);

% accuracy bar (vector) + permutation test
figBarV = plotGroupBar_perm_(earlyAcc_V, lateAcc_V, colEarlyShade, colLateShade, ...
    'Accuracy', 'Accuracy (event identity; neuron vector; early unit-dropped)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);

saveas(figBarV, outBarPng_V); savefig(figBarV, outBarFig_V);
try, print(figBarV, outBarSvg_V, '-dsvg'); catch ME, warning('SVG save failed (vector accuracy): %s', ME.message); end
close(figBarV);

% per-class recall bars (vector) + permutation tests
labels = {'Cue','Press','Lick'};
titlesLower = {'cue','press','lick'};
outPngs_V = {outCueBarPng_V, outPressBarPng_V, outLickBarPng_V};
outSvgs_V = {outCueBarSvg_V, outPressBarSvg_V, outLickBarSvg_V};
outFigs_V = {outCueBarFig_V, outPressBarFig_V, outLickBarFig_V};

for c = 1:3
    eR = sessRecall_V(earlySess, c); lR = sessRecall_V(lateSess, c);
    eR = eR(isfinite(eR)); lR = lR(isfinite(lR));

    figC = plotGroupBar_perm_(eR, lR, colEarlyShade, colLateShade, ...
        'Recall', sprintf('%s recall (event identity; neuron vector)', titlesLower{c}), ...
        labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
        scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);

    saveas(figC, outPngs_V{c}); savefig(figC, outFigs_V{c});
    try, print(figC, outSvgs_V{c}, '-dsvg'); catch ME, warning('SVG save failed (vector %s recall): %s', labels{c}, ME.message); end
    close(figC);
end

%% ============================================================
%  PART B: WARPED CONTINUOUS DECODING (WHOLE WARPED TRIAL, SLIDING WINDOWS)
%% ============================================================

%% ---- GLOBAL warp targets computed across ALL sessions together ----
all_pressLat = [];
all_lickLat  = [];

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end
        if requireCorrect
            if ~isfield(tr,'correct') || ~isfinite(double(tr.correct)) || double(tr.correct) ~= 1
                continue;
            end
        end
        if ~all(isfield(tr, {'cue','press','lick'})), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end
        lr = lick - cue;
        if lr <= rt, continue; end

        all_pressLat(end+1,1) = rt; %#ok<AGROW>
        all_lickLat(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat) && ~isempty(all_lickLat), 'No trials passed filters for warp targets.');
Tpress_global = median(all_pressLat);
Tlick_global  = median(all_lickLat);

fprintf('\n=== GLOBAL WARP TARGETS (ALL sessions together) ===\n');
fprintf('Tpress_global = %.1f ms\n', Tpress_global);
fprintf('Tlick_global  = %.1f ms\n', Tlick_global);

%% ---- template axis ----
winRight_template_ms = Tlick_global + postLickMs;
tgrid_ms = winLeft_ms:DT_MS:winRight_template_ms;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,             winLeft_ms, DT_MS, nT);
idxPT = timeToIdx_(Tpress_global, winLeft_ms, DT_MS, nT);
idxLT = timeToIdx_(Tlick_global,  winLeft_ms, DT_MS, nT);

% sliding windows on template
winBins  = max(2, round(SW_WIN_MS / DT_MS));
stepBins = max(1, round(SW_STEP_MS / DT_MS));
winStarts = 1:stepBins:(nT - winBins + 1);
nW = numel(winStarts);
winCenters_idx = winStarts + floor(winBins/2);
tCenters_ms = tgrid_ms(winCenters_idx);

% labels by warped-time segment
Ywin = ones(nW,1);
Ywin(tCenters_ms >= Tpress_global & tCenters_ms < Tlick_global) = 2;
Ywin(tCenters_ms >= Tlick_global) = 3;

fprintf('\n=== WARPED TEMPLATE ===\n');
fprintf('Template: [%d, %d] ms, dt=%d => nT=%d\n', winLeft_ms, round(winRight_template_ms), DT_MS, nT);
fprintf('Sliding: win=%dms (%db), step=%dms (%db) => nW=%d\n', SW_WIN_MS, winBins, SW_STEP_MS, stepBins, nW);

%% ---- build warped continuous datasets per session ----
Xw = cell(nSessions,1);
Yw = cell(nSessions,1);
Widx = cell(nSessions,1);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    keepMask = keepMaskVec{sIdx};
    if isempty(keepMask), keepMask = false(1, vecDim); end

    xs = [];
    ys = [];
    ws = [];

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end
        if requireCorrect
            if ~isfield(tr,'correct') || ~isfinite(double(tr.correct)) || double(tr.correct) ~= 1
                continue;
            end
        end

        if ~all(isfield(tr, {'cue','press','lick'})), continue; end
        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end
        lr = lick - cue;
        if lr <= rt, continue; end

        % bounds for clipping only
        tStart = -inf; tEnd = inf;
        if isfield(tr,'t_start') && isfinite(double(tr.t_start)), tStart = double(tr.t_start); end
        if isfield(tr,'t_end')   && isfinite(double(tr.t_end)),   tEnd   = double(tr.t_end);   end

        winRight_trial_ms = lr + postLickMs;
        tgrid_trial_ms = winLeft_ms:DT_MS:winRight_trial_ms;
        nT_trial = numel(tgrid_trial_ms);

        idx0 = timeToIdx_(0,   winLeft_ms, DT_MS, nT_trial);
        idxP = timeToIdx_(rt,  winLeft_ms, DT_MS, nT_trial);
        idxL = timeToIdx_(lr,  winLeft_ms, DT_MS, nT_trial);
        if ~(idx0 < idxP && idxP < idxL), continue; end

        Rtrial = zeros(vecDim, nT_trial, 'single');

        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end
            for u = 1:numel(uc)
                idxU = unitIndex_(ch, u, maxUperCh, cumU);
                if idxU < 1 || ~keepMask(idxU), continue; end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                spk_rel = spk_abs - cue;
                spk_rel = spk_rel(spk_rel >= winLeft_ms & spk_rel <= winRight_trial_ms);
                if isempty(spk_rel), continue; end

                jj = round((spk_rel - winLeft_ms)/DT_MS) + 1;
                jj = jj(jj>=1 & jj<=nT_trial);
                if isempty(jj), continue; end

                cts = accumarray(jj, 1, [nT_trial,1], @sum, 0);
                hz  = (cts' ./ (DT_MS/1000));
                Rtrial(idxU,:) = single(hz);
            end
        end

        % warp each unit to template nT
        Rwarp = zeros(vecDim, nT, 'single');
        kA = max(2, (idxPT - idx0T + 1));
        kB = max(2, (idxLT - idxPT + 1));
        kC = max(1, (nT   - idxLT + 1));

        for idxU = 1:vecDim
            if ~keepMask(idxU), continue; end
            y = double(Rtrial(idxU,:));
            if ~any(y), continue; end

            pre = y(1:idx0);

            segA = y(idx0:idxP);
            wA   = warpSegment_floorceilEqn_(segA, kA);

            segB = y(idxP:idxL);
            wB   = warpSegment_floorceilEqn_(segB, kB);

            segC = y(idxL:end);
            if numel(segC) < kC
                segC_fit = [segC, zeros(1, kC-numel(segC))];
            else
                segC_fit = segC(1:kC);
            end

            yy = [pre, wA(2:end), wB(2:end), segC_fit(2:end)];

            if numel(yy) < nT
                yy = [yy, zeros(1, nT-numel(yy))];
            elseif numel(yy) > nT
                yy = yy(1:nT);
            end

            Rwarp(idxU,:) = single(yy);
        end

        % baseline per unit in warped coordinates [-500,0]
        b0i = timeToIdx_(BASELINE_WIN_MS(1), winLeft_ms, DT_MS, nT);
        b1i = timeToIdx_(BASELINE_WIN_MS(2), winLeft_ms, DT_MS, nT);
        if b1i <= b0i, continue; end
        base_vec = mean(double(Rwarp(:, b0i:b1i)), 2, 'omitnan');

        % sliding windows -> samples
        for w = 1:nW
            i0 = winStarts(w);
            i1 = i0 + winBins - 1;

            feat = mean(double(Rwarp(:, i0:i1)), 2, 'omitnan') - base_vec;
            feat(~isfinite(feat)) = 0;

            xs(end+1,:) = feat(:)'; %#ok<AGROW>
            ys(end+1,1) = Ywin(w);  %#ok<AGROW>
            ws(end+1,1) = w;        %#ok<AGROW>
        end
    end

    if ~isempty(xs)
        Xw{sIdx} = xs;
        Yw{sIdx} = ys;
        Widx{sIdx} = ws;
        fprintf('Session %02d: warped samples=%d feat=%d\n', sIdx, size(xs,1), size(xs,2));
    end
end

earlyUseW = earlySess(~cellfun(@isempty, Xw(earlySess)));
lateUseW  = lateSess( ~cellfun(@isempty, Xw(lateSess)));

fprintf('\nUsable EARLY sessions (WARPED): %s\n', mat2str(earlyUseW(:)'));
fprintf('Usable LATE  sessions (WARPED): %s\n', mat2str(lateUseW(:)'));

assert(~isempty(earlyUseW) && ~isempty(lateUseW), 'No usable sessions for warped decoding.');

%% ---- LOSO within-group (WARPED CONTINUOUS) ----
accCurveEarly = nan(numel(earlyUseW), nW);
accCurveLate  = nan(numel(lateUseW),  nW);

sessMeanAccEarly = nan(numel(earlyUseW),1);
sessMeanAccLate  = nan(numel(lateUseW),1);

sessPeakAccEarly = nan(numel(earlyUseW),1);
sessPeakAccLate  = nan(numel(lateUseW),1);

fprintf('\n=== [WARPED CONTINUOUS] EARLY-ONLY LOSO ===\n');
for iH = 1:numel(earlyUseW)
    sHold = earlyUseW(iH);
    [accW, meanAcc, peakAcc] = losoWarpedContinuous_(Xw, Yw, Widx, earlyUseW, sHold, DO_ZSCORE, nW);
    accCurveEarly(iH,:) = accW;
    sessMeanAccEarly(iH) = meanAcc;
    sessPeakAccEarly(iH) = peakAcc;
    fprintf('Holdout EARLY %02d: mean=%.3f peak=%.3f\n', sHold, meanAcc, peakAcc);
end

fprintf('\n=== [WARPED CONTINUOUS] LATE-ONLY LOSO ===\n');
for iH = 1:numel(lateUseW)
    sHold = lateUseW(iH);
    [accW, meanAcc, peakAcc] = losoWarpedContinuous_(Xw, Yw, Widx, lateUseW, sHold, DO_ZSCORE, nW);
    accCurveLate(iH,:) = accW;
    sessMeanAccLate(iH) = meanAcc;
    sessPeakAccLate(iH) = peakAcc;
    fprintf('Holdout LATE  %02d: mean=%.3f peak=%.3f\n', sHold, meanAcc, peakAcc);
end

mEarly = mean(accCurveEarly, 1, 'omitnan');
mLate  = mean(accCurveLate,  1, 'omitnan');

seEarly = std(accCurveEarly, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(accCurveEarly),1)));
seLate  = std(accCurveLate,  0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(accCurveLate),1)));

%% ---- plot continuous curve (one continuous axis) ----
xPlot_s = (tCenters_ms/1000);

figCurve = figure('Color','w','Position',figPos);
ax = axes(figCurve); hold(ax,'on');

patch(ax, [xPlot_s, fliplr(xPlot_s)], [mEarly-seEarly, fliplr(mEarly+seEarly)], ...
    colEarlyLine, 'FaceAlpha', 0.12, 'EdgeColor','none');
patch(ax, [xPlot_s, fliplr(xPlot_s)], [mLate-seLate, fliplr(mLate+seLate)], ...
    colLateLine, 'FaceAlpha', 0.12, 'EdgeColor','none');

plot(ax, xPlot_s, mEarly, '-', 'Color', colEarlyLine, 'LineWidth', curveLW);
plot(ax, xPlot_s, mLate,  '-', 'Color', colLateLine,  'LineWidth', curveLW);

xline(ax, 0, '--', 'Color', colCue,   'LineWidth', axesTickLW);
xline(ax, Tpress_global/1000, '--', 'Color', colPress, 'LineWidth', axesTickLW);
xline(ax, Tlick_global/1000,  '--', 'Color', colLick,  'LineWidth', axesTickLW);

xlabel(ax, 'Warped time from cue (s)', 'FontSize', labelFontSize);
ylabel(ax, 'Accuracy', 'FontSize', labelFontSize);
title(ax, sprintf('Warped continuous stage decoding (win=%dms, step=%dms)', SW_WIN_MS, SW_STEP_MS), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

legend(ax, {'Early','Late'}, 'Location','northeast', 'FontSize', tickFontSize);

set(ax, 'FontSize', tickFontSize, 'TickDir','out', 'LineWidth', axesTickLW, 'Box','off', 'Layer','top');
xlim(ax, [winLeft_ms/1000, winRight_template_ms/1000]);

saveas(figCurve, outCurvePng); savefig(figCurve, outCurveFig);
try, print(figCurve, outCurveSvg, '-dsvg'); catch ME, warning('SVG save failed (curve): %s', ME.message); end
close(figCurve);

%% ---- barplots (mean/peak) with permutation tests ----
figMean = plotGroupBar_perm_(sessMeanAccEarly, sessMeanAccLate, colEarlyShade, colLateShade, ...
    'Mean accuracy', 'Mean accuracy (warped continuous)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figMean, outBarMeanPng); savefig(figMean, outBarMeanFig);
try, print(figMean, outBarMeanSvg, '-dsvg'); catch ME, warning('SVG save failed (bar mean): %s', ME.message); end
close(figMean);

figPeak = plotGroupBar_perm_(sessPeakAccEarly, sessPeakAccLate, colEarlyShade, colLateShade, ...
    'Peak accuracy', 'Peak accuracy (warped continuous)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figPeak, outBarPeakPng); savefig(figPeak, outBarPeakFig);
try, print(figPeak, outBarPeakSvg, '-dsvg'); catch ME, warning('SVG save failed (bar peak): %s', ME.message); end
close(figPeak);

fprintf('\nDONE.\nSaved in: %s\n', outDir);

%% ========================= HELPERS =========================

function [accW, meanAcc, peakAcc] = losoWarpedContinuous_(Xcell, Ycell, WidxCell, sessList, sHold, doZ, nW)
    xTr = []; yTr = [];
    for s = sessList(:)'
        if s == sHold, continue; end
        xTr = [xTr; Xcell{s}]; %#ok<AGROW>
        yTr = [yTr; Ycell{s}]; %#ok<AGROW>
    end
    xTe = Xcell{sHold};
    yTe = Ycell{sHold};
    wTe = WidxCell{sHold};

    accW = nan(1,nW);
    meanAcc = nan;
    peakAcc = nan;

    if isempty(xTr) || isempty(xTe), return; end

    okTr = all(isfinite(xTr),2) & isfinite(yTr);
    okTe = all(isfinite(xTe),2) & isfinite(yTe) & isfinite(wTe);
    xTr = xTr(okTr,:); yTr = yTr(okTr);
    xTe = xTe(okTe,:); yTe = yTe(okTe); wTe = wTe(okTe);

    if size(xTr,1) < 50 || size(xTe,1) < 50, return; end

    if doZ
        [xTr, mu, sig] = zscoreSafe_(xTr);
        xTe = (xTe - mu) ./ sig;
        xTe(~isfinite(xTe)) = 0;
    end

    t = templateLinear('Learner','logistic', 'Lambda',1e-4, 'Regularization','ridge');
    Mdl = fitcecoc(xTr, yTr, 'Learners', t, 'ClassNames', [1 2 3]);
    yHat = predict(Mdl, xTe);

    ok = isfinite(yHat) & isfinite(yTe) & isfinite(wTe);
    yHat = yHat(ok); yTe = yTe(ok); wTe = wTe(ok);

    if numel(yHat) < 50, return; end

    for w = 1:nW
        idx = (wTe == w);
        if nnz(idx) < 10
            accW(w) = nan;
        else
            accW(w) = mean(yHat(idx) == yTe(idx));
        end
    end

    meanAcc = mean(accW, 'omitnan');
    peakAcc = max(accW, [], 'omitnan');
end

function [acc, nTe, rec] = losoOneHoldout_(Xcell, Ycell, sessList, sHold, doZ)
    xTr = []; yTr = [];
    for s = sessList(:)'
        if s == sHold, continue; end
        xTr = [xTr; Xcell{s}]; %#ok<AGROW>
        yTr = [yTr; Ycell{s}]; %#ok<AGROW>
    end
    xTe = Xcell{sHold};
    yTe = Ycell{sHold};

    acc = nan; nTe = 0; rec = [nan nan nan];
    if isempty(xTr) || isempty(xTe), return; end

    okTr = all(isfinite(xTr),2) & isfinite(yTr);
    okTe = all(isfinite(xTe),2) & isfinite(yTe);
    xTr = xTr(okTr,:); yTr = yTr(okTr);
    xTe = xTe(okTe,:); yTe = yTe(okTe);

    if size(xTr,1) < 10 || size(xTe,1) < 10, return; end

    if doZ
        [xTr, mu, sig] = zscoreSafe_(xTr);
        xTe = (xTe - mu) ./ sig;
        xTe(~isfinite(xTe)) = 0;
    end

    t = templateLinear('Learner','logistic', 'Lambda',1e-4, 'Regularization','ridge');
    Mdl = fitcecoc(xTr, yTr, 'Learners', t, 'ClassNames', [1 2 3]);

    yHat = predict(Mdl, xTe);

    ok = isfinite(yHat) & isfinite(yTe);
    nTe = nnz(ok);
    if nTe < 10, return; end

    acc = mean(yHat(ok) == yTe(ok));
    rec = perSessionRecall_(yTe(ok), yHat(ok), [1 2 3]);
end

function [Xz, mu, sig] = zscoreSafe_(X)
    mu = mean(X, 1, 'omitnan');
    sig = std(X, 0, 1, 'omitnan');
    sig(sig <= 0 | ~isfinite(sig)) = 1;
    Xz = (X - mu) ./ sig;
    Xz(~isfinite(Xz)) = 0;
end

function Cnorm = withinGroupConfusion_(Xcell, Ycell, sessList, classOrder, doZ)
    Csum = zeros(numel(classOrder), numel(classOrder));
    for sHold = sessList(:)'
        xTr = []; yTr = [];
        for s = sessList(:)'
            if s == sHold, continue; end
            xTr = [xTr; Xcell{s}]; %#ok<AGROW>
            yTr = [yTr; Ycell{s}]; %#ok<AGROW>
        end
        xTe = Xcell{sHold};
        yTe = Ycell{sHold};
        if isempty(xTr) || isempty(xTe), continue; end

        okTr = all(isfinite(xTr),2) & isfinite(yTr);
        okTe = all(isfinite(xTe),2) & isfinite(yTe);
        xTr = xTr(okTr,:); yTr = yTr(okTr);
        xTe = xTe(okTe,:); yTe = yTe(okTe);

        if size(xTr,1) < 10 || size(xTe,1) < 10, continue; end

        if doZ
            [xTr, mu, sig] = zscoreSafe_(xTr);
            xTe = (xTe - mu) ./ sig;
            xTe(~isfinite(xTe)) = 0;
        end

        t = templateLinear('Learner','logistic', 'Lambda',1e-4, 'Regularization','ridge');
        Mdl = fitcecoc(xTr, yTr, 'Learners', t, 'ClassNames', classOrder);
        yHat = predict(Mdl, xTe);

        ok = isfinite(yHat) & isfinite(yTe);
        if nnz(ok) < 10, continue; end

        C = confusionmat(yTe(ok), yHat(ok), 'Order', classOrder);
        Csum = Csum + C;
    end
    Cnorm = Csum ./ max(1, sum(Csum,2));
end

function plotConf_(C, ttl, fsTitle, fsAx, lwAx, cmap, fsConfLab, fsConfTick, fsConfNum)
    imagesc(C);
    axis square;
    colormap(cmap);
    caxis([0 1]);
    cb = colorbar;
    cb.TickDirection = 'out';
    cb.LineWidth = lwAx;

    labels = {'Cue','Press','Lick'};
    xticks(1:3); yticks(1:3);
    xticklabels(labels); yticklabels(labels);
    xtickangle(45);

    set(gca, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    box off;
    title(ttl, 'FontSize', fsTitle, 'FontWeight','bold');

    for i = 1:3
        for j = 1:3
            v = C(i,j);
            if ~isfinite(v), v = 0; end
            if v > 0.5, col = [1 1 1]; else, col = [0 0 0]; end
            text(j, i, sprintf('%.2f', v), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize', fsConfNum, 'FontWeight','bold', 'Color', col);
        end
    end

    xlabel('Predicted', 'FontSize', fsConfLab);
    ylabel('True',      'FontSize', fsConfLab);
    set(gca, 'FontSize', fsConfTick);
end

function recallDiag = perSessionRecall_(yTrue, yPred, classOrder)
    C = confusionmat(yTrue, yPred, 'Order', classOrder);
    Cnorm = C ./ max(1, sum(C,2));
    recallDiag = diag(Cnorm)';
end

function fig = plotGroupBar_perm_(earlyVals, lateVals, colEarly, colLate, ylab, ttl, ...
    labelFS, titleFS, tickFS, axLW, ms, edgeLW, jit, greyDot, Nperm)

    earlyVals = earlyVals(isfinite(earlyVals));
    lateVals  = lateVals(isfinite(lateVals));

    fig = figure('Color','w','Position',[120,120,520,700]);
    ax = axes(fig); hold(ax,'on');

    mE = mean(earlyVals,'omitnan');
    mL = mean(lateVals,'omitnan');
    semE = std(earlyVals,'omitnan') / sqrt(max(1,numel(earlyVals)));
    semL = std(lateVals,'omitnan')  / sqrt(max(1,numel(lateVals)));

    xPos = [1 2];
    barW = 0.60;

    b = bar(ax, xPos, [mE mL], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axLW);
    b.CData(1,:) = colEarly; b.CData(2,:) = colLate;
    b.FaceAlpha = 0.25;

    errorbar(ax, xPos, [mE mL], [semE semL], 'k', 'LineStyle','none', 'LineWidth', axLW, 'CapSize', 18);

    plot(ax, 1 + (rand(size(earlyVals))-0.5)*2*jit, earlyVals, 'o', 'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', edgeLW);
    plot(ax, 2 + (rand(size(lateVals))-0.5)*2*jit,  lateVals,  'o', 'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', edgeLW);

    set(ax, 'XLim',[0.5 2.5], 'XTick',[1 2], 'XTickLabel',{'Early','Late'});
    ylabel(ax, ylab, 'FontSize', labelFS);
    title(ax, ttl, 'FontWeight','bold', 'FontSize', titleFS);
    set(ax, 'FontSize', tickFS, 'Box','off');
    set(ax, 'TickDir','out', 'LineWidth', axLW, 'Layer','top');

    % permutation test (label-shuffle) on Δ(Late-Early)
    if ~isempty(earlyVals) && ~isempty(lateVals)
        dObs = mean(lateVals,'omitnan') - mean(earlyVals,'omitnan');
        allv = [earlyVals(:); lateVals(:)];
        nE = numel(earlyVals);

        rng(0);
        dPerm = nan(Nperm,1);
        for i = 1:Nperm
            perm = allv(randperm(numel(allv)));
            e = perm(1:nE);
            l = perm(nE+1:end);
            dPerm(i) = mean(l,'omitnan') - mean(e,'omitnan');
        end
        p = (1 + nnz(abs(dPerm) >= abs(dObs))) / (Nperm + 1);

        fprintf('\n%s perm-test: Δ(Late-Early)=%.6g, p=%.6g (Nperm=%d)\n', ttl, dObs, p, Nperm);
        addSigIfNeeded_(ax, 1, 2, p, [mE mL], axLW, tickFS);
    end
end

function addSigIfNeeded_(ax, x1, x2, p, means, lw, fs)
    if ~isfinite(p) || p >= 0.05, return; end
    stars = pToStars_(p);

    yl = ylim(ax);
    yMaxData = max(means);
    yRange = yl(2) - yl(1);
    y = yMaxData + 0.08*yRange;
    h = 0.03*yRange;

    if y + h + 0.05*yRange > yl(2)
        ylim(ax, [yl(1), y + h + 0.08*yRange]);
        yl = ylim(ax); yRange = yl(2)-yl(1);
        y = yMaxData + 0.08*yRange;
        h = 0.03*yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01*yRange, stars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', fs, 'FontWeight','bold', 'Color', [0 0 0]);
end

function s = pToStars_(p)
    if p < 0.001, s='***';
    elseif p < 0.01, s='**';
    else, s='*';
    end
end

function hz = rateHz_(spk, t0, t1)
    if ~(isfinite(t0) && isfinite(t1)) || t1 <= t0
        hz = nan; return;
    end
    n = nnz(spk >= t0 & spk <= t1);
    hz = n / ((t1 - t0)/1000);
end

function [t0c, t1c] = clipWin_(t0, t1, tStart, tEnd)
    t0c = t0; t1c = t1;
    if isfinite(tStart), t0c = max(t0c, tStart); end
    if isfinite(tEnd),   t1c = min(t1c, tEnd);   end
end

function idx = timeToIdx_(t_ms, winLeft_ms, dt_ms, nT)
    idx = round((t_ms - winLeft_ms)/dt_ms) + 1;
    idx = max(1, min(nT, idx));
end

function ywarp = warpSegment_floorceilEqn_(y, k)
    y = y(:)';
    n = numel(y);

    if n < 2 || k < 2
        ywarp = zeros(1, k);
        if n >= 1 && k >= 1, ywarp(1) = y(1); end
        return;
    end

    s = (k-1)/(n-1);
    ywarp = zeros(1, k);

    for i = 1:k
        x = 1 + (i-1)/s;
        x0 = floor(x);
        x1 = ceil(x);
        x0 = max(1, min(n, x0));
        x1 = max(1, min(n, x1));

        if x0 == x1
            ywarp(i) = y(x0);
        else
            w1 = x - x0;
            w0 = 1 - w1;
            ywarp(i) = w0*y(x0) + w1*y(x1);
        end
    end
end

function tf = isGoodTrialsStruct_(session)
    tf = false;
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    if ~isstruct(session.trials), return; end
    req = {'cue','press','lick'};
    tf = all(isfield(session.trials, req));
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

function [maxCh, maxUperCh, unitCountPerSession] = inferMaxChanUnitsAndCounts_(beh, S)
    maxCh = 0;
    maxUperCh = zeros(0,1);
    nSessions = numel(beh);
    unitCountPerSession = nan(nSessions,1);

    for sIdx = 1:nSessions
        session = beh(sIdx);
        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        maxCh = max(maxCh, numel(spikes_by_ch));
        if numel(maxUperCh) < maxCh
            maxUperCh(maxCh,1) = 0; %#ok<AGROW>
        end

        nUnitsHere = 0;
        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc), continue; end
            maxUperCh(ch) = max(maxUperCh(ch), numel(uc));

            for u = 1:numel(uc)
                sabs = uc{u};
                if isempty(sabs) || ~isnumeric(sabs), continue; end
                nUnitsHere = nUnitsHere + 1;
            end
        end
        unitCountPerSession(sIdx) = nUnitsHere;
    end

    if isempty(maxUperCh), maxUperCh = 0; end
end

function idx = unitIndex_(ch, u, maxUperCh, cumU)
    if ch < 1 || ch > numel(maxUperCh), idx = -1; return; end
    if u < 1 || u > maxUperCh(ch), idx = -1; return; end
    idx = cumU(ch) + u;
end

function cmap = plasma(m)
    if nargin < 1 || isempty(m)
        m = size(get(gcf,'colormap'),1);
    end

    anchors = [ ...
        0.050383, 0.029803, 0.527975
        0.186213, 0.018803, 0.587228
        0.287076, 0.010855, 0.627295
        0.381047, 0.001814, 0.653068
        0.471457, 0.005678, 0.659897
        0.557243, 0.047331, 0.643443
        0.636008, 0.112092, 0.605205
        0.705673, 0.184019, 0.552295
        0.765690, 0.255477, 0.497050
        0.816363, 0.329727, 0.436033
        0.857763, 0.406008, 0.370668
        0.889436, 0.484898, 0.301467
        0.910513, 0.566949, 0.229412
        0.920049, 0.652764, 0.156346
        0.916242, 0.742065, 0.087714
        0.896091, 0.835793, 0.029491
        0.940015, 0.975158, 0.131326];

    xA = linspace(0,1,size(anchors,1));
    xQ = linspace(0,1,m);

    cmap = zeros(m,3);
    cmap(:,1) = interp1(xA, anchors(:,1), xQ, 'pchip');
    cmap(:,2) = interp1(xA, anchors(:,2), xQ, 'pchip');
    cmap(:,3) = interp1(xA, anchors(:,3), xQ, 'pchip');

    cmap = max(0, min(1, cmap));
end