%% ========= FIG3E_EVENT_MODULATION_WITHIN_GROUP__SQUARE__KERNEL_PLUS_HZ__EVENTCOLORS__SESSIONLEVEL.m =========
% PURPOSE:
%   Make TWO Fig3E-style boxplot+points figures (SQUARE):
%     (1) EARLY sessions only: x-axis = event type (Cue/Press/Lick)
%     (2) LATE  sessions only: x-axis = event type (Cue/Press/Lick)
%
%   Do this for BOTH metrics:
%     A) Kernel modulation (session-level; baseline=0):
%          dEvent_k_sess = mean(kernel(lags in event window))
%     B) Hz modulation (session-level; baseline-clean baseline):
%          dEvent_hz_sess = eventHz - baseHz
%
% STYLE:
%   - Same Fig3E boxchart+scatter style (minimal changes)
%   - EVENT COLORS for boxes/points: Cue (red), Press (blue), Lick (green)
%   - SQUARE FIGURES (width == height)
%   - Font sizes + line widths DOUBLED (per request)
%
% STATS:
%   - For each plot (Early or Late), compute within-group paired permutation tests (sign-flip):
%       Cue vs Press, Cue vs Lick, Press vs Lick   (paired by session)
%   - Draw significance bars + stars above the corresponding event pairs
%   - Print observed mean diff + p-value to command window
%
% SAVES:
%   4 figures total: Early(Kernel), Late(Kernel), Early(Hz), Late(Hz)
%     as PNG + SVG + FIG
%   plus a MAT with computed vectors and indices.
%
% NOTE (CHANGED REQUEST):
%   - Additional null plot now uses a NULL ΔHz (circularly shifted spike times),
%     NOT a null kernel refit.
%
% OPTION 2 (CHANGED REQUEST):
%   - Hz modulation is computed per-unit, then averaged across units (per session),
%     instead of pooling spikes across units.

clear; clc;

%% ---- USER SETTINGS ----
outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/GLM';

% GLM outputs used by your Option-A session script
glmSaveFile = fullfile(outDir, 'GLM_SESSIONLEVEL_OUTPUTS_for_window_validity.mat');
if exist(glmSaveFile,'file') ~= 2
    glmSaveFile_alt = fullfile(outDir, 'GLM_SESSIONLEVEL_OUTPUTS_for_window_validity_SIMPLIFIED.mat');
    if exist(glmSaveFile_alt,'file') == 2
        glmSaveFile = glmSaveFile_alt;
    end
end
assert(exist(glmSaveFile,'file')==2, 'Missing saved GLM outputs: %s', glmSaveFile);

% Behavior/spikes file used by your Option-A Hz reconstruction
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
assert(exist(matFile,'file')==2, 'Missing mat file: %s', matFile);

% Global windows (same as Option-A script)
winFile = fullfile(outDir, 'GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat');

% Baseline-clean settings (same as Option-A script)
MIN_ISI_FROM_PREV_LICK_MS = 1000;

% Defaults if not stored
dt_ms = 10;
BASELINE_PRE_CUE_MS = 500;

% Groups (EDIT if needed)
earlyRange = [1 8];
lateRange  = [29 36];

% ---- OUTPUT FILES ----
% Kernel
outPngK_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__KERNEL.png');
outSvgK_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__KERNEL.svg');
outFigK_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__KERNEL.fig');

outPngK_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__KERNEL.png');
outSvgK_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__KERNEL.svg');
outFigK_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__KERNEL.fig');

% Hz
outPngH_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__HZ.png');
outSvgH_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__HZ.svg');
outFigH_E = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__EARLY__HZ.fig');

outPngH_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__HZ.png');
outSvgH_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__HZ.svg');
outFigH_L = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__LATE__HZ.fig');

% MAT
outMat    = fullfile(outDir, 'FIG3E_EVENT_MODULATION_WITHIN_GROUP__KERNEL_PLUS_HZ.mat');

%% ---- PLOT STYLE (ONLY CHANGED: FONT SIZES + LINE WIDTHS DOUBLED) ----
figPos = [120, 120, 720, 720];  % square

fsTitle = 60;     % was 30
fsAx    = 60;     % was 30
lwAx    = 8.0;    % was 4.0

ptSize   = 120;   % was 60
boxLineW = 6.0;   % was 3.0

% stats bar styling
statLineW = 6.0;  % was 3.0

% EVENT COLORS (from your styled session script)
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

rng(0);

%% ---------------- LOAD GLM OUTPUTS ----------------
tmp = load(glmSaveFile);
assert(isfield(tmp,'GLM_OUTPUTS') && isstruct(tmp.GLM_OUTPUTS), 'GLM_OUTPUTS struct missing in: %s', glmSaveFile);
G = tmp.GLM_OUTPUTS;

sessKernels_cue   = G.sessKernels_cue;
sessKernels_press = G.sessKernels_press;
sessKernels_lick  = G.sessKernels_lick;

lags_cue_ms   = G.lags_cue_ms(:);
lags_press_ms = G.lags_press_ms(:);
lags_lick_ms  = G.lags_lick_ms(:);

if isfield(G,'dt_ms') && isfinite(G.dt_ms), dt_ms = G.dt_ms; end
if isfield(G,'BASELINE_PRE_CUE_MS') && isfinite(G.BASELINE_PRE_CUE_MS), BASELINE_PRE_CUE_MS = G.BASELINE_PRE_CUE_MS; end

%% ---------------- LOAD GLOBAL WINDOWS (if available) ----------------
useWindowsFromFile = (exist(winFile,'file')==2);
if useWindowsFromFile
    W = load(winFile);
    assert(isfield(W,'GLOBAL_WINDOWS') && isstruct(W.GLOBAL_WINDOWS), 'GLOBAL_WINDOWS missing in: %s', winFile);
    GLOBAL_WINDOWS = W.GLOBAL_WINDOWS;
    wCue_ms   = GLOBAL_WINDOWS.cue_ms;
    wPress_ms = GLOBAL_WINDOWS.press_ms;
    wLick_ms  = GLOBAL_WINDOWS.lick_ms;
else
    wCue_ms   = [0 300];
    wPress_ms = [-100 200];
    wLick_ms  = [-100 200];
end

%% ---------------- KERNEL MODULATION (PER SESSION) ----------------
Kcue   = sessKernels_cue;
Kpress = sessKernels_press;
Klick  = sessKernels_lick;

nSess_glm = size(Kcue,1);

keepCue   = any(isfinite(Kcue),2);
keepPress = any(isfinite(Kpress),2);
keepLick  = any(isfinite(Klick),2);

idxCueWin   = (lags_cue_ms   >= wCue_ms(1))   & (lags_cue_ms   <= wCue_ms(2));
idxPressWin = (lags_press_ms >= wPress_ms(1)) & (lags_press_ms <= wPress_ms(2));
idxLickWin  = (lags_lick_ms  >= wLick_ms(1))  & (lags_lick_ms  <= wLick_ms(2));

assert(any(idxCueWin),   'Cue window has no overlap with cue lags. Window=[%g %g] ms',     wCue_ms(1), wCue_ms(2));
assert(any(idxPressWin), 'Press window has no overlap with press lags. Window=[%g %g] ms', wPress_ms(1), wPress_ms(2));
assert(any(idxLickWin),  'Lick window has no overlap with lick lags. Window=[%g %g] ms',   wLick_ms(1), wLick_ms(2));

dCue_k_sess   = nan(nSess_glm,1);
dPress_k_sess = nan(nSess_glm,1);
dLick_k_sess  = nan(nSess_glm,1);

dCue_k_sess(keepCue)     = mean(Kcue(keepCue, idxCueWin),        2, 'omitnan');
dPress_k_sess(keepPress) = mean(Kpress(keepPress, idxPressWin),  2, 'omitnan');
dLick_k_sess(keepLick)   = mean(Klick(keepLick, idxLickWin),     2, 'omitnan');

%% ---------------- HZ MODULATION (PER SESSION; PER-UNIT THEN AVG; BASELINE-CLEAN) ----------------
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in: %s', matFile);
nSessions = numel(beh);

% IMPORTANT: assume GLM sessions match beh session order
if nSess_glm ~= nSessions
    warning('GLM has %d sessions but beh has %d sessions. Using nS=min(...).', nSess_glm, nSessions);
end
nS = min(nSess_glm, nSessions);

dCue_hz_sess   = nan(nS,1);
dPress_hz_sess = nan(nS,1);
dLick_hz_sess  = nan(nS,1);

for sIdx = 1:nS
    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials) || ~isstruct(session.trials)
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        continue;
    end

    [idxKeep, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask] = ...
        selectTrials_TrialBounds_forHz_(trials, true, BASELINE_PRE_CUE_MS, MIN_ISI_FROM_PREV_LICK_MS);
    if isempty(idxKeep)
        continue;
    end

    % OPTION 2: compute ΔHz per unit, then average across units for this session
    [dCue_hz_sess(sIdx), dPress_hz_sess(sIdx), dLick_hz_sess(sIdx)] = ...
        session_event_deltaHz_baselineClean_unitsAvg_(spikes_by_ch, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask, ...
            wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms);
end

%% ---------------- BUILD EARLY/LATE SESSION INDICES (CLAMPED) ----------------
earlyIdx = clampRange_(earlyRange, nS);
lateIdx  = clampRange_(lateRange,  nS);

% Clamp kernel vectors to nS as well
dCue_k_use   = dCue_k_sess(1:nS);
dPress_k_use = dPress_k_sess(1:nS);
dLick_k_use  = dLick_k_sess(1:nS);

epochNames = { ...
    sprintf('Cue  [%d,%d] ms',   round(wCue_ms(1)),   round(wCue_ms(2))), ...
    sprintf('Press  [%d,%d] ms', round(wPress_ms(1)), round(wPress_ms(2))), ...
    sprintf('Lick  [%d,%d] ms',  round(wLick_ms(1)),  round(wLick_ms(2)))};

%% ========================= PLOT: KERNEL (EARLY + LATE; x=event) =========================
plot_event_boxes_onegroup_with_tstats_( ...
    dCue_k_use, dPress_k_use, dLick_k_use, earlyIdx, epochNames, ...
    'Kernel Early', ...
    'Modulation (kernel)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPngK_E, outSvgK_E, outFigK_E);

plot_event_boxes_onegroup_with_tstats_( ...
    dCue_k_use, dPress_k_use, dLick_k_use, lateIdx, epochNames, ...
    'Kernel Late', ...
    'Modulation (kernel)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPngK_L, outSvgK_L, outFigK_L);

%% ========================= PLOT: HZ (EARLY + LATE; x=event) =========================
plot_event_boxes_onegroup_with_tstats_( ...
    dCue_hz_sess, dPress_hz_sess, dLick_hz_sess, earlyIdx, epochNames, ...
    'Hz Early', ...
    'Modulation (Hz \Delta)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPngH_E, outSvgH_E, outFigH_E);

plot_event_boxes_onegroup_with_tstats_( ...
    dCue_hz_sess, dPress_hz_sess, dLick_hz_sess, lateIdx, epochNames, ...
    'Hz Late', ...
    'Modulation (Hz \Delta)', ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, ...
    colCue, colPress, colLick, ...
    outPngH_L, outSvgH_L, outFigH_L);

%% ---- SAVE MAT ----
save(outMat, ...
    'earlyRange','lateRange','earlyIdx','lateIdx', ...
    'wCue_ms','wPress_ms','wLick_ms','dt_ms','BASELINE_PRE_CUE_MS','MIN_ISI_FROM_PREV_LICK_MS', ...
    'dCue_k_use','dPress_k_use','dLick_k_use', ...
    'dCue_hz_sess','dPress_hz_sess','dLick_hz_sess', ...
    'glmSaveFile','matFile','winFile');

fprintf('Saved MAT: %s\n', outMat);

%% ========================= ADDITIONAL PLOTS (REQUESTED) =========================
% (1) NULL ΔHz plot (circular shift spike times): observed vs null violin with dots (Cue/Press/Lick)
% (2) Covariate plot: regression coefficient (beta_stage ± CI) across epochs (one panel)

% ---- OUTPUT FILES: NULL ΔHz (CIRC SHIFT SPIKES) ----
outPngNullHz = fullfile(outDir, 'CIRCULARSHIFT_NULL__OBS_VS_NULL__HZ__SESSIONLEVEL.png');
outSvgNullHz = fullfile(outDir, 'CIRCULARSHIFT_NULL__OBS_VS_NULL__HZ__SESSIONLEVEL.svg');
outFigNullHz = fullfile(outDir, 'CIRCULARSHIFT_NULL__OBS_VS_NULL__HZ__SESSIONLEVEL.fig');

% ---- OUTPUT FILES: STAGE REGRESSION COEF (KERNEL) ----
outPngBetaK = fullfile(outDir, 'STAGE_BETA__ADJUSTED__KERNEL__SESSIONLEVEL.png');
outSvgBetaK = fullfile(outDir, 'STAGE_BETA__ADJUSTED__KERNEL__SESSIONLEVEL.svg');
outFigBetaK = fullfile(outDir, 'STAGE_BETA__ADJUSTED__KERNEL__SESSIONLEVEL.fig');

% ---- SETTINGS for additional plots ----
FAST_PREVIEW = true;

if FAST_PREVIEW
    nShuffles = 100;
    useSessForNull = unique([earlyIdx(1:min(4,numel(earlyIdx))); lateIdx(1:min(4,numel(lateIdx)))]);
else
    nShuffles = 100;
    useSessForNull = unique([earlyIdx(:); lateIdx(:)]);
end

useEarlyLateForNullDots = true;

% Observed dots for NULL plot: use OBSERVED ΔHz (same units as null)
obsCue_null   = dCue_hz_sess(:);
obsPress_null = dPress_hz_sess(:);
obsLick_null  = dLick_hz_sess(:);

nullCue_all   = [];
nullPress_all = [];
nullLick_all  = [];

fprintf('\n=== Building NULL ΔHz by circularly shifting spike times (per-unit then avg) ===\n');
fprintf('nSessions used for null = %d, nShuffles per session = %d\n', numel(useSessForNull), nShuffles);

for ii = 1:numel(useSessForNull)
    sIdx = useSessForNull(ii);

    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials) || ~isstruct(session.trials)
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        continue;
    end

    [idxKeep, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask] = ...
        selectTrials_TrialBounds_forHz_(trials, true, BASELINE_PRE_CUE_MS, MIN_ISI_FROM_PREV_LICK_MS);
    if isempty(idxKeep)
        continue;
    end

    sessT0 = min(tStartK);
    sessT1 = max(tEndK);
    if ~(isfinite(sessT0) && isfinite(sessT1) && sessT1 > sessT0)
        continue;
    end
    dur = sessT1 - sessT0;
    if ~(isfinite(dur) && dur > 0)
        continue;
    end

    for sh = 1:nShuffles
        % random non-zero circular shift (ms)
        kShift = randi([round(0.1*dur) max(round(0.1*dur)+1, round(0.9*dur))], 1, 1);

        % shift each unit individually (same shift applied to all units)
        spikes_by_ch_sh = circshift_spikes_by_ch_(spikes_by_ch, sessT0, dur, kShift);

        [dC, dP, dL] = session_event_deltaHz_baselineClean_unitsAvg_(spikes_by_ch_sh, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask, ...
            wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms);

        if isfinite(dC), nullCue_all(end+1,1)   = dC; end %#ok<AGROW>
        if isfinite(dP), nullPress_all(end+1,1) = dP; end %#ok<AGROW>
        if isfinite(dL), nullLick_all(end+1,1)  = dL; end %#ok<AGROW>
    end
end

plot_null_violin_obs_dots_( ...
    nullCue_all, nullPress_all, nullLick_all, ...
    obsCue_null, obsPress_null, obsLick_null, earlyIdx, lateIdx, useEarlyLateForNullDots, ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, ...
    colCue, colPress, colLick, ...
    'Null(circshift spikes)', ...
    'Modulation (Hz \Delta)', ...
    outPngNullHz, outSvgNullHz, outFigNullHz);

%% ---- Covariate plot: beta_stage ± CI across epochs ----
stage = nan(nS,1);
stage(earlyIdx) = 0;
stage(lateIdx)  = 1;

RT_press = nan(nS,1);
RT_lick  = nan(nS,1);
succRate = nan(nS,1);
nTrials  = nan(nS,1);

for sIdx = 1:nS
    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials) || ~isstruct(session.trials)
        continue;
    end
    tr = session.trials;
    nTrials(sIdx) = numel(tr);

    cue   = nan(numel(tr),1);
    press = nan(numel(tr),1);
    lick  = nan(numel(tr),1);
    vld   = nan(numel(tr),1);

    for k = 1:numel(tr)
        if isfield(tr(k),'cue'),   cue(k)   = double(tr(k).cue); end
        if isfield(tr(k),'press'), press(k) = double(tr(k).press); end
        if isfield(tr(k),'lick'),  lick(k)  = double(tr(k).lick); end
        if isfield(tr(k),'valid'), vld(k)   = double(tr(k).valid); end
    end

    kp = isfinite(cue) & isfinite(press);
    if any(kp)
        RT_press(sIdx) = median(press(kp) - cue(kp), 'omitnan');
    end
    kl = isfinite(press) & isfinite(lick);
    if any(kl)
        RT_lick(sIdx)  = median(lick(kl) - press(kl), 'omitnan');
    end
    if any(isfinite(vld))
        succRate(sIdx) = mean(vld, 'omitnan');
    end
end

covCue   = RT_press;
covPress = RT_press;
covLick  = RT_lick;

obsCue   = dCue_k_use(:);
obsPress = dPress_k_use(:);
obsLick  = dLick_k_use(:);

if FAST_PREVIEW
    bootB = 50;
else
    bootB = 500;
end

[betaStageCue,   ciCue]   = stage_beta_boot_(obsCue,   stage, covCue,   nTrials, succRate, bootB);
[betaStagePress, ciPress] = stage_beta_boot_(obsPress, stage, covPress, nTrials, succRate, bootB);
[betaStageLick,  ciLick]  = stage_beta_boot_(obsLick,  stage, covLick,  nTrials, succRate, bootB);

betas = [betaStageCue, betaStagePress, betaStageLick];
cis   = [ciCue(:)'; ciPress(:)'; ciLick(:)'];

plot_stage_beta_forest_( ...
    betas, cis, ...
    figPos, fsTitle, fsAx, lwAx, boxLineW, ...
    colCue, colPress, colLick, ...
    '\beta_{stage}', ...
    '\beta_{stage} (kernel)', ...
    outPngBetaK, outSvgBetaK, outFigBetaK);

fprintf('\nDONE.\n');

%% =============================== HELPERS ===============================

function plot_event_boxes_onegroup_with_tstats_(yCue, yPress, yLick, idxGroup, epochNames, ttl, ylab, ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, statLineW, colCue, colPress, colLick, outPng, outSvg, outFig)

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    vCue   = yCue(idxGroup);
    vPress = yPress(idxGroup);
    vLick  = yLick(idxGroup);

    y1 = vCue;   y1 = y1(isfinite(y1));
    y2 = vPress; y2 = y2(isfinite(y2));
    y3 = vLick;  y3 = y3(isfinite(y3));

    bc1 = boxchart(ax, 1*ones(size(y1)), y1, 'BoxFaceColor', colCue,   'MarkerStyle','none');
    bc2 = boxchart(ax, 2*ones(size(y2)), y2, 'BoxFaceColor', colPress, 'MarkerStyle','none');
    bc3 = boxchart(ax, 3*ones(size(y3)), y3, 'BoxFaceColor', colLick,  'MarkerStyle','none');

    bc1.LineWidth = boxLineW;
    bc2.LineWidth = boxLineW;
    bc3.LineWidth = boxLineW;

    jitter = 0.18;
    scatter(ax, 1 + (rand(size(y1))-0.5)*2*jitter, y1, ptSize, 'filled', ...
        'MarkerFaceColor', colCue,   'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 2 + (rand(size(y2))-0.5)*2*jitter, y2, ptSize, 'filled', ...
        'MarkerFaceColor', colPress, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 3 + (rand(size(y3))-0.5)*2*jitter, y3, ptSize, 'filled', ...
        'MarkerFaceColor', colLick,  'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);

    set(ax, 'XLim',[0.4 3.6], 'XTick',[1 2 3], ...
        'XTickLabel',{'Cue','Press','Lick'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out', 'Layer','top');
    xlabel(ax, 'Event', 'FontSize', fsAx);
    ylabel(ax, ylab,    'FontSize', fsAx);
    title(ax, ttl, 'FontSize', fsTitle, 'FontWeight','bold');
    box(ax,'off');
    ax.Position = [0.14 0.16 0.82 0.76];

    % ---- PERMUTATION SETTINGS (minimal, hard-coded) ----
    Nperm = 10000;

    fprintf('\n=== Paired sign-flip permutation tests within group: %s ===\n', ttl);

    yAll = [y1(:); y2(:); y3(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        yMin = 0; yMax = 1;
    else
        yMin = min(yAll);
        yMax = max(yAll);
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

        diffs = va2(:) - vb2(:);
        diffs = diffs(isfinite(diffs));

        Tobs = mean(diffs, 'omitnan');

        rng(0); % deterministic
        Tperm = nan(Nperm,1);
        for ip = 1:Nperm
            flip = (rand(size(diffs)) > 0.5) * 2 - 1; % +/-1
            Tperm(ip) = mean(diffs .* flip, 'omitnan');
        end
        p = (1 + nnz(abs(Tperm) >= abs(Tobs))) / (Nperm + 1);

        stars = '';
        if p < 0.001
            stars = '***';
        elseif p < 0.01
            stars = '**';
        elseif p < 0.05
            stars = '*';
        end

        fprintf('Events %d vs %d: mean(diff) = %.6g, p = %.6g, n = %d (Nperm=%d)\n', ...
            a, b, Tobs, p, numel(diffs), Nperm);

        if ~isempty(stars)
            yBar = baseBarY + (ci-1)*stepY;
            plot(ax, [a a b b], [yBar yBar+hBar yBar+hBar yBar], 'k-', 'LineWidth', statLineW);
            text(ax, mean([a b]), yBar+hBar+0.01*yRange, stars, ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontSize', fsTitle, 'FontWeight','bold');
        end
    end

    saveas(fig, outPng); fprintf('Saved: %s\n', outPng);
    saveas(fig, outSvg); fprintf('Saved: %s\n', outSvg);
    savefig(fig, outFig); fprintf('Saved: %s\n', outFig);
    close(fig);
end

function plot_null_violin_obs_dots_(nullCue, nullPress, nullLick, obsCue, obsPress, obsLick, earlyIdx, lateIdx, showEarlyLate, ...
    figPos, fsTitle, fsAx, lwAx, ptSize, boxLineW, colCue, colPress, colLick, ttl, ylab, outPng, outSvg, outFig)

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    fsTitleSmall = max(18, round(fsTitle*0.55));

    X = [1 2 3];
    nulls = {nullCue(:), nullPress(:), nullLick(:)};
    cols  = {colCue, colPress, colLick};

    width = 0.32;

    allForY = [nullCue(:); nullPress(:); nullLick(:); obsCue(:); obsPress(:); obsLick(:)];
    allForY = allForY(isfinite(allForY));

    for i = 1:3
        v = nulls{i};
        v = v(isfinite(v));
        if numel(v) < 5
            continue;
        end

        vCoreLo = prctile(v, 0.5);
        vCoreHi = prctile(v, 99.5);
        if ~isfinite(vCoreLo) || ~isfinite(vCoreHi) || vCoreHi <= vCoreLo
            vCoreLo = min(v);
            vCoreHi = max(v);
        end
        yiGrid = linspace(vCoreLo, vCoreHi, 240);

        sdv = std(v, 0, 'omitnan');
        if ~isfinite(sdv) || sdv <= 0
            sdv = 1e-6;
        end
        bw = max(1e-6, 0.20 * sdv);

        [f, yi] = ksdensity(v, yiGrid, 'Bandwidth', bw);
        if ~isempty(f) && max(f) > 0
            f = f / max(f) * width;
        end

        px = [X(i)-f, fliplr(X(i)+f)];
        py = [yi,     fliplr(yi)];
        p = patch(ax, px, py, [0.7 0.7 0.7], 'EdgeColor','none');
        p.FaceAlpha = 0.35;

        med = median(v, 'omitnan');
        lo  = prctile(v, 2.5);
        hi  = prctile(v, 97.5);
        plot(ax, [X(i)-0.18 X(i)+0.18], [med med], 'k-', 'LineWidth', boxLineW);
        plot(ax, [X(i) X(i)], [lo hi], 'k-', 'LineWidth', boxLineW/1.5);
    end

    jitter = 0.14;
    if showEarlyLate
        e = earlyIdx(:);
        l = lateIdx(:);

        yE = {obsCue(e), obsPress(e), obsLick(e)};
        yL = {obsCue(l), obsPress(l), obsLick(l)};

        for i = 1:3
            yy = yE{i}; yy = yy(isfinite(yy));
            scatter(ax, X(i)-0.12 + (rand(size(yy))-0.5)*2*jitter, yy, ptSize, 'filled', ...
                'MarkerFaceColor', cols{i}, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
            yy = yL{i}; yy = yy(isfinite(yy));
            scatter(ax, X(i)+0.12 + (rand(size(yy))-0.5)*2*jitter, yy, ptSize, 'filled', ...
                'MarkerFaceColor', cols{i}, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
        end
        text(ax, 0.55, 0.95, 'dots: early (left) / late (right)', 'Units','normalized', ...
            'FontSize', fsAx*0.55, 'Color',[0 0 0]);
    else
        yO = {obsCue(:), obsPress(:), obsLick(:)};
        for i = 1:3
            yy = yO{i}; yy = yy(isfinite(yy));
            scatter(ax, X(i) + (rand(size(yy))-0.5)*2*jitter, yy, ptSize, 'filled', ...
                'MarkerFaceColor', cols{i}, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
        end
    end

    set(ax, 'XLim',[0.4 3.6], 'XTick',[1 2 3], 'XTickLabel',{'Cue','Press','Lick'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out', 'Layer','top');
    xlabel(ax, 'Event', 'FontSize', fsAx);
    ylabel(ax, ylab,    'FontSize', fsAx);
    title(ax, ttl, 'FontSize', fsTitleSmall, 'FontWeight','bold');
    box(ax,'off');
    ax.Position = [0.14 0.20 0.82 0.68];

    if ~isempty(allForY)
        yLo = prctile(allForY, 0.5);
        yHi = prctile(allForY, 99.5);
        if ~isfinite(yLo) || ~isfinite(yHi) || yHi <= yLo
            yLo = min(allForY);
            yHi = max(allForY);
        end
        pad = 0.08 * max(eps, (yHi - yLo));
        ylim(ax, [yLo - pad, yHi + pad]);
    end

    saveas(fig, outPng); fprintf('Saved: %s\n', outPng);
    saveas(fig, outSvg); fprintf('Saved: %s\n', outSvg);
    savefig(fig, outFig); fprintf('Saved: %s\n', outFig);
    close(fig);
end

function plot_stage_beta_forest_(betas, cis, figPos, fsTitle, fsAx, lwAx, lineW, colCue, colPress, colLick, ttl, ylab, outPng, outSvg, outFig)
    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    fsTitleSmall = max(18, round(fsTitle*0.55));

    X = [1 2 3];
    cols = {colCue, colPress, colLick};

    plot(ax, [0.4 3.6], [0 0], 'k-', 'LineWidth', lineW/1.5);

    for i = 1:3
        b  = betas(i);
        lo = cis(i,1);
        hi = cis(i,2);
        if ~(isfinite(b) && isfinite(lo) && isfinite(hi))
            continue;
        end
        plot(ax, [X(i) X(i)], [lo hi], 'k-', 'LineWidth', lineW);
        scatter(ax, X(i), b, 220, 'filled', 'MarkerFaceColor', cols{i}, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeAlpha', 0);
    end

    set(ax, 'XLim',[0.4 3.6], 'XTick',[1 2 3], 'XTickLabel',{'Cue','Press','Lick'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out', 'Layer','top');
    xlabel(ax, 'Event', 'FontSize', fsAx);
    ylabel(ax, ylab,    'FontSize', fsAx);
    title(ax, ttl, 'FontSize', fsTitleSmall, 'FontWeight','bold');
    box(ax,'off');
    ax.Position = [0.14 0.20 0.82 0.68];

    saveas(fig, outPng); fprintf('Saved: %s\n', outPng);
    saveas(fig, outSvg); fprintf('Saved: %s\n', outSvg);
    savefig(fig, outFig); fprintf('Saved: %s\n', outFig);
    close(fig);
end

function [betaStage, ci] = stage_beta_boot_(y, stage, cov1, nTrials, succRate, B)
    if nargin < 6 || isempty(B), B = 500; end

    keep = isfinite(y) & isfinite(stage);
    y = y(keep);
    stage = stage(keep);
    cov1 = cov1(keep);
    nTrials = nTrials(keep);
    succRate = succRate(keep);

    X = [ones(size(y)), stage(:), cov1(:), nTrials(:), succRate(:)];
    good = all(isfinite(X),2) & isfinite(y);
    y = y(good);
    X = X(good,:);

    if size(X,1) < 6
        betaStage = nan;
        ci = [nan; nan];
        return;
    end

    for j = 3:size(X,2)
        mu = mean(X(:,j), 'omitnan');
        sd = std(X(:,j), 0, 'omitnan');
        if ~isfinite(sd) || sd==0, sd = 1; end
        X(:,j) = (X(:,j) - mu) / sd;
    end

    b = X \ y;
    betaStage = b(2);

    b2 = nan(B,1);
    n = size(X,1);
    for k = 1:B
        idx = randi(n, [n 1]);
        bb = X(idx,:) \ y(idx);
        b2(k) = bb(2);
    end
    ci = prctile(b2, [2.5 97.5])';
end

function idx = clampRange_(range12, nMax)
    a = max(1, min(nMax, range12(1)));
    b = max(1, min(nMax, range12(2)));
    if b < a, tmp = a; a = b; b = tmp; end
    idx = (a:b).';
end

function [idxKeep, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask] = ...
    selectTrials_TrialBounds_forHz_(trials, requireValid, BASELINE_PRE_CUE_MS, MIN_ISI_FROM_PREV_LICK_MS)

    keepTrial = false(numel(trials),1);
    cueAbs    = nan(numel(trials),1);
    pressAbs  = nan(numel(trials),1);
    lickAbs   = nan(numel(trials),1);
    tStart    = nan(numel(trials),1);
    tEnd      = nan(numel(trials),1);

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end

        req = {'cue','press','lick','t_start','t_end'};
        if ~all(isfield(tr, req)), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        ts    = double(tr.t_start);
        te    = double(tr.t_end);

        if any(~isfinite([cue press lick ts te])), continue; end
        if te <= ts, continue; end

        if ~(cue>=ts && cue<=te), continue; end
        if ~(press>=ts && press<=te), continue; end
        if ~(lick>=ts && lick<=te), continue; end

        if cue - BASELINE_PRE_CUE_MS < ts, continue; end

        keepTrial(k) = true;
        cueAbs(k)    = cue;
        pressAbs(k)  = press;
        lickAbs(k)   = lick;
        tStart(k)    = ts;
        tEnd(k)      = te;
    end

    idxKeep   = find(keepTrial);
    cueAbsK   = cueAbs(idxKeep);
    pressAbsK = pressAbs(idxKeep);
    lickAbsK  = lickAbs(idxKeep);
    tStartK   = tStart(idxKeep);
    tEndK     = tEnd(idxKeep);

    baseCleanMask = false(numel(idxKeep),1);
    if numel(idxKeep) >= 2
        prevLick = lickAbsK(1:end-1);
        thisCue  = cueAbsK(2:end);
        isiPrevLickToCue = thisCue - prevLick;
        baseCleanMask(2:end) = (isiPrevLickToCue >= MIN_ISI_FROM_PREV_LICK_MS);
    end
end

function [dCue, dPress, dLick] = session_event_deltaHz_baselineClean_(spk_all, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask, ...
                                                                     wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms)
    nTr = numel(cueAbsK);
    if nTr == 0
        dCue = nan; dPress = nan; dLick = nan;
        return;
    end

    baseHz_tr = nan(nTr,1);
    for tr = 1:nTr
        if ~baseCleanMask(tr)
            continue;
        end
        cue = cueAbsK(tr);
        ts  = tStartK(tr);
        te  = tEndK(tr);

        b0 = cue - BASELINE_PRE_CUE_MS;
        b1 = cue;
        if b0 < ts || b1 > te
            continue;
        end
        baseHz_tr(tr) = mean_rate_in_window_(spk_all, b0, b1, dt_ms);
    end
    baseHz = mean(baseHz_tr, 'omitnan');
    if ~isfinite(baseHz)
        dCue = nan; dPress = nan; dLick = nan;
        return;
    end

    cueHz_tr   = nan(nTr,1);
    pressHz_tr = nan(nTr,1);
    lickHz_tr  = nan(nTr,1);

    for tr = 1:nTr
        cue   = cueAbsK(tr);
        press = pressAbsK(tr);
        lick  = lickAbsK(tr);
        ts    = tStartK(tr);
        te    = tEndK(tr);

        if (cue + wCue_ms(1)   < ts) || (cue + wCue_ms(2)   > te)
            cueHz_tr(tr) = nan;
        else
            cueHz_tr(tr) = mean_rate_in_window_(spk_all, cue + wCue_ms(1), cue + wCue_ms(2), dt_ms);
        end

        if (press + wPress_ms(1) < ts) || (press + wPress_ms(2) > te)
            pressHz_tr(tr) = nan;
        else
            pressHz_tr(tr) = mean_rate_in_window_(spk_all, press + wPress_ms(1), press + wPress_ms(2), dt_ms);
        end

        if (lick + wLick_ms(1)  < ts) || (lick + wLick_ms(2)  > te)
            lickHz_tr(tr) = nan;
        else
            lickHz_tr(tr) = mean_rate_in_window_(spk_all, lick + wLick_ms(1), lick + wLick_ms(2), dt_ms);
        end
    end

    dCue   = mean(cueHz_tr,   'omitnan') - baseHz;
    dPress = mean(pressHz_tr, 'omitnan') - baseHz;
    dLick  = mean(lickHz_tr,  'omitnan') - baseHz;
end

function hz = mean_rate_in_window_(spk, t0, t1, dt_ms)
    if ~(isfinite(t0) && isfinite(t1)) || t1 <= t0
        hz = nan;
        return;
    end
    edges = t0:dt_ms:t1;
    if numel(edges) < 2
        hz = nan;
        return;
    end
    s = spk(spk >= t0 & spk <= t1);
    if isempty(s)
        hz = 0;
        return;
    end
    counts = histcounts(s, edges);
    hzBins = counts / (dt_ms/1000);
    hz = mean(hzBins, 'omitnan');
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

% ===================== OPTION 2 HELPERS (HZ per-unit then avg) =====================

function [dCue_sess, dPress_sess, dLick_sess] = session_event_deltaHz_baselineClean_unitsAvg_(spikes_by_ch, cueAbsK, pressAbsK, lickAbsK, ...
    tStartK, tEndK, baseCleanMask, wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms)

    dCue_u   = [];
    dPress_u = [];
    dLick_u  = [];

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end
        for u = 1:numel(uc)
            sabs = uc{u};
            if isempty(sabs) || ~isnumeric(sabs), continue; end
            spk = double(sabs(:));
            spk = spk(isfinite(spk));
            if isempty(spk), continue; end
            spk = sort(spk);

            [dC, dP, dL] = session_event_deltaHz_baselineClean_(spk, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask, ...
                wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms);

            dCue_u(end+1,1)   = dC; %#ok<AGROW>
            dPress_u(end+1,1) = dP; %#ok<AGROW>
            dLick_u(end+1,1)  = dL; %#ok<AGROW>
        end
    end

    dCue_sess   = mean(dCue_u,   'omitnan');
    dPress_sess = mean(dPress_u, 'omitnan');
    dLick_sess  = mean(dLick_u,  'omitnan');
end

function spikes_by_ch_sh = circshift_spikes_by_ch_(spikes_by_ch, sessT0, dur, kShift)
    spikes_by_ch_sh = spikes_by_ch;

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end
        for u = 1:numel(uc)
            sabs = uc{u};
            if isempty(sabs) || ~isnumeric(sabs), continue; end
            spk = double(sabs(:));
            spk = spk(isfinite(spk));
            if isempty(spk)
                spikes_by_ch_sh{ch}{u} = spk;
                continue;
            end
            spk_sh = sessT0 + mod((spk - sessT0) + kShift, dur);
            spikes_by_ch_sh{ch}{u} = spk_sh;
        end
    end
end