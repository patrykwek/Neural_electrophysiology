%% ========= FIG3E_EVENT_MODULATION_EARLY_VS_LATE__KERNEL_PLUS_HZ__GLOBALWINDOWS__SESSIONLEVEL.m =========
% PURPOSE (Fig3E-style boxplots, but with MODULATION):
%   Make the same Early vs Late boxplot+points layout as:
%     FIG3E_ABS_ZSCORE_EARLY_VS_LATE__GLOBALWARP__FIG2TRIALS_EXACT.m
%   but using the *session-level event modulation* metrics computed the SAME WAY as:
%     FIG_SESSION_EVENT_MODULATION__KERNEL_PLUS_HZ__STYLED__A_THESIS_FINAL.m
%
% METRICS:
%   A) Kernel modulation (per session, per event):
%        dEvent_k_sess = mean(kernel(lags in event window))   (baseline=0)
%   B) Hz modulation (per session, per event; population pooled):
%        baseHz = mean(baseline-window Hz) using baseline-clean trials only
%        eventHz = mean(event-window Hz) using ALL valid trials
%        dEvent_hz_sess = eventHz - baseHz
%
% PLOT:
%   - 2 figures (Kernel and Hz), each has 1x3 subplots: Cue / Press / Lick
%   - Each subplot: Early vs Late boxchart + jittered points + ranksum p + stars (same style)
%
% NOTE:
%   Assumes GLM session order matches beh session order. If not, align sessions explicitly.

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

% Output files
outPngK = fullfile(outDir, 'FIG3E_EVENT_MODULATION_EARLY_VS_LATE__KERNEL.png');
outSvgK = fullfile(outDir, 'FIG3E_EVENT_MODULATION_EARLY_VS_LATE__KERNEL.svg');
outPngH = fullfile(outDir, 'FIG3E_EVENT_MODULATION_EARLY_VS_LATE__HZ.png');
outSvgH = fullfile(outDir, 'FIG3E_EVENT_MODULATION_EARLY_VS_LATE__HZ.svg');
outMat  = fullfile(outDir, 'FIG3E_EVENT_MODULATION_EARLY_VS_LATE__KERNEL_PLUS_HZ.mat');

%% ---- PLOT STYLE (match the BOXPLOT script style) ----
figPos  = [120, 120, 1050, 520];

% (DOUBLED)
fsTitle = 40;     % was 20
fsAx    = 32;     % was 16
lwAx    = 5.0;    % was 2.5

% Colors (only for points/boxes) -- same as your boxplot script
colEarly = [0.10 0.60 0.15];
colLate  = [0.90 0.70 0.10];

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
    % fallback
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

%% ---------------- HZ MODULATION (PER SESSION; POPULATION POOL; BASELINE-CLEAN) ----------------
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in: %s', matFile);
nSessions = numel(beh);

% IMPORTANT: assume GLM sessions match beh sessions
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

    % concatenate spikes across units for this session (population firing)
    spk_all = [];
    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end
        for u = 1:numel(uc)
            sabs = uc{u};
            if isempty(sabs) || ~isnumeric(sabs), continue; end
            spk_all = [spk_all; double(sabs(:))]; %#ok<AGROW>
        end
    end
    if isempty(spk_all), continue; end
    spk_all = sort(spk_all);

    [dCue_hz_sess(sIdx), dPress_hz_sess(sIdx), dLick_hz_sess(sIdx)] = ...
        session_event_deltaHz_baselineClean_(spk_all, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK, baseCleanMask, ...
            wCue_ms, wPress_ms, wLick_ms, BASELINE_PRE_CUE_MS, dt_ms);
end

%% ---------------- BUILD EARLY/LATE SESSION INDICES (CLAMPED) ----------------
earlyIdx = clampRange_(earlyRange, nS);
lateIdx  = clampRange_(lateRange,  nS);

% Kernel vectors must also be clamped to nS (since we use nS everywhere)
dCue_k_use   = dCue_k_sess(1:nS);
dPress_k_use = dPress_k_sess(1:nS);
dLick_k_use  = dLick_k_sess(1:nS);

%% ========================= PLOT: KERNEL (Fig3E-like) =========================
figK = figure('Color','w','Position',figPos);

epochNamesK = { ...
    sprintf('Cue  [%d, %d] ms',        round(wCue_ms(1)),   round(wCue_ms(2))), ...
    sprintf('Lever press  [%d, %d] ms',round(wPress_ms(1)), round(wPress_ms(2))), ...
    sprintf('Reward lick  [%d, %d] ms',round(wLick_ms(1)),  round(wLick_ms(2)))};

% ---- PERMUTATION SETTINGS (minimal, hard-coded; same style as code 2) ----
Nperm = 10000;

fprintf('\n=== Permutation tests (Early vs Late; KERNEL modulation; session-level) ===\n');

Yk = {dCue_k_use, dPress_k_use, dLick_k_use};

for e = 1:3
    ax = subplot(1,3,e); hold(ax,'on');

    yE = Yk{e}(earlyIdx);
    yL = Yk{e}(lateIdx);

    yE = yE(isfinite(yE));
    yL = yL(isfinite(yL));

    boxchart(ax, ones(size(yE)), yE, 'BoxFaceColor', colEarly, 'MarkerStyle','none');
    boxchart(ax, 2*ones(size(yL)), yL, 'BoxFaceColor', colLate,  'MarkerStyle','none');

    jitter = 0.18;
    scatter(ax, 1 + (rand(size(yE))-0.5)*2*jitter, yE, 28, 'filled', ...
        'MarkerFaceColor', colEarly, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 2 + (rand(size(yL))-0.5)*2*jitter, yL, 28, 'filled', ...
        'MarkerFaceColor', colLate,  'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);

    % ---- UNPAIRED LABEL-PERMUTATION TEST (two-sided; diff-of-means) ----
    if numel(yE) < 2 || numel(yL) < 2
        p = nan;
        Tobs = nan;
        fprintf('%s: not enough samples (nE=%d, nL=%d)\n', epochNamesK{e}, numel(yE), numel(yL));
    else
        Tobs = mean(yE, 'omitnan') - mean(yL, 'omitnan');

        pool = [yE(:); yL(:)];
        nE = numel(yE);

        rng(0); % deterministic
        Tperm = nan(Nperm,1);
        for ip = 1:Nperm
            rp = pool(randperm(numel(pool)));
            Tperm(ip) = mean(rp(1:nE), 'omitnan') - mean(rp(nE+1:end), 'omitnan');
        end
        p = (1 + nnz(abs(Tperm) >= abs(Tobs))) / (Nperm + 1);

        fprintf('%s: mean(diff)=%.6g (Early-Late), p=%.6g, nE=%d, nL=%d (Nperm=%d)\n', ...
            epochNamesK{e}, Tobs, p, numel(yE), numel(yL), Nperm);
    end

    stars = '';
    if isfinite(p)
        if p < 0.001
            stars = '***';
        elseif p < 0.01
            stars = '**';
        elseif p < 0.05
            stars = '*';
        end
    end

    if ~isempty(stars)
        yMax = max([yE(:); yL(:)]);
        yMin = min([yE(:); yL(:)]);
        yRange = yMax - yMin;
        if ~isfinite(yRange) || yRange == 0, yRange = 1; end

        yBar = yMax + 0.08*yRange;
        hBar = 0.03*yRange;

        plot(ax, [1 1 2 2], [yBar yBar+hBar yBar+hBar yBar], 'k-', 'LineWidth', 4); % was 2 (DOUBLED)
        text(ax, 1.5, yBar+hBar+0.01*yRange, stars, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', fsTitle, 'FontWeight','bold');
    end

    title(ax, epochNamesK{e}, 'FontSize', fsTitle, 'FontWeight','bold');
    set(ax, 'XLim',[0.4 2.6], 'XTick',[1 2], 'XTickLabel',{'Early','Late'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    ylabel(ax, 'Response (kernel)');
    box(ax,'off');
end

sgtitle(sprintf('Event modulation (KERNEL): Early (%d-%d) vs Late (%d-%d)', ...
    earlyRange(1), earlyRange(2), lateRange(1), lateRange(2)), ...
    'FontSize', fsTitle+4, 'FontWeight','bold'); % was fsTitle+2 (DOUBLED offset)

saveas(figK, outPngK); fprintf('Saved: %s\n', outPngK);
saveas(figK, outSvgK); fprintf('Saved: %s\n', outSvgK);

%% ========================= PLOT: HZ (Fig3E-like) =========================
figH = figure('Color','w','Position',figPos);

epochNamesH = epochNamesK; % same window labels
fprintf('\n=== Permutation tests (Early vs Late; HZ modulation; session-level) ===\n');

Yh = {dCue_hz_sess, dPress_hz_sess, dLick_hz_sess};

for e = 1:3
    ax = subplot(1,3,e); hold(ax,'on');

    yE = Yh{e}(earlyIdx);
    yL = Yh{e}(lateIdx);

    yE = yE(isfinite(yE));
    yL = yL(isfinite(yL));

    boxchart(ax, ones(size(yE)), yE, 'BoxFaceColor', colEarly, 'MarkerStyle','none');
    boxchart(ax, 2*ones(size(yL)), yL, 'BoxFaceColor', colLate,  'MarkerStyle','none');

    jitter = 0.18;
    scatter(ax, 1 + (rand(size(yE))-0.5)*2*jitter, yE, 28, 'filled', ...
        'MarkerFaceColor', colEarly, 'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);
    scatter(ax, 2 + (rand(size(yL))-0.5)*2*jitter, yL, 28, 'filled', ...
        'MarkerFaceColor', colLate,  'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0);

    % ---- UNPAIRED LABEL-PERMUTATION TEST (two-sided; diff-of-means) ----
    if numel(yE) < 2 || numel(yL) < 2
        p = nan;
        Tobs = nan;
        fprintf('%s: not enough samples (nE=%d, nL=%d)\n', epochNamesH{e}, numel(yE), numel(yL));
    else
        Tobs = mean(yE, 'omitnan') - mean(yL, 'omitnan');

        pool = [yE(:); yL(:)];
        nE = numel(yE);

        rng(0); % deterministic
        Tperm = nan(Nperm,1);
        for ip = 1:Nperm
            rp = pool(randperm(numel(pool)));
            Tperm(ip) = mean(rp(1:nE), 'omitnan') - mean(rp(nE+1:end), 'omitnan');
        end
        p = (1 + nnz(abs(Tperm) >= abs(Tobs))) / (Nperm + 1);

        fprintf('%s: mean(diff)=%.6g (Early-Late), p=%.6g, nE=%d, nL=%d (Nperm=%d)\n', ...
            epochNamesH{e}, Tobs, p, numel(yE), numel(yL), Nperm);
    end

    stars = '';
    if isfinite(p)
        if p < 0.001
            stars = '***';
        elseif p < 0.01
            stars = '**';
        elseif p < 0.05
            stars = '*';
        end
    end

    if ~isempty(stars)
        yMax = max([yE(:); yL(:)]);
        yMin = min([yE(:); yL(:)]);
        yRange = yMax - yMin;
        if ~isfinite(yRange) || yRange == 0, yRange = 1; end

        yBar = yMax + 0.08*yRange;
        hBar = 0.03*yRange;

        plot(ax, [1 1 2 2], [yBar yBar+hBar yBar+hBar yBar], 'k-', 'LineWidth', 4); % was 2 (DOUBLED)
        text(ax, 1.5, yBar+hBar+0.01*yRange, stars, 'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom', 'FontSize', fsTitle, 'FontWeight','bold');
    end

    title(ax, epochNamesH{e}, 'FontSize', fsTitle, 'FontWeight','bold');
    set(ax, 'XLim',[0.4 2.6], 'XTick',[1 2], 'XTickLabel',{'Early','Late'}, ...
        'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    ylabel(ax, 'Response (Hz \Delta)');
    box(ax,'off');
end

sgtitle(sprintf('Event modulation (Hz \x0394): Early (%d-%d) vs Late (%d-%d)', ...
    earlyRange(1), earlyRange(2), lateRange(1), lateRange(2)), ...
    'FontSize', fsTitle+4, 'FontWeight','bold'); % was fsTitle+2 (DOUBLED offset)

saveas(figH, outPngH); fprintf('Saved: %s\n', outPngH);
saveas(figH, outSvgH); fprintf('Saved: %s\n', outSvgH);

%% ---- SAVE MAT ----
save(outMat, ...
    'earlyRange','lateRange','earlyIdx','lateIdx', ...
    'wCue_ms','wPress_ms','wLick_ms','dt_ms','BASELINE_PRE_CUE_MS','MIN_ISI_FROM_PREV_LICK_MS', ...
    'dCue_k_use','dPress_k_use','dLick_k_use', ...
    'dCue_hz_sess','dPress_hz_sess','dLick_hz_sess', ...
    'glmSaveFile','matFile');

fprintf('Saved MAT: %s\n', outMat);
fprintf('\nDONE.\n');

%% =============================== HELPERS ===============================

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

    % ---- BASELINE: ONLY baseline-clean trials ----
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

    % ---- EVENTS: ALL trials ----
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