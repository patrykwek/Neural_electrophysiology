%% ========= FIG_SESSION_EVENT_MODULATION__KERNEL_PLUS_HZ__STYLED__A_THESIS_FINAL.m =========
% PURPOSE:
%   Plot the *session-level event modulation values* computed in:
%     ADDON_OPTIONA_EVENT_MODULATION_DOTPLOT__SESSIONLEVEL_TTEST_STARS_PLUS_HZ.m
%   but rendered in the *exact style* of:
%     FIG_SESSION_EPOCH_FIRINGRATE_HZ__GLOBALWINDOWS_STYLED__A_THESIS_FINAL.m
%
% WHAT IS PLOTTED (same computations as your Option-A script above):
%   1) Kernel metric (per session, per event):
%        deltaEvent_sess = mean(kernel(lags in event window))
%      Baseline corresponds to 0 (kernel already baseline-subtracted by construction).
%
%   2) Hz metric (per session, per event; population spikes pooled):
%        baseHz = mean( baseline-window Hz ) using baseline-clean trials only
%        eventHz = mean( event-window Hz ) using ALL valid trials
%        deltaHz_sess = eventHz - baseHz
%
% PLOT STYLE (copied from your styled per-session FR script):
%   - x-axis: Session
%   - 3 lines: Cue / Press / Lick
%   - thick lines + filled circular markers
%   - manual tick marks (every 5 sessions labeled)
%   - early/late background shading based on sStar
%   - optional transition line at sStar+0.5
%
% SAVES:
%   PNG/SVG/MAT into outDir (GLM folder by default).
%
% NOTE:
%   This script assumes the GLM session-level outputs contain session kernels in
%   session order matching beh struct order. If not, you must align sessions explicitly.

clear; clc;

%% ---------------- USER SETTINGS ----------------
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
matFileDefault = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
assert(exist(matFileDefault,'file')==2, 'Missing mat file: %s', matFileDefault);

% Global windows (same as Option-A script)
winFile = fullfile(outDir, 'GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat');

% Baseline-clean settings (same as Option-A script)
MIN_ISI_FROM_PREV_LICK_MS = 1000;

% Defaults if not stored
dt_ms = 10;
BASELINE_PRE_CUE_MS = 500;

% Learning transition session (for shading/line; set to [] to disable)
sStar = 8;

% Output files
outPng = fullfile(outDir, 'SESSION_EVENT_MODULATION_KERNEL_PLUS_HZ_STYLED.png');
outSvg = fullfile(outDir, 'SESSION_EVENT_MODULATION_KERNEL_PLUS_HZ_STYLED.svg');
outMat = fullfile(outDir, 'SESSION_EVENT_MODULATION_KERNEL_PLUS_HZ_STYLED.mat');

rng(0);

%% ---------------- STYLE (EXACTLY like FIG_SESSION_EPOCH_FIRINGRATE_HZ__GLOBALWINDOWS_STYLED__A_THESIS_FINAL.m) ----------------
figPos = [120, 120, 900, 520];

dataLW      = 5;    % line width
dataMS      = 10;   % marker size
markerLW    = 2.5;  % marker edge width
eventLineLW = 5;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

axesTickLW  = 4.0;

% ---- TICK / LABEL GEOMETRY ----
majorLenFrac      = 0.060;   % tick length as fraction of y-span
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;  % move tick numbers down (fraction of y-span)
xLabelDownFrac    = 0.08;    % move "Session" further down (fraction of y-span)

% ---- COLORS (use the ones from your Option-A script) ----
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

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

%% ---------------- KERNEL METRIC (PER SESSION) ----------------
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

dCue_k_sess(keepCue)     = mean(Kcue(keepCue, idxCueWin),   2, 'omitnan');
dPress_k_sess(keepPress) = mean(Kpress(keepPress, idxPressWin), 2, 'omitnan');
dLick_k_sess(keepLick)   = mean(Klick(keepLick, idxLickWin), 2, 'omitnan');

%% ---------------- HZ METRIC (PER SESSION; POPULATION POOL; BASELINE-CLEAN) ----------------
S = load(matFileDefault);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in: %s', matFileDefault);
nSessions = numel(beh);

% IMPORTANT: we assume GLM sessions match beh sessions.
% If they don't, you must align by date/session id and reorder one side.
if nSess_glm ~= nSessions
    warning('GLM has %d sessions but beh has %d sessions. Plotting will use nSessions=min(...).', nSess_glm, nSessions);
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

%% ========================= PLOT (STYLED) =========================
xAll = (1:nS)';

% ---- FIGURE 1: Kernel metric (per session) ----
plot_session_lines_styled_( ...
    xAll, dCue_k_sess(1:nS), dPress_k_sess(1:nS), dLick_k_sess(1:nS), ...
    colCue, colPress, colLick, figPos, dataLW, dataMS, markerLW, eventLineLW, ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, ...
    majorLenFrac, minorLenFrac, tickLabelDownFrac, xLabelDownFrac, ...
    sStar, ...
    'Session-level event modulation (Kernel metric; window mean; baseline=0)', ...
    'Response (kernel)', outPng, outSvg);

% ---- FIGURE 2: Hz metric (per session) ----
outPng2 = strrep(outPng, '.png', '_HZ.png');
outSvg2 = strrep(outSvg, '.svg', '_HZ.svg');

plot_session_lines_styled_( ...
    xAll, dCue_hz_sess, dPress_hz_sess, dLick_hz_sess, ...
    colCue, colPress, colLick, figPos, dataLW, dataMS, markerLW, eventLineLW, ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, ...
    majorLenFrac, minorLenFrac, tickLabelDownFrac, xLabelDownFrac, ...
    sStar, ...
    'Session-level event modulation (Hz metric; event Hz minus baseline Hz)', ...
    'Response (Hz \Delta)', outPng2, outSvg2);

%% ---- SAVE MAT ----
save(outMat, ...
    'dCue_k_sess','dPress_k_sess','dLick_k_sess', ...
    'dCue_hz_sess','dPress_hz_sess','dLick_hz_sess', ...
    'wCue_ms','wPress_ms','wLick_ms', ...
    'dt_ms','BASELINE_PRE_CUE_MS','MIN_ISI_FROM_PREV_LICK_MS', ...
    'glmSaveFile','matFileDefault','winFile');

fprintf('\nSaved MAT: %s\n', outMat);
fprintf('Saved KERNEL PNG/SVG: %s | %s\n', outPng, outSvg);
fprintf('Saved HZ     PNG/SVG: %s | %s\n', outPng2, outSvg2);
fprintf('\nDONE.\n');

%% =============================== HELPERS ===============================

function plot_session_lines_styled_( ...
    xAll, yCue, yPress, yLick, ...
    colCue, colPress, colLick, figPos, dataLW, dataMS, markerLW, eventLineLW, ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, ...
    majorLenFrac, minorLenFrac, tickLabelDownFrac, xLabelDownFrac, ...
    sStar, ttl, ylab, outPng, outSvg)

    nSessions = numel(xAll);

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    % y-limits based on data (clean, stable) (same logic)
    yAll = [yCue(:); yPress(:); yLick(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll)
        yl = [0 1];
    else
        yMin = min(yAll);
        yMax = max(yAll);
        pad  = 0.08 * max(eps, (yMax - yMin));
        yl   = [yMin - pad, yMax + pad];
        if yl(2) <= yl(1), yl = [yl(1)-1, yl(1)+1]; end
    end

    set(ax, 'XLim',[1 nSessions], 'YLim',yl);

    % -------------------- BACKGROUND SHADING --------------------
    if ~isempty(sStar) && isfinite(sStar)
        yl_patch = yl;
        patch(ax, [1 (sStar+0.5) (sStar+0.5) 1], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
            [1 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % yellowish
        patch(ax, [(sStar+0.5) nSessions nSessions (sStar+0.5)], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
            [0 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % greenish
    end
    % ------------------------------------------------------------

    % --- Lines ---
    plot(ax, xAll, yCue,   '-', 'LineWidth', dataLW, 'Color', colCue);
    plot(ax, xAll, yPress, '-', 'LineWidth', dataLW, 'Color', colPress);
    plot(ax, xAll, yLick,  '-', 'LineWidth', dataLW, 'Color', colLick);

    % --- Filled circles ---
    plot(ax, xAll, yCue,   'o', 'LineStyle','none', 'MarkerSize', dataMS, 'MarkerFaceColor', colCue,   'MarkerEdgeColor', colCue,   'LineWidth', markerLW);
    plot(ax, xAll, yPress, 'o', 'LineStyle','none', 'MarkerSize', dataMS, 'MarkerFaceColor', colPress, 'MarkerEdgeColor', colPress, 'LineWidth', markerLW);
    plot(ax, xAll, yLick,  'o', 'LineStyle','none', 'MarkerSize', dataMS, 'MarkerFaceColor', colLick,  'MarkerEdgeColor', colLick,  'LineWidth', markerLW);

    % --- Transition line ---
    if ~isempty(sStar) && isfinite(sStar)
        xline(ax, sStar+0.5, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);
    end

    hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
    ylabel(ax, ylab, 'FontSize', labelFontSize);
    title(ax, ttl, 'FontWeight','bold', 'FontSize', titleFontSize);

    set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSessions], 'YLim',yl);
    set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

    %% ---- TICKS: label every 5th; manual ticks + labels (same as your styled script) ----
    ax.XTick = 1:nSessions;
    labIdx = 5:5:nSessions;

    ax.XTickLabel = repmat({''}, 1, nSessions);
    xtickangle(ax, 0);

    set(ax, 'TickLength', [0 0]);

    y0 = ax.YLim(1);
    yr = diff(ax.YLim);

    majorLen = majorLenFrac * yr;
    minorLen = minorLenFrac * yr;

    tickCol = [0 0 0];

    for s = 1:nSessions
        if ismember(s, labIdx)
            L  = majorLen;
            lw = axesTickLW;
        else
            L  = minorLen;
            lw = max(1.5, axesTickLW * 0.55);
        end
        line(ax, [s s], [y0, y0 - L], 'Color', tickCol, 'LineWidth', lw, 'Clipping','off');
    end

    yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
    for s = labIdx
        text(ax, s, yTickText, num2str(s), ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontSize', tickFontSize, ...
            'Color', [0 0 0], ...
            'Clipping','off');
    end

    xCenter = mean(ax.XLim);
    yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
    set(hxlab, 'Units','data', 'Position',[xCenter, yXlab, 0]);

    ax.Position = [0.12 0.22 0.84 0.70];
    set(ax, 'XLim',[1 nSessions], 'YLim',yl);

    %% ---- SAVE ----
    saveas(fig, outPng);
    saveas(fig, outSvg);
    close(fig);

    fprintf('Saved: %s\n', outPng);
    fprintf('Saved: %s\n', outSvg);
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
        spikes = S.spikes;
    end
end