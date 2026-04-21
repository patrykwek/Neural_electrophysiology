%% ========= FIG_SESSION_WHOLETRIAL_FIRINGRATE_HZ__STYLED__A_THESIS_FINAL.m =========
% Same plot style as:
%   PLOT_PERCENT_CORRECT_WEIBULL_FIT_ASYMPTOTE_STYLED__A_THESIS_FINAL.m
% but computes ONE firing-rate value per session:
%   Whole-trial firing rate (Hz) over fixedWin = [-1000, 6000] ms relative to cue
% (cue-aligned SDF axis, no warping).
%
% Notes:
% - FR is computed per trial as mean SDF (Hz) across the entire fixedWin window
% - Per-unit FR summary = median(FR_trial_wholewin) across trials
% - Per-session FR summary = median(per-unit median(FR_trial_wholewin)) across units
%
% Saves:
%   PNG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.png
%   SVG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.svg
%   MAT: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.mat
%
% ADDED:
%   Early vs Late session-level bar graph (mean +/- SEM) with grey dot scatter:
%     early = sessions 1-8
%     late  = sessions 29-36
%   Y axis = FR (Hz)
%   X axis = {Early, Late}
%   Includes two-sample t-test; prints results to command window
%   If significant, annotate with stars
%   PNG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_WHOLETRIAL_FIRINGRATE_HZ_EARLY_LATE_BAR_STYLED.png
%   SVG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_WHOLETRIAL_FIRINGRATE_HZ_EARLY_LATE_BAR_STYLED.svg
%
% MAT contains:
%   medFR, medFR_s, nUnits, baseFR, perfPct, fixedWin, dt_ms, gaussSigmaMs, minTrialsPerUnit, sStar
%
% NEW (REQUESTED):
%   Scatter of performance vs neural metric (session medFR), each dot=session,
%   overlay STRAIGHT LINE (least-squares) and compute Spearman correlation.
%   Print Spearman results to command window. Save PNG + SVG.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE';

outPng  = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.png');
outSvg  = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.svg');
outMat  = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_STYLED.mat');

% --- ADDED outputs (early vs late bar) ---
outBarPng = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_EARLY_LATE_BAR_STYLED.png');
outBarSvg = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_EARLY_LATE_BAR_STYLED.svg');

% --- NEW outputs (performance vs FR scatter) ---
outScatPng = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_VS_PERFORMANCE_SCATTER.png');
outScatSvg = fullfile(baseOut, 'SESSION_WHOLETRIAL_FIRINGRATE_HZ_VS_PERFORMANCE_SCATTER.svg');

% >>> REQUIRED: load classification output <<<
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skipped)

% Analysis window used to BUILD SDFs (kept; we ensure coverage below)
preCueMs   = 500;
postLickMs = 3000;

% Whole-trial window (relative to cue)
fixedWin = [-1000, 6000];

% Trial filters
MIN_RT_MS    = 100;
requireValid = true;

% PER-UNIT inclusion
minTrialsPerUnit = 10;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Learning transition session (from Weibull)
sStar = 8;

% Baseline window (relative cue) -- computed only for sanity, not plotted
baselineWin_relCue = [-500, -100];  % ms

rng(0);

%% ---- STYLE (matched to Weibull-style script) ----
figPos = [120, 120, 900, 520];

dataLW      = 5;    % line width
dataMS      = 10;   % marker size
markerLW    = 2.5;  % marker edge width
eventLineLW = 5;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

axesTickLW  = 4.0;

%% ---- TICK / LABEL GEOMETRY ----
majorLenFrac      = 0.060;   % tick length as fraction of y-span
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;  % move tick numbers down (fraction of y-span)
xLabelDownFrac    = 0.08;    % move "Session" further down (fraction of y-span)

%% ---- COLOR (single series) ----
colFR = [0 0.45 0.75];  % match the "data" blue from your Weibull script

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

%% ---- PRECOMPUTE GAUSS KERNEL (unit area; Hz after dividing by dt) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- BUILD A SINGLE CUE-ALIGNED TIME AXIS THAT GUARANTEES WINDOW COVERAGE ----
% For whole-trial FR we only need fixedWin (plus baseline window for sanity).
minLag = min([baselineWin_relCue(1), fixedWin(1)]);
maxLag = max([baselineWin_relCue(2), fixedWin(2)]);

winLeft  = fixedWin(1) + minLag;   % conservative; ensures baseline coverage too
winRight = fixedWin(2) + maxLag;

% keep same safeguards as your other scripts
winLeft  = min(winLeft, -preCueMs);
winRight = max(winRight, fixedWin(2));

nT = numel(winLeft:dt_ms:winRight);
fprintf('\nSDF axis: [%d, %d] ms (dt=%d ms, nT=%d)\n', round(winLeft), round(winRight), dt_ms, nT);

%% ---- COMPUTE PER-SESSION METRICS (Hz; whole trial) ----
% Outputs:
% medFR(sess,1) : per-session median across units of per-unit median(FR_trial_whole) (Hz)
% nUnits(sess)  : number of included units
% baseFR(sess)  : mean across units of per-unit mean baseline FR (Hz) (sanity check)
% perfPct(sess) : optional, percent-correct estimate if detectable, else NaN
[medFR, nUnits, baseFR, perfPct] = ...
    computeSessionWholeTrialFR_Hz_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        baselineWin_relCue);

%% ---- SMOOTHING (ADDITION ONLY) ----
SMOOTH_WIN = 5;  % sessions (odd recommended)
medFR_s = medFR;
medFR_s(:) = smoothdata(medFR(:), 'movmean', SMOOTH_WIN);

%% ---- PLOT (STYLED) ----
xAll = (1:nSessions)';

fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

% y-limits based on data
yAll = medFR_s(:);
yAll = yAll(isfinite(yAll));
if isempty(yAll)
    yl = [0 1];
else
    yMin = min(yAll);
    yMax = max(yAll);
    pad  = 0.08 * max(eps, (yMax - yMin));
    yl   = [max(0, yMin - pad), yMax + pad];
    if yl(2) <= yl(1), yl = [max(0, yl(1)-1), yl(1)+1]; end
end

set(ax, 'XLim',[1 nSessions], 'YLim',yl);

% -------------------- BACKGROUND SHADING --------------------
% early (yellowish): [1, sStar+0.5]
% late  (greenish):  [sStar+0.5, nSessions]
yl_patch = yl;

patch(ax, [1 (sStar+0.5) (sStar+0.5) 1], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
    [1 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % yellowish

patch(ax, [(sStar+0.5) nSessions nSessions (sStar+0.5)], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
    [0 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % greenish
% ------------------------------------------------------------

% --- Line ---
plot(ax, xAll, medFR_s, '-', 'LineWidth', dataLW, 'Color', colFR);

% --- Markers ---
plot(ax, xAll, medFR_s, 'o', 'LineStyle','none', ...
    'MarkerSize', dataMS, 'MarkerFaceColor', colFR, 'MarkerEdgeColor', colFR, 'LineWidth', markerLW);

% --- Event line (learning transition) ---
xline(ax, sStar+0.5, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Whole-trial firing rate (Hz)', 'FontSize', labelFontSize);
title(ax, 'Per-session whole-trial firing rate', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSessions], 'YLim',yl);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

%% ---- TICKS: label every 5th; draw ticks + labels manually (Weibull style) ----
ax.XTick = 1:nSessions;
labIdx = 5:5:nSessions;

% blank built-in labels; we draw them manually
ax.XTickLabel = repmat({''}, 1, nSessions);
xtickangle(ax, 0);

% Hide axis' own tick marks
set(ax, 'TickLength', [0 0]);

% Manual ticks (DATA units)
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

% Manual tick labels (moved down)
yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
for s = labIdx
    text(ax, s, yTickText, num2str(s), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize', tickFontSize, ...
        'Color', [0 0 0], ...
        'Clipping','off');
end

% ---- MOVE X-LABEL DOWN BELOW TICK NUMBERS ----
xCenter = mean(ax.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlab, 'Units','data', 'Position',[xCenter, yXlab, 0]);

% Give more bottom room so nothing clips
ax.Position = [0.12 0.22 0.84 0.70];

% Ensure limits stay fixed
set(ax, 'XLim',[1 nSessions], 'YLim',yl);

%% ---- SAVE OUTPUTS ----
saveas(fig, outPng);
saveas(fig, outSvg);

%% ================= NEW: PERFORMANCE vs FR SCATTER + straight line overlay + Spearman =================
% IMPORTANT: Match Fano script performance computation exactly:
% performance = % correct from Hit or trials.Hit (NO perfPct fallback logic).

perfObs = nan(nSessions,1);
for sIdx = 1:nSessions
    sess = beh(sIdx);
    hit = [];
    if isfield(sess,'Hit') && ~isempty(sess.Hit)
        hit = sess.Hit;
    elseif isfield(sess,'trials') && isstruct(sess.trials) && isfield(sess.trials,'Hit')
        try
            hit = [sess.trials.Hit];
        catch
            hit = [];
        end
    end
    if isempty(hit), continue; end
    hit = double(hit(:));
    hit = hit(isfinite(hit));
    if isempty(hit), continue; end
    hit = hit(hit==0 | hit==1);
    if isempty(hit), continue; end
    perfObs(sIdx) = 100 * mean(hit==1);
end

% Use unsmoothed session metric for stats/scatter
metricY = medFR(:);

okSc = isfinite(perfObs) & isfinite(metricY);
xPerf = perfObs(okSc);
yFR   = metricY(okSc);

% Spearman correlation (print to command window)
if numel(xPerf) >= 3
    [rhoS, pS] = corr(xPerf, yFR, 'Type', 'Spearman', 'Rows', 'complete');
else
    rhoS = nan; pS = nan;
end

fprintf('\n===== Spearman: performance vs session medFR (unsmoothed) =====\n');
fprintf('n=%d sessions, rho=%.4f, p=%.6g\n', numel(xPerf), rhoS, pS);

% Straight line overlay (least-squares)
hasLine = numel(xPerf) >= 2;
if hasLine
    pLin  = polyfit(xPerf, yFR, 1);
    xLine = linspace(min(xPerf), max(xPerf), 200);
    yLine = polyval(pLin, xLine);
end

% Plot
figPosSc = figPos;
figPosSc(3) = round(figPos(3) * 0.55);  % modest width
figS = figure('Color','w','Position',figPosSc);
axS = axes(figS); hold(axS,'on');

% scatter
plot(axS, xPerf, yFR, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.6 0.6 0.6], ...
    'MarkerEdgeColor', [0.6 0.6 0.6], ...
    'LineWidth', 1.0);

% straight line overlay
if hasLine
    plot(axS, xLine, yLine, '-', 'LineWidth', axesTickLW, 'Color', colFR);
end

xlabel(axS, 'Performance (% correct)', 'FontSize', labelFontSize);
ylabel(axS, 'Whole-trial firing rate (Hz)', 'FontSize', labelFontSize);
title(axS, sprintf('Performance vs firing rate (Spearman \\rho=%.2f, p=%.3g)', rhoS, pS), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axS, 'FontSize', tickFontSize, 'Box','off');
set(axS, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Limits with padding
if ~isempty(xPerf)
    xMin = min(xPerf); xMax = max(xPerf);
    xPad = 0.06 * max(eps, (xMax-xMin));
    set(axS, 'XLim', [xMin-xPad, xMax+xPad]);
end
if ~isempty(yFR)
    yMin = min(yFR); yMax = max(yFR);
    yPad = 0.08 * max(eps, (yMax-yMin));
    set(axS, 'YLim', [max(0, yMin-yPad), yMax+yPad]);
end

axS.Position = [0.15 0.18 0.80 0.74];

saveas(figS, outScatPng);
saveas(figS, outScatSvg);

fprintf('Saved: %s\n', outScatPng);
fprintf('Saved: %s\n', outScatSvg);

%% ================= ADDED: EARLY vs LATE BAR (mean +/- SEM) + grey scatter + t-test =================
% Early: sessions 1-8
% Late : sessions 29-36
earlySess = 1:8;
lateSess  = 29:36;

% clip to available sessions
earlySess = earlySess(earlySess >= 1 & earlySess <= nSessions);
lateSess  = lateSess(lateSess  >= 1 & lateSess  <= nSessions);

% Use session-level values (unsmoothed medFR) for the comparison
earlyFR = medFR(earlySess);
lateFR  = medFR(lateSess);

% remove NaNs
earlyFR = earlyFR(isfinite(earlyFR));
lateFR  = lateFR(isfinite(lateFR));

% ------------------- CHANGED: permutation test instead of ttest2 -------------------
Nperm = 10000;
Tobs  = mean(lateFR, 'omitnan') - mean(earlyFR, 'omitnan');

pooled = [earlyFR(:); lateFR(:)];
nE = numel(earlyFR);
nL = numel(lateFR);

Tperm = nan(Nperm,1);
for ip = 1:Nperm
    permIdx = randperm(nE + nL);
    eIdx = permIdx(1:nE);
    lIdx = permIdx(nE+1:end);
    Tperm(ip) = mean(pooled(lIdx), 'omitnan') - mean(pooled(eIdx), 'omitnan');
end

p_t = (1 + sum(abs(Tperm) >= abs(Tobs))) / (1 + Nperm);

% keep original variable names for downstream save()
h_t = double(p_t < 0.05);
ci_t = [nan nan];
stats_t = struct('tstat', nan, 'df', nan, 'sd', nan);

% print results to command window
fprintf('\n===== Early vs Late permutation test (session-level medFR) =====\n');
fprintf('Early sessions: %s (n=%d)\n', mat2str(earlySess), numel(earlyFR));
fprintf('Late  sessions: %s (n=%d)\n', mat2str(lateSess),  numel(lateFR));
fprintf('Permutation test (two-sided, N=%d): p=%.6g\n', Nperm, p_t);
% -------------------------------------------------------------------------------

% significance stars
if ~isfinite(p_t)
    sigStars = '';
elseif p_t < 1e-3
    sigStars = '***';
elseif p_t < 1e-2
    sigStars = '**';
elseif p_t < 0.05
    sigStars = '*';
else
    sigStars = '';
end

% colors matched to your background shading
colEarly = [1 1 0];  % yellowish
colLate  = [0 1 0];  % greenish

% mean +/- SEM
mEarly = mean(earlyFR, 'omitnan');
mLate  = mean(lateFR,  'omitnan');

semEarly = std(earlyFR, 'omitnan') / sqrt(max(1, numel(earlyFR)));
semLate  = std(lateFR,  'omitnan') / sqrt(max(1, numel(lateFR)));

% ---- CHANGED ONLY: make BAR FIGURE 40% of previous width ----
figPosBar = figPos;
figPosBar(3) = round(figPos(3) * 0.40);  % 40% width, same height

figB = figure('Color','w','Position',figPosBar);
axB = axes(figB); hold(axB,'on');

% bar positions
xPos = [1 2];
barW = 0.60;

% bars
b = bar(axB, xPos, [mEarly mLate], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axesTickLW);
b.CData(1,:) = colEarly;
b.CData(2,:) = colLate;
b.FaceAlpha  = 0.25;

% errorbars (SEM)
errorbar(axB, xPos, [mEarly mLate], [semEarly semLate], ...
    'k', 'LineStyle','none', 'LineWidth', axesTickLW, 'CapSize', 18);

% scatter: grey dots (jittered)
jit = 0.10;
greyDot = [0.6 0.6 0.6];

if ~isempty(earlyFR)
    xj = 1 + (rand(size(earlyFR))-0.5)*2*jit;
    plot(axB, xj, earlyFR, 'o', 'MarkerSize', 8, 'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end
if ~isempty(lateFR)
    xj = 2 + (rand(size(lateFR))-0.5)*2*jit;
    plot(axB, xj, lateFR, 'o', 'MarkerSize', 8, 'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end

% Axes labels/ticks: X={Early,Late}, Y=FR
set(axB, 'XLim', [0.5 2.5]);
axB.XTick = [1 2];
axB.XTickLabel = {'Early','Late'};

ylabel(axB, 'Whole-trial firing rate (Hz)', 'FontSize', labelFontSize);
title(axB, 'Early vs late session firing rates', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axB, 'FontSize', tickFontSize, 'Box','off');
set(axB, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Y-limits with padding
yB = [earlyFR(:); lateFR(:); (mEarly+semEarly); (mLate+semLate)];
yB = yB(isfinite(yB));
if isempty(yB)
    ylB = [0 1];
else
    yMinB = min(yB);
    yMaxB = max(yB);
    padB  = 0.12 * max(eps, (yMaxB - yMinB));
    ylB   = [max(0, yMinB - padB), yMaxB + padB];
    if ylB(2) <= ylB(1), ylB = [max(0, ylB(1)-1), ylB(1)+1]; end
end
set(axB, 'YLim', ylB);

% significance bar + stars (only if significant)
if ~isempty(sigStars)
    yTop = ylB(2);
    yrB  = diff(ylB);
    yBar = yTop - 0.12*yrB;
    barH = 0.03*yrB;

    line(axB, [1 1 2 2], [yBar yBar+barH yBar+barH yBar], ...
        'Color', [0 0 0], 'LineWidth', axesTickLW);

    text(axB, 1.5, yBar+barH + 0.02*yrB, sigStars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', tickFontSize, 'FontWeight','bold', 'Color', [0 0 0]);
end

% similar margins
axB.Position = [0.12 0.18 0.84 0.74];

saveas(figB, outBarPng);
saveas(figB, outBarSvg);

fprintf('\nSaved: %s\n', outBarPng);
fprintf('Saved: %s\n', outBarSvg);

%% ---- SAVE MAT (include early/late vectors + test result) ----
save(outMat, ...
    'medFR','medFR_s','nUnits','baseFR','perfPct', ...
    'fixedWin','dt_ms','gaussSigmaMs','minTrialsPerUnit','sStar', ...
    'MIN_RT_MS','requireValid','baselineWin_relCue', ...
    'winLeft','winRight','nT','SMOOTH_WIN', ...
    'earlySess','lateSess','earlyFR','lateFR','p_t', ...
    'perfObs','rhoS','pS');

fprintf('\nSaved: %s\n', outPng);
fprintf('Saved: %s\n', outSvg);
fprintf('Saved: %s\n', outMat);

%% ================= LOCAL HELPERS =================

function [medFR, nUnits, baseFR, perfPct] = ...
    computeSessionWholeTrialFR_Hz_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        baselineWin_relCue)

    nSessions = numel(beh);

    medFR   = nan(nSessions,1);
    nUnits  = nan(nSessions,1);
    baseFR  = nan(nSessions,1);
    perfPct = nan(nSessions,1);

    for sIdx = 1:nSessions
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session)
            continue;
        end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
            continue;
        end

        % ---- Behavior-only keepTrial (same logic as your MI code) ----
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
        if isempty(idxKeep)
            continue;
        end

        cueAbsK   = cueAbs(idxKeep);
        pressLatK = pressLat(idxKeep); %#ok<NASGU>
        lickLatK  = lickLat(idxKeep);  %#ok<NASGU>

        % Optional: attempt to compute %correct if there is an obvious field
        perfPct(sIdx) = tryComputePerfPct_(trials(idxKeep));

        % Per-unit summaries
        unitMedFR  = nan(0,1);
        unitBaseFR = nan(0,1);

        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end

            for u = 1:numel(uc)
                % only classified units
                ct = safe_get_celltype_(cell_type, sIdx, ch, u);
                if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                    continue;
                end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                % Active-trial filter: >=1 spike in [cue+fixedWin]
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
                nTrU = numel(cueU);

                FR_whole_tr = nan(nTrU,1);
                baseFR_tr   = nan(nTrU,1);

                for iTr = 1:nTrU
                    cue0 = cueU(iTr);

                    % SDF (Hz)
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

                    % Indices for baseline (relative cue) and whole-trial fixedWin
                    idxBase  = msWindowToIdx_(baselineWin_relCue, winLeft, dt_ms, nT);
                    idxWhole = msWindowToIdx_(fixedWin,            winLeft, dt_ms, nT);

                    frBase  = mean(y(idxBase(1):idxBase(2)),   'omitnan');
                    frWhole = mean(y(idxWhole(1):idxWhole(2)), 'omitnan');

                    baseFR_tr(iTr)   = frBase;
                    FR_whole_tr(iTr) = frWhole;
                end

                % Per-unit summary (Hz): median across trials
                mWhole = median(FR_whole_tr, 'omitnan');

                unitMedFR(end+1,1)  = mWhole; %#ok<AGROW>
                unitBaseFR(end+1,1) = mean(baseFR_tr, 'omitnan'); %#ok<AGROW>
            end
        end

        if isempty(unitMedFR)
            continue;
        end

        nU = numel(unitMedFR);
        nUnits(sIdx) = nU;

        % Per-session summary: median across units
        medFR(sIdx,1) = median(unitMedFR, 'omitnan');
        baseFR(sIdx)  = mean(unitBaseFR, 'omitnan');
    end
end

function pct = tryComputePerfPct_(trials)
    % Optional best-effort: returns %correct if there is a clearly named boolean field.
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