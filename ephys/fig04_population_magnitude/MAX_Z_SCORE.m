%% ========= FIG_MAXTRIALZ_BY_SESSION__GLOBALWARP__WHOLETRIALZ__CLASSIFIEDONLY__WEIBULLAXISSTYLE__PLUS_BAR.m =========
% PURPOSE:
%   (1) Session-by-session plot of "Max trial z modulation" using SAME method as your
%       z-score modulation script:
%         - SDF per trial
%         - GLOBAL warp (median RT and lick latency across ALL sessions)
%         - Whole-trial z-score per trial
%         - Trial metric = max(z(t)) across time
%         - Unit metric = mean(trial metric) across trials
%         - Session metric = mean(unit metric) across units
%       Plot style: Weibull ticks (label every 5th starting at 5), early/late shading swapped
%       BUT line + dots are BLUE like your Weibull script.
%
%   (2) Early vs Late bar plot (mean +/- SEM) + grey scatter (session-level values),
%       with two-sample t-test like your FR script.
%
%   (3) NEW (REQUESTED): Performance vs max-z scatter:
%       - x-axis: performance (% correct; observed)
%       - y-axis: sessMaxZ (session metric)
%       - dot = session
%       - overlay straight line (least-squares)
%       - Spearman rho and p printed to command window
%       - saves PNG + SVG
%
% SAVES:
%   - Session curve: PNG + SVG + FIG + MAT
%   - Bar plot:      PNG + SVG + FIG
%   - Perf scatter:  PNG + SVG
%
% NOTES:
%   - Uses ONLY classified units (cell type 1/2/3), and requires >=minTrialsPerUnit
%     active trials per unit (active = spike in [cue+fixedWin(1), cue+lickLat]).
%   - Session-level values used in t-test/bar plot are sessMaxZ (unsmoothed).

clear; clc;

%% ---- USER SETTINGS ----
matFile      = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE';
if ~exist(baseOut,'dir'), mkdir(baseOut); end

% outputs (session curve)
outPng = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_STYLED.png');
outSvg = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_STYLED.svg');
outFig = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_STYLED.fig');
outMat = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_STYLED.mat');

% outputs (early vs late bar)
outBarPng = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_EARLY_LATE_BAR_STYLED.png');
outBarSvg = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_EARLY_LATE_BAR_STYLED.svg');
outBarFig = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_EARLY_LATE_BAR_STYLED.fig');

% --- NEW outputs (performance vs maxZ scatter) ---
outScatPng = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_VS_PERFORMANCE_SCATTER.png');
outScatSvg = fullfile(baseOut, 'SESSION_MAXTRIALZ_MODULATION_VS_PERFORMANCE_SCATTER.svg');

% Windowing (match your z-trace script)
preCueMs   = 500;
postLickMs = 3000;

% Used ONLY for defining per-unit spike window (activeTrials check)
fixedWin = [-1000, 6000];   % ms relative to cue (fixedWin(2) not used for inclusion)

% Trial filters (match your z-trace script)
MIN_RT_MS         = 100;
minTrialsPerUnit  = 10;
requireValid      = true;

% SDF settings
dt_ms        = 10;

% Smoothing kernel (same function as your script)
gaussSigmaMs = 50;          % keep consistent with your current script

% Post-warp smoothing to reduce interpolation edge sharpening
doPostWarpSmooth = true;

% Session groups
earlyGroup = 1:8;
lateGroup  = 29:36;
boundaryAt = 8.5;

rng(0);

%% ---- STYLE (match Weibull axis style) ----
figPos   = [120, 120, 1200, 650];

dataLW      = 5;
dataMS      = 10;
markerLW    = 2.5;
eventLineLW = 5;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;
axesTickLW    = 4.0;

majorLenFrac      = 0.060;
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;
xLabelDownFrac    = 0.08;

% Shading colors (swapped)
colEarlyShade = [1 1 0];     % Early = YELLOW
colLateShade  = [0 1 0];     % Late  = GREEN
bgAlpha   = 0.10;
dashColor = [0.2 0.2 0.2];

% Blue data color (like Weibull script)
colBlue = [0 0.45 0.75];

%% ---- LOAD BEHSTRUCT ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = numel(beh);

%% ---- LOAD CELL TYPES ----
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

%% ---- PRECOMPUTE SDF KERNEL (unit area) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (shared across all sessions) ----
all_pressLat_global = [];
all_lickLat_global  = [];

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

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

        lr = lick  - cue;

        all_pressLat_global(end+1,1) = rt; %#ok<AGROW>
        all_lickLat_global(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat_global) && ~isempty(all_lickLat_global), ...
    'No trials passed global filters (cannot set shared warp targets).');

Tpress_global = median(all_pressLat_global);
Tlick_global  = median(all_lickLat_global);

fprintf('\n=== GLOBAL targets: Tpress=%.1f ms, Tlick=%.1f ms ===\n', ...
    Tpress_global, Tlick_global);

%% ---- TEMPLATE AXIS (same as your z-trace script) ----
winLeft  = -preCueMs;
winRight_template = Tlick_global + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,            winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress_global, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick_global,  winLeft, dt_ms, nT);

%% ============================================================
%  COMPUTE SESSION METRIC: avg over units of (avg over trials of max(z(t)))
%% ============================================================
sessMaxZ   = nan(nSessions,1);
sessNUnits = zeros(nSessions,1);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Behavior-only keepTrial: valid, RT>=MIN_RT ----
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

        lr = lick - cue;

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

    unitVals = []; % per-unit average max-trial-z

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)

            % ONLY KEEP CLASSIFIED UNITS (1/2/3)
            ct = safeCellType_(cell_type, sIdx, ch, u);
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                continue;
            end

            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
            spk_abs = double(spk_abs(:));

            % PER-UNIT activeTrials: spike in [cue + fixedWin(1), cue + lickLat]
            activeTrials = false(numel(idxKeep),1);
            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                lr0  = lickLatK(iTr);

                w0 = cue0 + fixedWin(1);
                w1 = cue0 + lr0;

                if any(spk_abs >= w0 & spk_abs <= w1)
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

            trialMax = nan(nTrU,1);

            for iTr = 1:nTrU
                cue0   = cueU(iTr);
                tPress = rtU(iTr);
                tLick  = lkU(iTr);

                winRight_trial = tLick + postLickMs;
                nT_trial = numel(winLeft:dt_ms:winRight_trial);

                idx0 = timeToIdx_(0,      winLeft, dt_ms, nT_trial);
                idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

                % SDF (Hz)
                spk_rel = spk_abs - cue0;
                spk_rel = spk_rel(spk_rel >= winLeft & spk_rel <= winRight_trial);

                counts = zeros(1, nT_trial);
                if ~isempty(spk_rel)
                    jj = round((spk_rel - winLeft)/dt_ms) + 1;
                    jj = jj(jj >= 1 & jj <= nT_trial);
                    for q = 1:numel(jj)
                        counts(jj(q)) = counts(jj(q)) + 1;
                    end
                end

                sm = conv(counts, g, 'same');
                y  = sm / (dt_ms/1000);

                % Warp cue->press, press->lick; keep pre-cue and post-lick unwarped
                pre = y(1:idx0);

                segA = y(idx0:idxP);
                kA   = max(2, (idxPT - idx0T + 1));
                wA   = warpSegment_floorceilEqn_(segA, kA);

                segB = y(idxP:idxL);
                kB   = max(2, (idxLT - idxPT + 1));
                wB   = warpSegment_floorceilEqn_(segB, kB);

                segC = y(idxL:end);
                kC_target = max(1, (nT - idxLT + 1));
                if numel(segC) < kC_target
                    segC_fit = [segC, zeros(1, kC_target - numel(segC))];
                else
                    segC_fit = segC(1:kC_target);
                end

                yy = [pre, wA(2:end), wB(2:end), segC_fit(2:end)];

                if numel(yy) < nT
                    yy = [yy, zeros(1, nT-numel(yy))];
                elseif numel(yy) > nT
                    yy = yy(1:nT);
                end

                if doPostWarpSmooth
                    yy = conv(yy, g, 'same');
                end

                z = zscoreTrial_(yy);
                trialMax(iTr) = max(z);
            end

            unitVals(end+1,1) = mean(trialMax, 'omitnan'); %#ok<AGROW>
        end
    end

    if ~isempty(unitVals)
        sessMaxZ(sIdx)   = mean(unitVals, 'omitnan');
        sessNUnits(sIdx) = numel(unitVals);
    end

    fprintf('Session %02d: nUnits=%d, avgMaxTrialZ=%.4f\n', ...
        sIdx, sessNUnits(sIdx), sessMaxZ(sIdx));
end

%% ================= NEW: PERFORMANCE vs maxZ SCATTER + straight line overlay + Spearman =================
% Performance = observed % correct from Hit/session.trials.Hit (same logic as your Weibull script).
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

okSc = isfinite(perfObs) & isfinite(sessMaxZ);
xPerf = perfObs(okSc);
yZ    = sessMaxZ(okSc);

if numel(xPerf) >= 3
    [rhoS, pS] = corr(xPerf, yZ, 'Type', 'Spearman', 'Rows', 'complete');
else
    rhoS = nan; pS = nan;
end

fprintf('\n===== Spearman: performance vs session avgMaxTrialZ =====\n');
fprintf('n=%d sessions, rho=%.4f, p=%.6g\n', numel(xPerf), rhoS, pS);

hasLine = numel(xPerf) >= 2;
if hasLine
    pLin  = polyfit(xPerf, yZ, 1);
    xLine = linspace(min(xPerf), max(xPerf), 200);
    yLine = polyval(pLin, xLine);
end

figPosSc = figPos;
figPosSc(3) = round(figPos(3) * 0.55); % modest width
figS = figure('Color','w','Position',figPosSc);
axS = axes(figS); hold(axS,'on');

plot(axS, xPerf, yZ, 'o', ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.6 0.6 0.6], ...
    'MarkerEdgeColor', [0.6 0.6 0.6], ...
    'LineWidth', 1.0);

if hasLine
    plot(axS, xLine, yLine, '-', 'LineWidth', axesTickLW, 'Color', colBlue);
end

xlabel(axS, 'Performance (% correct)', 'FontSize', labelFontSize);
ylabel(axS, 'Max trial z modulation', 'FontSize', labelFontSize);
title(axS, sprintf('Performance vs max trial z (Spearman \\rho=%.2f, p=%.3g)', rhoS, pS), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axS, 'FontSize', tickFontSize, 'Box','off');
set(axS, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Limits with padding
if ~isempty(xPerf)
    xMin = min(xPerf); xMax = max(xPerf);
    xPad = 0.06 * max(eps, (xMax-xMin));
    set(axS, 'XLim', [xMin-xPad, xMax+xPad]);
end
if ~isempty(yZ)
    yMin = min(yZ); yMax = max(yZ);
    yPad = 0.08 * max(eps, (yMax-yMin));
    set(axS, 'YLim', [yMin-yPad, yMax+yPad]);
end

axS.Position = [0.15 0.18 0.80 0.74];

saveas(figS, outScatPng);
saveas(figS, outScatSvg);

fprintf('Saved: %s\n', outScatPng);
fprintf('Saved: %s\n', outScatSvg);

%% ============================================================
%  (1) SESSION CURVE PLOT (Weibull style ticks; BLUE line+dots)
%% ============================================================
nSess = nSessions;
x = 1:nSess;

fig1 = figure('Color','w','Position',figPos);
ax = axes(fig1); hold(ax,'on');

finiteY = sessMaxZ(isfinite(sessMaxZ));
if isempty(finiteY)
    yl = [-1 1];
else
    pad = 0.05 * max(1, range(finiteY));
    yl = [min(finiteY)-pad, max(finiteY)+pad];
end
set(ax, 'XLim',[1 nSess], 'YLim', yl);

% background shading (swapped)
patch(ax, [1 boundaryAt boundaryAt 1], [yl(1) yl(1) yl(2) yl(2)], ...
    colEarlyShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(ax, [boundaryAt nSess nSess boundaryAt], [yl(1) yl(1) yl(2) yl(2)], ...
    colLateShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');

% data (BLUE)
plot(ax, x, sessMaxZ, '-', 'LineWidth', dataLW, 'Color', colBlue);
plot(ax, x, sessMaxZ, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', dataMS, ...
    'MarkerFaceColor', colBlue, ...
    'MarkerEdgeColor', colBlue, ...
    'LineWidth', markerLW);

% boundary
xline(ax, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Max trial z modulation', 'FontSize', labelFontSize);
title(ax, sprintf('Avg max trial z modulation by session (GLOBAL warp); minTrials/unit=%d', minTrialsPerUnit), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSess], 'YLim', yl);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% manual ticks: label every 5th starting at 5
ax.XTick = 1:nSess;
labIdx = 5:5:nSess;
ax.XTickLabel = repmat({''}, 1, nSess);
xtickangle(ax, 0);
set(ax, 'TickLength', [0 0]);

y0 = ax.YLim(1);
yr = diff(ax.YLim);
majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

for s = 1:nSess
    if ismember(s, labIdx)
        L  = majorLen;
        lw = axesTickLW;
    else
        L  = minorLen;
        lw = max(1.5, axesTickLW * 0.55);
    end
    line(ax, [s s], [y0, y0 - L], 'Color', [0 0 0], 'LineWidth', lw, 'Clipping','off');
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
set(ax, 'XLim',[1 nSess], 'YLim', yl);

% save curve fig
saveas(fig1, outPng);
saveas(fig1, outSvg);
savefig(fig1, outFig);

%% ============================================================
%  (2) EARLY vs LATE BAR PLOT + t-test (like your FR script)
%% ============================================================
earlySess = earlyGroup;
lateSess  = lateGroup;

earlySess = earlySess(earlySess >= 1 & earlySess <= nSessions);
lateSess  = lateSess(lateSess  >= 1 & lateSess  <= nSessions);

earlyVals = sessMaxZ(earlySess);
lateVals  = sessMaxZ(lateSess);

earlyVals = earlyVals(isfinite(earlyVals));
lateVals  = lateVals(isfinite(lateVals));

% ------------------- CHANGED: permutation test instead of ttest2 -------------------
Nperm = 10000;
Tobs  = mean(lateVals, 'omitnan') - mean(earlyVals, 'omitnan');

pooled = [earlyVals(:); lateVals(:)];
nE = numel(earlyVals);
nL = numel(lateVals);

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

fprintf('\n===== Early vs Late permutation test (session-level avgMaxTrialZ) =====\n');
fprintf('Early sessions: %s (n=%d)\n', mat2str(earlySess), numel(earlyVals));
fprintf('Late  sessions: %s (n=%d)\n', mat2str(lateSess),  numel(lateVals));
fprintf('Permutation test (two-sided, N=%d): p=%.6g\n', Nperm, p_t);
% -------------------------------------------------------------------------------

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

mEarly = mean(earlyVals, 'omitnan');
mLate  = mean(lateVals,  'omitnan');
semEarly = std(earlyVals, 'omitnan') / sqrt(max(1, numel(earlyVals)));
semLate  = std(lateVals,  'omitnan') / sqrt(max(1, numel(lateVals)));

% Bar figure position (use same convention as your FR script)
figPosBar = figPos;
figPosBar(3) = round(figPos(3) * 0.40);

fig2 = figure('Color','w','Position',figPosBar);
axB = axes(fig2); hold(axB,'on');

xPos = [1 2];
barW = 0.60;

b = bar(axB, xPos, [mEarly mLate], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axesTickLW);
b.CData(1,:) = colEarlyShade;
b.CData(2,:) = colLateShade;
b.FaceAlpha  = 0.25;

errorbar(axB, xPos, [mEarly mLate], [semEarly semLate], ...
    'k', 'LineStyle','none', 'LineWidth', axesTickLW, 'CapSize', 18);

% grey scatter (session points)
jit = 0.10;
greyDot = [0.6 0.6 0.6];

if ~isempty(earlyVals)
    xj = 1 + (rand(size(earlyVals))-0.5)*2*jit;
    plot(axB, xj, earlyVals, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end
if ~isempty(lateVals)
    xj = 2 + (rand(size(lateVals))-0.5)*2*jit;
    plot(axB, xj, lateVals, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end

set(axB, 'XLim', [0.5 2.5]);
axB.XTick = [1 2];
axB.XTickLabel = {'Early','Late'};

ylabel(axB, 'Max trial z modulation', 'FontSize', labelFontSize);
title(axB, 'Early vs late max trial z modulation', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axB, 'FontSize', tickFontSize, 'Box','off');
set(axB, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% y-lims
yB = [earlyVals(:); lateVals(:); (mEarly+semEarly); (mLate+semLate)];
yB = yB(isfinite(yB));
if isempty(yB)
    ylB = [-1 1];
else
    yMinB = min(yB);
    yMaxB = max(yB);
    padB  = 0.12 * max(eps, (yMaxB - yMinB));
    ylB   = [yMinB - padB, yMaxB + padB];
    if ylB(2) <= ylB(1), ylB = [ylB(1)-1, ylB(1)+1]; end
end
set(axB, 'YLim', ylB);

% significance annotation
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

axB.Position = [0.12 0.18 0.84 0.74];

saveas(fig2, outBarPng);
saveas(fig2, outBarSvg);
savefig(fig2, outBarFig);

%% ---- SAVE MAT ----
save(outMat, ...
    'sessMaxZ','sessNUnits', ...
    'earlySess','lateSess','earlyVals','lateVals', ...
    'h_t','p_t','ci_t','stats_t', ...
    'matFile','cellTypeFile','preCueMs','postLickMs','fixedWin','MIN_RT_MS','minTrialsPerUnit', ...
    'dt_ms','gaussSigmaMs','doPostWarpSmooth','Tpress_global','Tlick_global');

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n  %s\n  %s\n', outPng, outSvg, outFig, outBarPng, outBarSvg, outMat);
fprintf('\nDONE.\n');

%% ================= HELPERS =================
function z = zscoreTrial_(x)
    x = x(:)';
    mu = mean(x);
    sd = std(x);
    if ~isfinite(sd) || sd <= 0, sd = 1; end
    z = (x - mu) / sd;
end

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
    halfWidth = ceil(5*sigmaMs/dtMs);
    x = (-halfWidth:halfWidth) * dtMs;
    g = exp(-0.5*(x./sigmaMs).^2);
    g = g / sum(g);
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

function ct = safeCellType_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
            ch <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
            u <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end