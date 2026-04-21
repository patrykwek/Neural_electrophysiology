%% ========= FIG3X_COMPARE_CUE_vs_PRESS_DECODING__LATESESSIONS__A_THESIS_FINAL.m =========
% Compare incorrect-trial decoding accuracy:
%   cue decoder vs press decoder
% using saved outputs from the matched TRAINALL_TESTINCORRECT scripts.
%
% PLOTS
%   - Late sessions only
%   - Paired barplots:
%       * overall mean decoding accuracy
%       * cue-segment mean decoding accuracy
%       * press-segment mean decoding accuracy
%       * lick-segment mean decoding accuracy
%       * individual sessions shown as dots
%       * the same session connected across cue vs press by a line
%   - One session-level paired comparison scatter:
%       * one dot per late session
%       * x-axis = cue decoder accuracy
%       * y-axis = press decoder accuracy
%       * unity line
%   - One warped-trial trace:
%       * cue decoder mean +/- SEM across late sessions
%       * press decoder mean +/- SEM across late sessions
%       * cue / press / lick event lines
%       * paired cluster-based permutation test across time
%
% INPUTS REQUIRED
%   - Cue-decoder results MAT from:
%       FIG3X_SESSION_CUETYPE_DECODING_HEATMAP__TRAINALL_TESTINCORRECT__...
%   - Press-decoder results MAT from:
%       FIG3X_SESSION_PRESSTYPE_DECODING_HEATMAP__TRAINALL_TESTINCORRECT__...
%
% NOTES
%   - Uses session-level summary fields and AccMat from each MAT file
%   - Uses the intersection of sessions with finite values in both decoders
%   - Defaults to late sessions only (sessions 29:36, clipped to available session count)
%   - Styled to match your existing figures
%   - PCA of cue-space across sessions cannot be built from these saved outputs alone,
%     because trial-level / class-level neural feature data are not present in the MATs

clear; clc;

%% ---- USER SETTINGS ----
cueMat = '/Volumes/WD_BLACK/A_THESIS_FINAL/DECODING_CUETYPE_HEATMAP_TRAINALL_TESTINCORRECT/Sessions01_36_CueTypeDecodingResults_TRAINALL_TESTINCORRECT_GLOBALWARP_BIN200ms_train80_repeats050.mat';
pressMat = '/Volumes/WD_BLACK/A_THESIS_FINAL/DECODING_PRESSTYPE_HEATMAP_TRAINALL_TESTINCORRECT/Sessions01_36_PressTypeDecodingResults_TRAINALL_TESTINCORRECT_GLOBALWARP_BIN200ms_train80_repeats050.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL';
outDir  = fullfile(baseOut, 'DECODING_CUE_vs_PRESS_COMPARISON');
if ~exist(outDir,'dir'), mkdir(outDir); end

useLateOnly = true;   % requested focus
manualSess  = [];     % leave empty to use lateSess from files

% Figure style
figPos = [120, 120, 520, 700];
figPosTrace = [120, 120, 650, 300];

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;
cbFontSize    = 22;

axesTickLW  = 4.0;
axesTickLen = [0.02 0.02];
traceLW     = 4;
eventLineLW = 5.0;

scatterMS     = 8;
scatterEdgeLW = 1.0;
jit           = 0.06;
greyDot       = [0.6 0.6 0.6];

% EXACT SAME COLORS AS PROVIDED STYLE SCRIPT
colCueBar   = [1 0 0];
colPressBar = [0.10 0.55 0.95];
colLine     = [0.45 0.45 0.45];

colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

barFaceAlpha = 0.25;

Nperm = 10000;  % paired sign-flip permutation test for bars
rngSeed = 0;
rng(rngSeed);

% trace cluster permutation settings
Nperm_cluster = 10000;
clusterAlpha  = 0.05;

xUnit = "s";
forceIntegerSecondsTicks = true;

%% ---- LOAD RESULTS ----
assert(exist(cueMat,'file')==2,   'Cue MAT not found: %s', cueMat);
assert(exist(pressMat,'file')==2, 'Press MAT not found: %s', pressMat);

C = load(cueMat);
P = load(pressMat);

reqFields = {'sessMeanAcc','sessCueAcc','sessPressAcc','sessLickAcc','AccMat','nSessions','xPlot'};
for i = 1:numel(reqFields)
    assert(isfield(C, reqFields{i}), 'Cue MAT missing field: %s', reqFields{i});
    assert(isfield(P, reqFields{i}), 'Press MAT missing field: %s', reqFields{i});
end

assert(numel(C.sessMeanAcc)   == C.nSessions, 'Cue sessMeanAcc length mismatch.');
assert(numel(C.sessCueAcc)    == C.nSessions, 'Cue sessCueAcc length mismatch.');
assert(numel(C.sessPressAcc)  == C.nSessions, 'Cue sessPressAcc length mismatch.');
assert(numel(C.sessLickAcc)   == C.nSessions, 'Cue sessLickAcc length mismatch.');

assert(numel(P.sessMeanAcc)   == P.nSessions, 'Press sessMeanAcc length mismatch.');
assert(numel(P.sessCueAcc)    == P.nSessions, 'Press sessCueAcc length mismatch.');
assert(numel(P.sessPressAcc)  == P.nSessions, 'Press sessPressAcc length mismatch.');
assert(numel(P.sessLickAcc)   == P.nSessions, 'Press sessLickAcc length mismatch.');

nSessions = min(C.nSessions, P.nSessions);

cueOverallAll   = C.sessMeanAcc(1:nSessions);
pressOverallAll = P.sessMeanAcc(1:nSessions);

cueCueSegAll    = C.sessCueAcc(1:nSessions);
pressCueSegAll  = P.sessCueAcc(1:nSessions);

cuePressSegAll   = C.sessPressAcc(1:nSessions);
pressPressSegAll = P.sessPressAcc(1:nSessions);

cueLickSegAll   = C.sessLickAcc(1:nSessions);
pressLickSegAll = P.sessLickAcc(1:nSessions);

%% ---- CHOOSE SESSIONS ----
if ~isempty(manualSess)
    sessUse = manualSess(:)';
elseif useLateOnly
    if isfield(C,'lateSess') && ~isempty(C.lateSess)
        sessUse = C.lateSess(:)';
    elseif isfield(P,'lateSess') && ~isempty(P.lateSess)
        sessUse = P.lateSess(:)';
    else
        sessUse = 29:36;
    end
else
    sessUse = 1:nSessions;
end

sessUse = sessUse(sessUse >= 1 & sessUse <= nSessions);
assert(~isempty(sessUse), 'No sessions selected.');

fprintf('Requested sessions: %s\n', mat2str(sessUse));

%% ---- OVERALL BARPLOT ----
plotPairedBar_( ...
    cueOverallAll, pressOverallAll, sessUse, ...
    'Late sessions: overall incorrect-trial decoding', ...
    'Mean incorrect-trial decoding accuracy', ...
    fullfile(outDir, 'CUE_vs_PRESS_DECODING__OVERALL_BAR__LateSessions'), ...
    figPos, titleFontSize, labelFontSize, tickFontSize, ...
    axesTickLW, axesTickLen, scatterMS, scatterEdgeLW, jit, greyDot, ...
    colCueBar, colPressBar, colLine, barFaceAlpha, Nperm, rngSeed);

%% ---- CUE-SEGMENT BARPLOT ----
plotPairedBar_( ...
    cueCueSegAll, pressCueSegAll, sessUse, ...
    'Late sessions: cue-segment incorrect-trial decoding', ...
    'Cue-segment accuracy', ...
    fullfile(outDir, 'CUE_vs_PRESS_DECODING__CUESEGMENT_BAR__LateSessions'), ...
    figPos, titleFontSize, labelFontSize, tickFontSize, ...
    axesTickLW, axesTickLen, scatterMS, scatterEdgeLW, jit, greyDot, ...
    colCueBar, colPressBar, colLine, barFaceAlpha, Nperm, rngSeed);

%% ---- PRESS-SEGMENT BARPLOT ----
plotPairedBar_( ...
    cuePressSegAll, pressPressSegAll, sessUse, ...
    'Late sessions: press-segment incorrect-trial decoding', ...
    'Press-segment accuracy', ...
    fullfile(outDir, 'CUE_vs_PRESS_DECODING__PRESSSEGMENT_BAR__LateSessions'), ...
    figPos, titleFontSize, labelFontSize, tickFontSize, ...
    axesTickLW, axesTickLen, scatterMS, scatterEdgeLW, jit, greyDot, ...
    colCueBar, colPressBar, colLine, barFaceAlpha, Nperm, rngSeed);

%% ---- LICK-SEGMENT BARPLOT ----
plotPairedBar_( ...
    cueLickSegAll, pressLickSegAll, sessUse, ...
    'Late sessions: lick-segment incorrect-trial decoding', ...
    'Lick-segment accuracy', ...
    fullfile(outDir, 'CUE_vs_PRESS_DECODING__LICKSEGMENT_BAR__LateSessions'), ...
    figPos, titleFontSize, labelFontSize, tickFontSize, ...
    axesTickLW, axesTickLen, scatterMS, scatterEdgeLW, jit, greyDot, ...
    colCueBar, colPressBar, colLine, barFaceAlpha, Nperm, rngSeed);

%% ---- SESSION-LEVEL PAIRED SCATTER (OVERALL) ----
plotPairedScatter_( ...
    cueOverallAll, pressOverallAll, sessUse, ...
    'Late sessions: session-level cue vs press decoding', ...
    'Cue decoder accuracy', ...
    'Press decoder accuracy', ...
    fullfile(outDir, 'CUE_vs_PRESS_DECODING__OVERALL_SCATTER__LateSessions'), ...
    figPos, titleFontSize, labelFontSize, tickFontSize, ...
    axesTickLW, axesTickLen, scatterMS, scatterEdgeLW, ...
    colCueBar, colPressBar, greyDot);

%% ---- WARPED-TRIAL TRACE (LATE SESSIONS ONLY) ----
assert(isfield(C,'AccMat') && isfield(P,'AccMat'), 'AccMat missing in one of the MAT files.');
assert(isfield(C,'xPlot')  && isfield(P,'xPlot'),  'xPlot missing in one of the MAT files.');

nBins = min(size(C.AccMat,1), size(P.AccMat,1));
xPlot = C.xPlot(:)';
xPlot = xPlot(1:nBins);

sessTrace = sessUse;
sessTrace = sessTrace(sessTrace >= 1 & sessTrace <= size(C.AccMat,2) & sessTrace <= size(P.AccMat,2));

cueTraceMat   = C.AccMat(1:nBins, sessTrace)';   % nSess x nBins
pressTraceMat = P.AccMat(1:nBins, sessTrace)';   % nSess x nBins

validSessTrace = any(isfinite(cueTraceMat),2) & any(isfinite(pressTraceMat),2);
cueTraceMat   = cueTraceMat(validSessTrace,:);
pressTraceMat = pressTraceMat(validSessTrace,:);
sessTrace     = sessTrace(validSessTrace);

assert(~isempty(sessTrace), 'No overlapping finite late sessions for trace plot.');

mCueTrace   = mean(cueTraceMat,   1, 'omitnan');
mPressTrace = mean(pressTraceMat, 1, 'omitnan');

seCueTrace   = std(cueTraceMat,   0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(cueTraceMat),1)));
sePressTrace = std(pressTraceMat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(pressTraceMat),1)));

if isfield(C,'cueLine')
    cueLine = C.cueLine;
elseif isfield(P,'cueLine')
    cueLine = P.cueLine;
else
    cueLine = 0;
end

if isfield(C,'pressLine')
    pressLine = C.pressLine;
elseif isfield(P,'pressLine')
    pressLine = P.pressLine;
else
    pressLine = NaN;
end

if isfield(C,'lickLine')
    lickLine = C.lickLine;
elseif isfield(P,'lickLine')
    lickLine = P.lickLine;
else
    lickLine = NaN;
end

if isfield(C,'xPlot')
    xlab = 'Warped time from Cue (s)';
else
    xlab = 'Warped time';
end

%% ---- PAIRED CLUSTER-BASED PERMUTATION TEST ON TRACE ----
obsT = nan(1, nBins);
obsP = nan(1, nBins);

for b = 1:nBins
    c = cueTraceMat(:,b);
    p = pressTraceMat(:,b);
    keep = isfinite(c) & isfinite(p);
    c = c(keep);
    p = p(keep);

    if numel(c) >= 2
        [~, pval, ~, stats] = ttest(c, p);  % paired
        obsT(b) = stats.tstat;
        obsP(b) = pval;
    end
end

clusterMaskObs = isfinite(obsP) & (obsP < clusterAlpha);
obsClusters = findClusters1D_(clusterMaskObs);

obsClusterMass = nan(numel(obsClusters),1);
for iC = 1:numel(obsClusters)
    idx = obsClusters{iC};
    obsClusterMass(iC) = nansum(abs(obsT(idx)));
end

diffMat = pressTraceMat - cueTraceMat;  % paired difference, nSess x nBins
nSessTrace = size(diffMat,1);

maxClusterMassNull = zeros(Nperm_cluster,1);

rng(rngSeed);
for iPerm = 1:Nperm_cluster
    flips = ((rand(nSessTrace,1) > 0.5) * 2) - 1;  % +/-1 sign flip
    permMat = diffMat .* flips;

    permT = nan(1, nBins);
    permP = nan(1, nBins);

    for b = 1:nBins
        d = permMat(:,b);
        d = d(isfinite(d));

        if numel(d) >= 2
            [~, pval, ~, stats] = ttest(d, 0);
            permT(b) = stats.tstat;
            permP(b) = pval;
        end
    end

    permMask = isfinite(permP) & (permP < clusterAlpha);
    permClusters = findClusters1D_(permMask);

    if isempty(permClusters)
        maxClusterMassNull(iPerm) = 0;
    else
        permMasses = zeros(numel(permClusters),1);
        for jC = 1:numel(permClusters)
            idx = permClusters{jC};
            permMasses(jC) = nansum(abs(permT(idx)));
        end
        maxClusterMassNull(iPerm) = max(permMasses);
    end
end

obsClusterP = nan(numel(obsClusters),1);
for iC = 1:numel(obsClusters)
    obsClusterP(iC) = (1 + nnz(maxClusterMassNull >= obsClusterMass(iC))) / (Nperm_cluster + 1);
end

fprintf('\n=== PAIRED CLUSTER-BASED PERMUTATION TEST (CUE vs PRESS TRACE, LATE) ===\n');
fprintf('nLate=%d, clusterAlpha=%.3f, Nperm=%d\n', size(cueTraceMat,1), clusterAlpha, Nperm_cluster);
if isempty(obsClusters)
    fprintf('No observed supra-threshold clusters.\n');
else
    for iC = 1:numel(obsClusters)
        idx = obsClusters{iC};
        fprintf('Cluster %d: bins %d-%d | time %.3f to %.3f s | mass=%.6f | p=%.6g\n', ...
            iC, idx(1), idx(end), xPlot(idx(1)), xPlot(idx(end)), obsClusterMass(iC), obsClusterP(iC));
    end
end

%% ---- PLOT TRACE ----
figTrace = figure('Color','w','Position',figPosTrace);
ax2 = axes(figTrace); hold(ax2,'on');

patch(ax2, [xPlot, fliplr(xPlot)], [mCueTrace-seCueTrace, fliplr(mCueTrace+seCueTrace)], ...
    colCue, 'FaceAlpha', 0.15, 'EdgeColor','none');
patch(ax2, [xPlot, fliplr(xPlot)], [mPressTrace-sePressTrace, fliplr(mPressTrace+sePressTrace)], ...
    colPress, 'FaceAlpha', 0.15, 'EdgeColor','none');

plot(ax2, xPlot, mCueTrace,   '-', 'Color', colCue,   'LineWidth', traceLW);
plot(ax2, xPlot, mPressTrace, '-', 'Color', colPress, 'LineWidth', traceLW);

xline(ax2, cueLine,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
xline(ax2, pressLine, '--', 'Color', colPress, 'LineWidth', eventLineLW);
xline(ax2, lickLine,  '--', 'Color', colLick,  'LineWidth', eventLineLW);

xlabel(ax2, xlab, 'FontSize', labelFontSize);
ylabel(ax2, 'Accuracy', 'FontSize', labelFontSize);
title(ax2, 'Late sessions: cue vs press decoding across warped trial', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

legend(ax2, {'Cue decoder','Press decoder'}, 'Location','northeast', 'FontSize', tickFontSize);

set(ax2, 'FontSize', tickFontSize, ...
    'TickDir','out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box','off', ...
    'Layer','top');

if xUnit == "s" && forceIntegerSecondsTicks
    xl = xlim(ax2);
    ts = ceil(xl(1));
    te = floor(xl(2));
    if te >= ts
        xticks(ax2, ts:te);
    end
end

% draw significant clusters above the traces
yl = ylim(ax2);
yRange = yl(2) - yl(1);
yBase = yl(2) + 0.05*yRange;
yTop  = yl(2) + 0.12*yRange;

sigCount = 0;
for iC = 1:numel(obsClusters)
    if obsClusterP(iC) < 0.05
        idx = obsClusters{iC};
        x1 = xPlot(idx(1));
        x2 = xPlot(idx(end));
        sigCount = sigCount + 1;

        plot(ax2, [x1 x2], [yBase yBase], 'k-', 'LineWidth', axesTickLW, 'Clipping','off');
        text(ax2, mean([x1 x2]), yTop, '*', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', titleFontSize, 'FontWeight','bold', 'Color', [0 0 0], 'Clipping','off');
    end
end

if sigCount > 0
    ylim(ax2, [yl(1), yl(2) + 0.18*yRange]);
end

outTraceBase = fullfile(outDir, 'CUE_vs_PRESS_DECODING__TRACE__LateSessions');

set(figTrace, 'PaperPositionMode', 'auto');
exportgraphics(figTrace, [outTraceBase '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figTrace, [outTraceBase '.fig']);
try
    saveas(figTrace, [outTraceBase '.svg']);
catch
    print(figTrace, [outTraceBase '.svg'], '-dsvg');
end

fprintf('Saved:\n  %s.png\n  %s.fig\n  %s.svg\n', outTraceBase, outTraceBase, outTraceBase);

%% ---- PCA SECTION ----
% IMPORTANT:
% The saved MAT outputs used here contain session-level decoding summaries only.
% They do NOT contain trial-level neural features or class centroids needed to
% compute a 3D PCA plot with three cue points per session.
%
% To build that PCA faithfully, you would need, for each late session:
%   - trial x unit feature matrices (or class centroids) at a defined time/interval
%   - cue labels for those trials
% and then compute PCA on those neural features.
%
% So this section prints a note and exits cleanly instead of producing a
% misleading plot from insufficient data.

fprintf('\n=== PCA NOTE ===\n');
fprintf(['3D PCA cue-space plot was NOT generated.\n' ...
         'Reason: the saved decoder output MAT files do not contain trial-level or class-level\n' ...
         'neural feature coordinates needed to place 3 cue points per session in PCA space.\n' ...
         'You would need the underlying trial x unit matrices (or per-cue centroids) saved from\n' ...
         'the original decoding pipeline to make that figure faithfully.\n']);

%% ================= HELPERS =================
function plotPairedBar_(cueAll, pressAll, sessUse, ttl, ylab, outBase, ...
    figPos, titleFS, labelFS, tickFS, ...
    axLW, axTickLen, ms, edgeLW, jit, greyDot, ...
    colCueBar, colPressBar, colLine, barFaceAlpha, Nperm, rngSeed)

    cueVals   = cueAll(sessUse);
    pressVals = pressAll(sessUse);

    valid = isfinite(cueVals) & isfinite(pressVals);
    sessPlot = sessUse(valid);
    cueVals   = cueVals(valid);
    pressVals = pressVals(valid);

    assert(~isempty(sessPlot), 'No overlapping finite sessions found for plot: %s', ttl);

    fprintf('\n=== %s ===\n', ttl);
    fprintf('Using %d sessions\n', numel(sessPlot));
    fprintf('Sessions used: %s\n', mat2str(sessPlot));

    mCue   = mean(cueVals,   'omitnan');
    mPress = mean(pressVals, 'omitnan');

    semCue   = std(cueVals,   'omitnan') / sqrt(max(1, numel(cueVals)));
    semPress = std(pressVals, 'omitnan') / sqrt(max(1, numel(pressVals)));

    % paired sign-flip permutation test
    dObs = mean(pressVals - cueVals, 'omitnan');

    rng(rngSeed);
    dPerm = nan(Nperm,1);
    for i = 1:Nperm
        flips = (rand(numel(cueVals),1) > 0.5)*2 - 1;   % +/-1
        dPerm(i) = mean((pressVals(:) - cueVals(:)) .* flips, 'omitnan');
    end
    pPair = (1 + nnz(abs(dPerm) >= abs(dObs))) / (Nperm + 1);

    fprintf('Mean cue   = %.6f\n', mCue);
    fprintf('Mean press = %.6f\n', mPress);
    fprintf('Δ(Press-Cue)= %.6f\n', dObs);
    fprintf('Paired sign-flip permutation p = %.6g (Nperm=%d)\n', pPair, Nperm);

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    xPos = [1 2];
    barW = 0.60;

    % bars
    b = bar(ax, xPos, [mCue mPress], barW, ...
        'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axLW);
    b.CData(1,:) = colCueBar;
    b.CData(2,:) = colPressBar;
    b.FaceAlpha = barFaceAlpha;

    % error bars
    errorbar(ax, xPos, [mCue mPress], [semCue semPress], 'k', ...
        'LineStyle','none', 'LineWidth', axLW, 'CapSize', 18);

    % paired lines + dots
    rng(rngSeed);
    jitCue   = 1 + (rand(size(cueVals))-0.5)*2*jit;
    jitPress = 2 + (rand(size(pressVals))-0.5)*2*jit;

    for i = 1:numel(cueVals)
        plot(ax, [jitCue(i) jitPress(i)], [cueVals(i) pressVals(i)], '-', ...
            'Color', colLine, 'LineWidth', 1.5);
    end

    plot(ax, jitCue, cueVals, 'o', ...
        'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, ...
        'MarkerEdgeColor', greyDot, ...
        'LineWidth', edgeLW);

    plot(ax, jitPress, pressVals, 'o', ...
        'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, ...
        'MarkerEdgeColor', greyDot, ...
        'LineWidth', edgeLW);

    set(ax, 'XLim',[0.5 2.5], ...
        'XTick',[1 2], ...
        'XTickLabel',{'Cue decoder','Press decoder'});

    ylabel(ax, ylab, 'FontSize', labelFS);
    title(ax, ttl, 'FontWeight','bold', 'FontSize', titleFS);

    set(ax, 'FontSize', tickFS, ...
        'TickDir','out', ...
        'LineWidth', axLW, ...
        'TickLength', axTickLen, ...
        'Box','off', ...
        'Layer','top');

    % significance bracket
    addSigIfNeeded_(ax, 1, 2, pPair, [mCue mPress], axLW, tickFS);

    % save
    outPng = [outBase '.png'];
    outSvg = [outBase '.svg'];
    outFig = [outBase '.fig'];

    set(fig, 'PaperPositionMode', 'auto');
    exportgraphics(fig, outPng, 'Resolution', 300, 'BackgroundColor', 'white');
    savefig(fig, outFig);
    try
        saveas(fig, outSvg);
    catch
        print(fig, outSvg, '-dsvg');
    end

    fprintf('Saved:\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg);
end

function plotPairedScatter_(cueAll, pressAll, sessUse, ttl, xlab, ylab, outBase, ...
    figPos, titleFS, labelFS, tickFS, ...
    axLW, axTickLen, ms, edgeLW, ...
    colCueBar, colPressBar, greyDot)

    cueVals   = cueAll(sessUse);
    pressVals = pressAll(sessUse);

    valid = isfinite(cueVals) & isfinite(pressVals);
    sessPlot = sessUse(valid);
    cueVals   = cueVals(valid);
    pressVals = pressVals(valid);

    assert(~isempty(sessPlot), 'No overlapping finite sessions found for scatter: %s', ttl);

    fprintf('\n=== %s ===\n', ttl);
    fprintf('Using %d sessions\n', numel(sessPlot));
    fprintf('Sessions used: %s\n', mat2str(sessPlot));

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    allVals = [cueVals(:); pressVals(:)];
    lo = min(allVals);
    hi = max(allVals);
    pad = 0.05 * max(hi - lo, eps);

    xlimVals = [lo - pad, hi + pad];
    ylimVals = [lo - pad, hi + pad];

    % unity line
    plot(ax, xlimVals, xlimVals, '--', 'Color', greyDot, 'LineWidth', 2);

    % points
    scatter(ax, cueVals, pressVals, ms^2, ...
        'MarkerFaceColor', greyDot, ...
        'MarkerEdgeColor', greyDot, ...
        'LineWidth', edgeLW);

    % session labels
    for i = 1:numel(sessPlot)
        text(ax, cueVals(i) + 0.005*(xlimVals(2)-xlimVals(1)), ...
            pressVals(i), sprintf('%d', sessPlot(i)), ...
            'FontSize', max(10, round(0.45*tickFS)), ...
            'Color', [0 0 0], ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','middle');
    end

    set(ax, 'XLim', xlimVals, 'YLim', ylimVals);
    axis(ax, 'square');

    xlabel(ax, xlab, 'FontSize', labelFS);
    ylabel(ax, ylab, 'FontSize', labelFS);
    title(ax, ttl, 'FontWeight','bold', 'FontSize', titleFS);

    set(ax, 'FontSize', tickFS, ...
        'TickDir','out', ...
        'LineWidth', axLW, ...
        'TickLength', axTickLen, ...
        'Box','off', ...
        'Layer','top');

    outPng = [outBase '.png'];
    outSvg = [outBase '.svg'];
    outFig = [outBase '.fig'];

    set(fig, 'PaperPositionMode', 'auto');
    exportgraphics(fig, outPng, 'Resolution', 300, 'BackgroundColor', 'white');
    savefig(fig, outFig);
    try
        saveas(fig, outSvg);
    catch
        print(fig, outSvg, '-dsvg');
    end

    fprintf('Saved:\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg);
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
        yl = ylim(ax);
        yRange = yl(2) - yl(1);
        y = yMaxData + 0.08*yRange;
        h = 0.03*yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01*yRange, stars, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize', fs, ...
        'FontWeight','bold', ...
        'Color', [0 0 0]);
end

function clusters = findClusters1D_(mask)
    clusters = {};
    if isempty(mask) || ~any(mask), return; end

    d = diff([false, mask(:)', false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for i = 1:numel(starts)
        clusters{i} = starts(i):stops(i); %#ok<AGROW>
    end
end

function s = pToStars_(p)
    if p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    else
        s = '*';
    end
end