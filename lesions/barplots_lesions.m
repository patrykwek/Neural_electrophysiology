function plot_accuracy_per_100_trials_all_animals_with_weibull__DMS_FILES()
%% Plot accuracy per 100 subsequent trials for all animals on one plot
% and calculate Weibull fits
%
% INPUT FILES:
%   Control (black):
%     /Volumes/WD_BLACK/AB/AB/B1Bay_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat
%     /Volumes/WD_BLACK/AB/AB/F4Fig2_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/J3Jelly_ratBEHstruct_CUEDNAMES_MIN100.mat
%
%   DMS lesion (red):
%     /Volumes/WD_BLACK/AB/AB/J5Joy_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/J7Jasmine_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/L3Lychee_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/T8Truffle_ratBEHstruct_CUEDNAMES_MIN100.mat
%
% WHAT THIS SCRIPT DOES:
%   1) Loads the behavior struct files
%   2) Concatenates valid Hit values across sessions for each animal
%   3) Computes accuracy per every 100 subsequent trials
%   4) Limits all animals to the minimum number of 100-trial bins across files
%   5) Plots all animals on ONE plot
%        - Control animals in black
%        - DMS lesion animals in red
%   6) Adds a legend indicating Control vs DMS lesion
%   7) Calculates Weibull fits:
%        - separately for each animal
%        - separately for each group (Control and DMS lesion)
%      but DOES NOT plot the Weibull curves
%   8) Saves the figure and Weibull results in the same directory
%   9) At the end, makes a paired barplot of Early vs Late accuracy
%        - Early = first 1000 valid trials for each animal
%        - Late  = last 1000 valid trials taken from the COMMON plotted end,
%                  i.e. from the first (minNBins * BIN_SIZE) trials only
%  10) Also makes a delta plot:
%        - Delta = Late accuracy - Early accuracy
%        - compares Control vs DMS lesion deltas
%  11) Also makes a Control mean ± SEM and DMS lesion mean ± SEM plot
%  12) Also makes latency plots:
%        - cue-to-press latency Early vs Late barplot
%        - latency delta plot
%        - using the last 1000 valid latency trials from the common plotted end
%
% OUTPUTS SAVED IN:
%   /Volumes/WD_BLACK/AB/AB
%
% OUTPUT FILES:
%   accuracy_per_100trials_all_animals.png
%   accuracy_per_100trials_all_animals.svg
%   accuracy_per_100trials_weibull_results.mat
%   accuracy_per_100trials_weibull_results.txt
%   early_vs_late_accuracy_barplot.png
%   early_vs_late_accuracy_barplot.svg
%   early_vs_late_accuracy_barplot.fig
%   early_vs_late_accuracy_barplot_results.mat
%   delta_accuracy_plot.png
%   delta_accuracy_plot.svg
%   delta_accuracy_plot.fig
%   delta_accuracy_plot_results.mat
%   accuracy_per_100trials_group_mean_sem.png
%   accuracy_per_100trials_group_mean_sem.svg
%   accuracy_per_100trials_group_mean_sem.fig
%   accuracy_per_100trials_group_mean_sem_results.mat
%   latency_early_vs_late_barplot.png
%   latency_early_vs_late_barplot.svg
%   latency_early_vs_late_barplot.fig
%   latency_early_vs_late_barplot_results.mat
%   delta_latency_plot.png
%   delta_latency_plot.svg
%   delta_latency_plot.fig
%   delta_latency_plot_results.mat

clear; clc;

%% ---------------- USER SETTINGS ----------------
matFiles = { ...
    '/Volumes/WD_BLACK/AB/AB/B1Bay_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat', ...
    '/Volumes/WD_BLACK/AB/AB/F4Fig2_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/AB/AB/J3Jelly_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/AB/AB/J5Joy_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/AB/AB/J7Jasmine_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/AB/AB/L3Lychee_ratBEHstruct_CUEDNAMES_MIN100.mat', ...
    '/Volumes/WD_BLACK/AB/AB/T8Truffle_ratBEHstruct_CUEDNAMES_MIN100.mat'};

groupLabels = { ...
    'Control', ...
    'Control', ...
    'Control', ...
    'Control', ...
    'DMS lesion', ...
    'DMS lesion', ...
    'DMS lesion', ...
    'DMS lesion'};

outDir = '/Volumes/WD_BLACK/AB/AB';

outPng   = fullfile(outDir, 'accuracy_per_100trials_all_animals.png');
outSvg   = fullfile(outDir, 'accuracy_per_100trials_all_animals.svg');
outMat   = fullfile(outDir, 'accuracy_per_100trials_weibull_results.mat');
outTxt   = fullfile(outDir, 'accuracy_per_100trials_weibull_results.txt');

BIN_SIZE = 100;
EARLY_LATE_NTRIALS = 1000;

% Weibull bootstrap iterations for A CI
nBoot = 1000;
rng(0);

%% ---------------- STYLE ----------------
figPos = [120, 120, 980, 560];

dataLW      = 4;
fitLW       = 5;
eventLineLW = 5;

titleFontSize  = 24;
labelFontSize  = 24;
tickFontSize   = 20;
legendFontSize = 16;

axesTickLW = 3.0;

majorLenFrac      = 0.060;
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;
xLabelDownFrac    = 0.08;

ctrlColor = [0 0 0];   % black
lesColor  = [1 0 0];   % red

%% ---------------- LOAD ALL FILES / COMPUTE ACCURACY PER 100 TRIALS ----------------
nFiles = numel(matFiles);
animalData = struct([]);

for iFile = 1:nFiles
    matFile = matFiles{iFile};
    assert(exist(matFile,'file') == 2, 'File not found: %s', matFile);

    S   = load(matFile);
    beh = pickBehStruct_(S);
    assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in: %s', matFile);

    [~, baseName, ~] = fileparts(matFile);

    allHits = [];

    for sIdx = 1:numel(beh)
        sess = beh(sIdx);

        if ~isfield(sess, 'Hit') || isempty(sess.Hit)
            continue;
        end

        hitVals = sess.Hit;
        hitVals = double(hitVals(:));
        hitVals = hitVals(isfinite(hitVals));
        hitVals = hitVals(hitVals == 0 | hitVals == 1);

        if isempty(hitVals)
            continue;
        end

        allHits = [allHits; hitVals];
    end

    nTotalTrials = numel(allHits);
    nBins = floor(nTotalTrials / BIN_SIZE);

    pctCorrect = nan(nBins,1);
    nUsed      = zeros(nBins,1);

    for bIdx = 1:nBins
        idx1 = (bIdx - 1) * BIN_SIZE + 1;
        idx2 = bIdx * BIN_SIZE;
        binHits = allHits(idx1:idx2);

        nUsed(bIdx) = numel(binHits);
        pctCorrect(bIdx) = 100 * sum(binHits == 1) / nUsed(bIdx);
    end

    if strcmp(groupLabels{iFile}, 'Control')
        thisColor = ctrlColor;
    else
        thisColor = lesColor;
    end

    animalData(iFile).matFile       = matFile;
    animalData(iFile).baseName      = baseName;
    animalData(iFile).group         = groupLabels{iFile};
    animalData(iFile).color         = thisColor;
    animalData(iFile).allHits       = allHits;
    animalData(iFile).nTotalTrials  = nTotalTrials;
    animalData(iFile).nBins         = nBins;
    animalData(iFile).pctCorrect    = pctCorrect;
    animalData(iFile).nUsed         = nUsed;

    fprintf('%s: total valid Hit trials = %d, complete %d-trial bins = %d\n', ...
        baseName, nTotalTrials, BIN_SIZE, nBins);
end

%% ---------------- LIMIT TO MINIMUM BIN COUNT ----------------
allNBins = [animalData.nBins];
minNBins = min(allNBins);

fprintf('Minimum number of %d-trial bins across files: %d\n', BIN_SIZE, minNBins);

for iFile = 1:nFiles
    animalData(iFile).pctCorrect_use = animalData(iFile).pctCorrect(1:minNBins);
    animalData(iFile).nUsed_use      = animalData(iFile).nUsed(1:minNBins);
    animalData(iFile).x_use          = (1:minNBins)';
    animalData(iFile).allHits_use    = animalData(iFile).allHits(1:(minNBins * BIN_SIZE));
end

%% ---------------- WEIBULL MODEL ----------------
% y(b) = A - (A-B) * exp(-(b/lambda)^k)
% p = [B, A, lambda, k]
weibullFun = @(p, x) p(2) - (p(2)-p(1)) .* exp( - (max(x,eps)./p(3)).^p(4) );

LB = [0,   0,   0.1,  0.2];
UB = [100, 100, 1e3,  10];
opts = optimset('Display','off', 'MaxFunEvals', 5e4, 'MaxIter', 5e4);

%% ---------------- FIT EACH ANIMAL SEPARATELY ----------------
animalFits = repmat(initEmptyFitResult_(), 1, nFiles);

for iFile = 1:nFiles
    x = animalData(iFile).x_use;
    y = animalData(iFile).pctCorrect_use;
    w = double(animalData(iFile).nUsed_use);

    ok = isfinite(y) & isfinite(w) & (w > 0);
    x = x(ok);
    y = y(ok);
    w = w(ok);

    fitRes = initEmptyFitResult_();
    fitRes.name  = animalData(iFile).baseName;
    fitRes.group = animalData(iFile).group;
    fitRes.nSessionsUsed = numel(x);

    if numel(x) < 4
        fprintf('[WARN] Not enough bins to fit Weibull for %s\n', animalData(iFile).baseName);
        animalFits(iFile) = fitRes;
        continue;
    end

    [pHat, A_CI, asymSess] = fitWeibullWeighted_(x, y, w, weibullFun, LB, UB, opts, nBoot);

    fitRes.success   = true;
    fitRes.B         = pHat(1);
    fitRes.A         = pHat(2);
    fitRes.lambda    = pHat(3);
    fitRes.k         = pHat(4);
    fitRes.T10       = pHat(3) * (-log(1-0.10))^(1/pHat(4));
    fitRes.T50       = pHat(3) * (-log(1-0.50))^(1/pHat(4));
    fitRes.T90       = pHat(3) * (-log(1-0.90))^(1/pHat(4));
    fitRes.A_CI      = A_CI;
    fitRes.asymSess  = asymSess;

    animalFits(iFile) = fitRes;

    fprintf('\n=== ANIMAL FIT: %s ===\n', fitRes.name);
    fprintf('Group=%s\n', fitRes.group);
    fprintf('B=%.3f  A=%.3f  lambda=%.3f  k=%.3f\n', fitRes.B, fitRes.A, fitRes.lambda, fitRes.k);
    fprintf('T10=%.3f  T50=%.3f  T90=%.3f\n', fitRes.T10, fitRes.T50, fitRes.T90);
    fprintf('A 95%% CI = [%.3f, %.3f]\n', fitRes.A_CI(1), fitRes.A_CI(2));
    if isfinite(fitRes.asymSess)
        fprintf('Asymptote reached at bin %d\n', fitRes.asymSess);
    else
        fprintf('Asymptote not reached by criterion\n');
    end
end

%% ---------------- FIT GROUPS ----------------
groupNames = {'Control', 'DMS lesion'};
groupFits = repmat(initEmptyFitResult_(), 1, numel(groupNames));

for g = 1:numel(groupNames)
    grp = groupNames{g};
    idx = strcmp({animalData.group}, grp);

    Y = nan(sum(idx), minNBins);
    W = zeros(sum(idx), minNBins);

    tmp = animalData(idx);
    for j = 1:numel(tmp)
        Y(j,:) = tmp(j).pctCorrect_use(:)';
        W(j,:) = tmp(j).nUsed_use(:)';
    end

    groupPct = nan(minNBins,1);
    groupW   = zeros(minNBins,1);

    for b = 1:minNBins
        valid = isfinite(Y(:,b)) & (W(:,b) > 0);
        if ~any(valid)
            continue;
        end

        totalTrials  = sum(W(valid,b));
        totalCorrect = sum((Y(valid,b) ./ 100) .* W(valid,b));

        groupPct(b) = 100 * (totalCorrect / totalTrials);
        groupW(b)   = totalTrials;
    end

    x = (1:minNBins)';
    ok = isfinite(groupPct) & isfinite(groupW) & (groupW > 0);

    fitRes = initEmptyFitResult_();
    fitRes.name  = grp;
    fitRes.group = grp;
    fitRes.nSessionsUsed = sum(ok);
    fitRes.groupPct = groupPct;
    fitRes.groupW   = groupW;

    if sum(ok) >= 4
        [pHat, A_CI, asymSess] = fitWeibullWeighted_(x(ok), groupPct(ok), groupW(ok), weibullFun, LB, UB, opts, nBoot);

        fitRes.success   = true;
        fitRes.B         = pHat(1);
        fitRes.A         = pHat(2);
        fitRes.lambda    = pHat(3);
        fitRes.k         = pHat(4);
        fitRes.T10       = pHat(3) * (-log(1-0.10))^(1/pHat(4));
        fitRes.T50       = pHat(3) * (-log(1-0.50))^(1/pHat(4));
        fitRes.T90       = pHat(3) * (-log(1-0.90))^(1/pHat(4));
        fitRes.A_CI      = A_CI;
        fitRes.asymSess  = asymSess;
    else
        fprintf('[WARN] Not enough bins to fit Weibull for group %s\n', grp);
    end

    groupFits(g) = fitRes;

    fprintf('\n=== GROUP FIT: %s ===\n', fitRes.name);
    if fitRes.success
        fprintf('B=%.3f  A=%.3f  lambda=%.3f  k=%.3f\n', fitRes.B, fitRes.A, fitRes.lambda, fitRes.k);
        fprintf('T10=%.3f  T50=%.3f  T90=%.3f\n', fitRes.T10, fitRes.T50, fitRes.T90);
        fprintf('A 95%% CI = [%.3f, %.3f]\n', fitRes.A_CI(1), fitRes.A_CI(2));
        if isfinite(fitRes.asymSess)
            fprintf('Asymptote reached at bin %d\n', fitRes.asymSess);
        else
            fprintf('Asymptote not reached by criterion\n');
        end
    end
end

%% ---------------- PLOT ALL ANIMALS ON ONE FIGURE ----------------
fig = figure('Color','w', 'Position', figPos);
ax = axes(fig); hold(ax,'on');

yl = [0 100];
set(ax, 'XLim', [1 minNBins], 'YLim', yl);

for iFile = 1:nFiles
    x = animalData(iFile).x_use;
    y = animalData(iFile).pctCorrect_use;

    plot(ax, x, y, '-', ...
        'Color', animalData(iFile).color, ...
        'LineWidth', dataLW);
end

hxlab = xlabel(ax, sprintf('Consecutive %d-trial bin', BIN_SIZE), 'FontSize', labelFontSize);
ylabel(ax, 'Accuracy (% correct)', 'FontSize', labelFontSize);
title(ax, sprintf('Accuracy per %d subsequent trials for each animal', BIN_SIZE), ...
    'FontWeight', 'bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box', 'off', 'XLim', [1 minNBins], 'YLim', [0 100]);
set(ax, 'TickDir', 'out', 'LineWidth', axesTickLW, 'Layer', 'top');

%% ---- TICKS: label every 10th; draw ticks + labels manually ----
ax.XTick = 1:minNBins;
labIdx = 10:10:minNBins;

ax.XTickLabel = repmat({''}, 1, minNBins);
xtickangle(ax, 0);
set(ax, 'TickLength', [0 0]);

y0 = ax.YLim(1);
yr = diff(ax.YLim);

majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

tickCol = [0 0 0];

for s = 1:minNBins
    if ismember(s, labIdx)
        L  = majorLen;
        lw = axesTickLW;
    else
        L  = minorLen;
        lw = max(1.5, axesTickLW * 0.55);
    end
    line(ax, [s s], [y0, y0 - L], 'Color', tickCol, 'LineWidth', lw, 'Clipping', 'off');
end

yTickText = y0 - (majorLen + tickLabelDownFrac * yr);
for s = labIdx
    text(ax, s, yTickText, num2str(s), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickFontSize, ...
        'Color', [0 0 0], ...
        'Clipping', 'off');
end

xCenter = mean(ax.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac * yr);
set(hxlab, 'Units', 'data', 'Position', [xCenter, yXlab, 0]);

ax.Position = [0.12 0.22 0.72 0.70];

ctrlHandle = plot(ax, nan, nan, '-', 'Color', ctrlColor, 'LineWidth', dataLW);
lesHandle  = plot(ax, nan, nan, '-', 'Color', lesColor,  'LineWidth', dataLW);

legend(ax, [ctrlHandle, lesHandle], {'Control', 'DMS lesion'}, ...
    'Location', 'eastoutside', ...
    'Box', 'off', ...
    'FontSize', legendFontSize);

grid(ax, 'off');
set(ax, 'XLim', [1 minNBins], 'YLim', [0 100]);

saveas(fig, outPng);
saveas(fig, outSvg);

fprintf('\nSaved figure:\n%s\n%s\n', outPng, outSvg);

%% ---------------- GROUP MEAN ± SEM PLOT ----------------
groupMeanOutBase = fullfile(outDir, 'accuracy_per_100trials_group_mean_sem');
groupMeanOutPng  = [groupMeanOutBase '.png'];
groupMeanOutSvg  = [groupMeanOutBase '.svg'];
groupMeanOutFig  = [groupMeanOutBase '.fig'];
groupMeanOutMat  = [groupMeanOutBase '_results.mat'];

ctrlIdx = strcmp({animalData.group}, 'Control');
lesIdx  = strcmp({animalData.group}, 'DMS lesion');

ctrlMat = nan(sum(ctrlIdx), minNBins);
lesMat  = nan(sum(lesIdx), minNBins);

tmpCtrl = animalData(ctrlIdx);
for i = 1:numel(tmpCtrl)
    ctrlMat(i,:) = tmpCtrl(i).pctCorrect_use(:)';
end

tmpLes = animalData(lesIdx);
for i = 1:numel(tmpLes)
    lesMat(i,:) = tmpLes(i).pctCorrect_use(:)';
end

ctrlMean = mean(ctrlMat, 1, 'omitnan');
lesMean  = mean(lesMat, 1, 'omitnan');

ctrlSEM = std(ctrlMat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(ctrlMat), 1)));
lesSEM  = std(lesMat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(lesMat), 1)));

xPlot = 1:minNBins;

figMean = figure('Color','w', 'Position', figPos);
axMean = axes(figMean); hold(axMean,'on');

yl = [0 100];
set(axMean, 'XLim', [1 minNBins], 'YLim', yl);

patch(axMean, [xPlot fliplr(xPlot)], [ctrlMean-ctrlSEM fliplr(ctrlMean+ctrlSEM)], ...
    ctrlColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

patch(axMean, [xPlot fliplr(xPlot)], [lesMean-lesSEM fliplr(lesMean+lesSEM)], ...
    lesColor, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

plot(axMean, xPlot, ctrlMean, '-', 'Color', ctrlColor, 'LineWidth', dataLW);
plot(axMean, xPlot, lesMean,  '-', 'Color', lesColor,  'LineWidth', dataLW);

hxlabMean = xlabel(axMean, sprintf('Consecutive %d-trial bin', BIN_SIZE), 'FontSize', labelFontSize);
ylabel(axMean, 'Accuracy (% correct)', 'FontSize', labelFontSize);
title(axMean, sprintf('Accuracy per %d subsequent trials: group mean \\pm SEM', BIN_SIZE), ...
    'FontWeight', 'bold', 'FontSize', titleFontSize);

set(axMean, 'FontSize', tickFontSize, 'Box', 'off', 'XLim', [1 minNBins], 'YLim', [0 100]);
set(axMean, 'TickDir', 'out', 'LineWidth', axesTickLW, 'Layer', 'top');

axMean.XTick = 1:minNBins;
labIdx = 10:10:minNBins;

axMean.XTickLabel = repmat({''}, 1, minNBins);
xtickangle(axMean, 0);
set(axMean, 'TickLength', [0 0]);

y0 = axMean.YLim(1);
yr = diff(axMean.YLim);

majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

for s = 1:minNBins
    if ismember(s, labIdx)
        L  = majorLen;
        lw = axesTickLW;
    else
        L  = minorLen;
        lw = max(1.5, axesTickLW * 0.55);
    end
    line(axMean, [s s], [y0, y0 - L], 'Color', tickCol, 'LineWidth', lw, 'Clipping', 'off');
end

yTickText = y0 - (majorLen + tickLabelDownFrac * yr);
for s = labIdx
    text(axMean, s, yTickText, num2str(s), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickFontSize, ...
        'Color', [0 0 0], ...
        'Clipping', 'off');
end

xCenter = mean(axMean.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac * yr);
set(hxlabMean, 'Units', 'data', 'Position', [xCenter, yXlab, 0]);

axMean.Position = [0.12 0.22 0.72 0.70];

ctrlHandleMean = plot(axMean, nan, nan, '-', 'Color', ctrlColor, 'LineWidth', dataLW);
lesHandleMean  = plot(axMean, nan, nan, '-', 'Color', lesColor,  'LineWidth', dataLW);

legend(axMean, [ctrlHandleMean, lesHandleMean], {'Control mean \pm SEM', 'DMS lesion mean \pm SEM'}, ...
    'Location', 'eastoutside', ...
    'Box', 'off', ...
    'FontSize', legendFontSize);

grid(axMean, 'off');
set(axMean, 'XLim', [1 minNBins], 'YLim', [0 100]);

groupMeanResults = struct();
groupMeanResults.matFiles    = matFiles;
groupMeanResults.groupLabels = groupLabels;
groupMeanResults.binSize     = BIN_SIZE;
groupMeanResults.minNBins    = minNBins;
groupMeanResults.ctrlMat     = ctrlMat;
groupMeanResults.lesMat      = lesMat;
groupMeanResults.ctrlMean    = ctrlMean;
groupMeanResults.lesMean     = lesMean;
groupMeanResults.ctrlSEM     = ctrlSEM;
groupMeanResults.lesSEM      = lesSEM;

save(groupMeanOutMat, 'groupMeanResults', '-v7.3');

set(figMean, 'PaperPositionMode', 'auto');
exportgraphics(figMean, groupMeanOutPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figMean, groupMeanOutFig);
try
    saveas(figMean, groupMeanOutSvg);
catch
    print(figMean, groupMeanOutSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', groupMeanOutPng, groupMeanOutFig, groupMeanOutSvg, groupMeanOutMat);

%% ---------------- SAVE RESULTS ----------------
results = struct();
results.matFiles    = matFiles;
results.groupLabels = groupLabels;
results.binSize     = BIN_SIZE;
results.minNBins    = minNBins;
results.animalData  = animalData;
results.animalFits  = animalFits;
results.groupFits   = groupFits;

save(outMat, 'results', '-v7.3');
fprintf('Saved results MAT: %s\n', outMat);

fid = fopen(outTxt, 'w');
assert(fid ~= -1, 'Could not open text output for writing: %s', outTxt);

fprintf(fid, 'Accuracy per %d-trial bin Weibull results\n', BIN_SIZE);
fprintf(fid, '=========================================\n\n');
fprintf(fid, 'Minimum number of %d-trial bins used across all files: %d\n\n', BIN_SIZE, minNBins);

fprintf(fid, 'ANIMAL FITS\n');
fprintf(fid, '-----------\n');
for iFile = 1:numel(animalFits)
    fr = animalFits(iFile);
    fprintf(fid, '%s\n', fr.name);
    fprintf(fid, '  Group: %s\n', fr.group);
    fprintf(fid, '  Success: %d\n', fr.success);
    fprintf(fid, '  Bins used: %d\n', fr.nSessionsUsed);
    if fr.success
        fprintf(fid, '  B=%.6f\n', fr.B);
        fprintf(fid, '  A=%.6f\n', fr.A);
        fprintf(fid, '  lambda=%.6f\n', fr.lambda);
        fprintf(fid, '  k=%.6f\n', fr.k);
        fprintf(fid, '  T10=%.6f\n', fr.T10);
        fprintf(fid, '  T50=%.6f\n', fr.T50);
        fprintf(fid, '  T90=%.6f\n', fr.T90);
        fprintf(fid, '  A_CI_low=%.6f\n', fr.A_CI(1));
        fprintf(fid, '  A_CI_high=%.6f\n', fr.A_CI(2));
        fprintf(fid, '  asymBin=%g\n', fr.asymSess);
    end
    fprintf(fid, '\n');
end

fprintf(fid, 'GROUP FITS\n');
fprintf(fid, '----------\n');
for g = 1:numel(groupFits)
    fr = groupFits(g);
    fprintf(fid, '%s\n', fr.name);
    fprintf(fid, '  Group: %s\n', fr.group);
    fprintf(fid, '  Success: %d\n', fr.success);
    fprintf(fid, '  Bins used: %d\n', fr.nSessionsUsed);
    if fr.success
        fprintf(fid, '  B=%.6f\n', fr.B);
        fprintf(fid, '  A=%.6f\n', fr.A);
        fprintf(fid, '  lambda=%.6f\n', fr.lambda);
        fprintf(fid, '  k=%.6f\n', fr.k);
        fprintf(fid, '  T10=%.6f\n', fr.T10);
        fprintf(fid, '  T50=%.6f\n', fr.T50);
        fprintf(fid, '  T90=%.6f\n', fr.T90);
        fprintf(fid, '  A_CI_low=%.6f\n', fr.A_CI(1));
        fprintf(fid, '  A_CI_high=%.6f\n', fr.A_CI(2));
        fprintf(fid, '  asymBin=%g\n', fr.asymSess);
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Saved results TXT: %s\n', outTxt);

%% ---------------- EARLY vs LATE BARPLOT ----------------
barOutBase = fullfile(outDir, 'early_vs_late_accuracy_barplot');
barOutPng  = [barOutBase '.png'];
barOutSvg  = [barOutBase '.svg'];
barOutFig  = [barOutBase '.fig'];
barOutMat  = [barOutBase '_results.mat'];

figPosBar = [120, 120, 640, 760];

barTitleFontSize = 30;
barLabelFontSize = 30;
barTickFontSize  = 26;

barAxesTickLW  = 4.0;
barAxesTickLen = [0.02 0.02];

scatterMS     = 10;
scatterEdgeLW = 1.0;
jit           = 0.06;

colCtrlEarly = [0.75 0.75 0.75];
colCtrlLate  = [0.35 0.35 0.35];
colLesEarly  = [1.00 0.70 0.70];
colLesLate   = [1.00 0.00 0.00];

colLine      = [0.45 0.45 0.45];
greyDot      = [0.60 0.60 0.60];
barFaceAlpha = 0.25;

for iFile = 1:nFiles
    nUseTrials = numel(animalData(iFile).allHits_use);
    assert(nUseTrials >= EARLY_LATE_NTRIALS, ...
        'Animal %s has only %d plotted trials; need at least %d.', ...
        animalData(iFile).baseName, nUseTrials, EARLY_LATE_NTRIALS);

    earlyHits = animalData(iFile).allHits_use(1:EARLY_LATE_NTRIALS);
    lateHits  = animalData(iFile).allHits_use(end-EARLY_LATE_NTRIALS+1:end);

    animalData(iFile).earlyAcc_1000 = 100 * sum(earlyHits == 1) / EARLY_LATE_NTRIALS;
    animalData(iFile).lateAcc_1000  = 100 * sum(lateHits  == 1) / EARLY_LATE_NTRIALS;
    animalData(iFile).deltaAcc_1000 = animalData(iFile).lateAcc_1000 - animalData(iFile).earlyAcc_1000;
end

isCtrl = strcmp({animalData.group}, 'Control');
isLes  = strcmp({animalData.group}, 'DMS lesion');

ctrlEarly = [animalData(isCtrl).earlyAcc_1000]';
ctrlLate  = [animalData(isCtrl).lateAcc_1000]';
ctrlDelta = [animalData(isCtrl).deltaAcc_1000]';

lesEarly  = [animalData(isLes).earlyAcc_1000]';
lesLate   = [animalData(isLes).lateAcc_1000]';
lesDelta  = [animalData(isLes).deltaAcc_1000]';

mCtrlEarly = mean(ctrlEarly, 'omitnan');
mCtrlLate  = mean(ctrlLate,  'omitnan');
mLesEarly  = mean(lesEarly,  'omitnan');
mLesLate   = mean(lesLate,   'omitnan');

semCtrlEarly = std(ctrlEarly, 'omitnan') / sqrt(max(1, numel(ctrlEarly)));
semCtrlLate  = std(ctrlLate,  'omitnan') / sqrt(max(1, numel(ctrlLate)));
semLesEarly  = std(lesEarly,  'omitnan') / sqrt(max(1, numel(lesEarly)));
semLesLate   = std(lesLate,   'omitnan') / sqrt(max(1, numel(lesLate)));

%% ---------------- STATISTICAL TESTS ----------------
[pCtrl_signrank, ~, statsCtrl_signrank] = signrank(ctrlEarly, ctrlLate);
[pLes_signrank,  ~, statsLes_signrank]  = signrank(lesEarly, lesLate);
[pDelta_ranksum, ~, statsDelta_ranksum] = ranksum(ctrlDelta, lesDelta);

fprintf('\n=== EARLY vs LATE SUMMARY (USING COMMON PLOTTED END) ===\n');
fprintf('Control early mean = %.6f | late mean = %.6f | Delta(late-early) = %.6f\n', ...
    mCtrlEarly, mCtrlLate, mean(ctrlDelta, 'omitnan'));
fprintf('DMS lesion early mean = %.6f | late mean = %.6f | Delta(late-early) = %.6f\n', ...
    mLesEarly, mLesLate, mean(lesDelta, 'omitnan'));

fprintf('\n=== STATISTICAL TESTS ===\n');
fprintf('Control Early vs Late: Wilcoxon signed-rank p = %.6g', pCtrl_signrank);
if isfield(statsCtrl_signrank, 'signedrank')
    fprintf(' | signed-rank = %g\n', statsCtrl_signrank.signedrank);
else
    fprintf('\n');
end

fprintf('DMS lesion Early vs Late: Wilcoxon signed-rank p = %.6g', pLes_signrank);
if isfield(statsLes_signrank, 'signedrank')
    fprintf(' | signed-rank = %g\n', statsLes_signrank.signedrank);
else
    fprintf('\n');
end

fprintf('Control vs DMS lesion Delta: Wilcoxon rank-sum p = %.6g', pDelta_ranksum);
if isfield(statsDelta_ranksum, 'ranksum')
    fprintf(' | rank-sum = %g\n', statsDelta_ranksum.ranksum);
else
    fprintf('\n');
end

fprintf('\nControl deltas:\n');
disp(ctrlDelta(:)');
fprintf('DMS lesion deltas:\n');
disp(lesDelta(:)');

figBar = figure('Color', 'w', 'Position', figPosBar);
axBar = axes(figBar); hold(axBar, 'on');

xPos = [1 2 4 5];
barW = 0.60;

barVals = [mCtrlEarly, mCtrlLate, mLesEarly, mLesLate];
semVals = [semCtrlEarly, semCtrlLate, semLesEarly, semLesLate];

b = bar(axBar, xPos, barVals, barW, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', barAxesTickLW);

b.CData(1,:) = colCtrlEarly;
b.CData(2,:) = colCtrlLate;
b.CData(3,:) = colLesEarly;
b.CData(4,:) = colLesLate;
b.FaceAlpha = barFaceAlpha;

errorbar(axBar, xPos, barVals, semVals, 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', barAxesTickLW, ...
    'CapSize', 18);

rng(0);
jitCtrlEarly = 1 + (rand(size(ctrlEarly)) - 0.5) * 2 * jit;
jitCtrlLate  = 2 + (rand(size(ctrlLate))  - 0.5) * 2 * jit;

for i = 1:numel(ctrlEarly)
    plot(axBar, [jitCtrlEarly(i) jitCtrlLate(i)], [ctrlEarly(i) ctrlLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(axBar, jitCtrlEarly, ctrlEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axBar, jitCtrlLate, ctrlLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

rng(1);
jitLesEarly = 4 + (rand(size(lesEarly)) - 0.5) * 2 * jit;
jitLesLate  = 5 + (rand(size(lesLate))  - 0.5) * 2 * jit;

for i = 1:numel(lesEarly)
    plot(axBar, [jitLesEarly(i) jitLesLate(i)], [lesEarly(i) lesLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(axBar, jitLesEarly, lesEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axBar, jitLesLate, lesLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

set(axBar, 'XLim', [0.3 5.7], ...
    'XTick', xPos, ...
    'XTickLabel', {'Control Early', 'Control Late', 'DMS lesion Early', 'DMS lesion Late'}, ...
    'YLim', [0 100]);

ylabel(axBar, 'Accuracy (% correct)', 'FontSize', barLabelFontSize);
title(axBar, sprintf('Early vs Late accuracy (%d trials each)', EARLY_LATE_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', barTitleFontSize);

set(axBar, 'FontSize', barTickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', barAxesTickLW, ...
    'TickLength', barAxesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

xtickangle(axBar, 45);

addSigIfNeeded_(axBar, 1, 2, pCtrl_signrank, [mCtrlEarly mCtrlLate], barAxesTickLW, barTickFontSize);
addSigIfNeeded_(axBar, 4, 5, pLes_signrank,  [mLesEarly  mLesLate],  barAxesTickLW, barTickFontSize);

barResults = struct();
barResults.matFiles           = matFiles;
barResults.groupLabels        = groupLabels;
barResults.nTrials            = EARLY_LATE_NTRIALS;
barResults.binSize            = BIN_SIZE;
barResults.minNBins           = minNBins;
barResults.commonEndTrials    = minNBins * BIN_SIZE;
barResults.animalData         = animalData;
barResults.ctrlEarly          = ctrlEarly;
barResults.ctrlLate           = ctrlLate;
barResults.ctrlDelta          = ctrlDelta;
barResults.lesEarly           = lesEarly;
barResults.lesLate            = lesLate;
barResults.lesDelta           = lesDelta;
barResults.meanVals           = barVals;
barResults.semVals            = semVals;
barResults.pCtrl_signrank     = pCtrl_signrank;
barResults.pLes_signrank      = pLes_signrank;
barResults.pDelta_ranksum     = pDelta_ranksum;
barResults.statsCtrl_signrank = statsCtrl_signrank;
barResults.statsLes_signrank  = statsLes_signrank;
barResults.statsDelta_ranksum = statsDelta_ranksum;

save(barOutMat, 'barResults', '-v7.3');

set(figBar, 'PaperPositionMode', 'auto');
exportgraphics(figBar, barOutPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figBar, barOutFig);
try
    saveas(figBar, barOutSvg);
catch
    print(figBar, barOutSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', barOutPng, barOutFig, barOutSvg, barOutMat);

%% ---------------- DELTA PLOT ----------------
deltaOutBase = fullfile(outDir, 'delta_accuracy_plot');
deltaOutPng  = [deltaOutBase '.png'];
deltaOutSvg  = [deltaOutBase '.svg'];
deltaOutFig  = [deltaOutBase '.fig'];
deltaOutMat  = [deltaOutBase '_results.mat'];

figDelta = figure('Color', 'w', 'Position', [120, 120, 520, 700]);
axDelta = axes(figDelta); hold(axDelta, 'on');

mCtrlDelta = mean(ctrlDelta, 'omitnan');
mLesDelta  = mean(lesDelta,  'omitnan');

semCtrlDelta = std(ctrlDelta, 'omitnan') / sqrt(max(1, numel(ctrlDelta)));
semLesDelta  = std(lesDelta,  'omitnan') / sqrt(max(1, numel(lesDelta)));

xDelta = [1 2];
deltaBarW = 0.60;

bd = bar(axDelta, xDelta, [mCtrlDelta mLesDelta], deltaBarW, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', barAxesTickLW);
bd.CData(1,:) = [0.45 0.45 0.45];
bd.CData(2,:) = lesColor;
bd.FaceAlpha = barFaceAlpha;

errorbar(axDelta, xDelta, [mCtrlDelta mLesDelta], [semCtrlDelta semLesDelta], 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', barAxesTickLW, ...
    'CapSize', 18);

rng(2);
jitCtrlDelta = 1 + (rand(size(ctrlDelta)) - 0.5) * 2 * jit;
jitLesDelta  = 2 + (rand(size(lesDelta))  - 0.5) * 2 * jit;

plot(axDelta, jitCtrlDelta, ctrlDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axDelta, jitLesDelta, lesDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

yline(axDelta, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

set(axDelta, 'XLim', [0.5 2.5], ...
    'XTick', [1 2], ...
    'XTickLabel', {'Control', 'DMS lesion'});

ylabel(axDelta, '\Delta Accuracy (Late - Early)', 'FontSize', barLabelFontSize);
title(axDelta, sprintf('Change from Early to Late (%d trials each)', EARLY_LATE_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', barTitleFontSize);

set(axDelta, 'FontSize', barTickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', barAxesTickLW, ...
    'TickLength', barAxesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

addSigIfNeeded_(axDelta, 1, 2, pDelta_ranksum, [mCtrlDelta mLesDelta], barAxesTickLW, barTickFontSize);

deltaResults = struct();
deltaResults.matFiles           = matFiles;
deltaResults.groupLabels        = groupLabels;
deltaResults.nTrials            = EARLY_LATE_NTRIALS;
deltaResults.binSize            = BIN_SIZE;
deltaResults.minNBins           = minNBins;
deltaResults.commonEndTrials    = minNBins * BIN_SIZE;
deltaResults.ctrlDelta          = ctrlDelta;
deltaResults.lesDelta           = lesDelta;
deltaResults.meanCtrlDelta      = mCtrlDelta;
deltaResults.meanLesDelta       = mLesDelta;
deltaResults.semCtrlDelta       = semCtrlDelta;
deltaResults.semLesDelta        = semLesDelta;
deltaResults.pDelta_ranksum     = pDelta_ranksum;
deltaResults.statsDelta_ranksum = statsDelta_ranksum;

save(deltaOutMat, 'deltaResults', '-v7.3');

set(figDelta, 'PaperPositionMode', 'auto');
exportgraphics(figDelta, deltaOutPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figDelta, deltaOutFig);
try
    saveas(figDelta, deltaOutSvg);
catch
    print(figDelta, deltaOutSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', deltaOutPng, deltaOutFig, deltaOutSvg, deltaOutMat);

%% ---------------- LATENCY ANALYSIS ----------------
latBarOutBase = fullfile(outDir, 'latency_early_vs_late_barplot');
latBarOutPng  = [latBarOutBase '.png'];
latBarOutSvg  = [latBarOutBase '.svg'];
latBarOutFig  = [latBarOutBase '.fig'];
latBarOutMat  = [latBarOutBase '_results.mat'];

latDeltaOutBase = fullfile(outDir, 'delta_latency_plot');
latDeltaOutPng  = [latDeltaOutBase '.png'];
latDeltaOutSvg  = [latDeltaOutBase '.svg'];
latDeltaOutFig  = [latDeltaOutBase '.fig'];
latDeltaOutMat  = [latDeltaOutBase '_results.mat'];

for iFile = 1:nFiles
    matFile = matFiles{iFile};
    S   = load(matFile);
    beh = pickBehStruct_(S);

    allLatencies = [];

    for sIdx = 1:numel(beh)
        latVals = extractLatenciesFromSession_(beh(sIdx));
        if isempty(latVals)
            continue;
        end
        allLatencies = [allLatencies; latVals(:)];
    end

    animalData(iFile).allLatencies = allLatencies;
    animalData(iFile).nLatencyTrials = numel(allLatencies);

    fprintf('%s: total valid latency trials = %d\n', animalData(iFile).baseName, numel(allLatencies));
end

allNLatencyTrials = [animalData.nLatencyTrials];
commonEndLatencyTrials = min(allNLatencyTrials);

fprintf('Common plotted end for latency across files: %d valid trials\n', commonEndLatencyTrials);

for iFile = 1:nFiles
    nUseTrials = numel(animalData(iFile).allLatencies);
    assert(nUseTrials >= EARLY_LATE_NTRIALS, ...
        'Animal %s has only %d valid latency trials; need at least %d.', ...
        animalData(iFile).baseName, nUseTrials, EARLY_LATE_NTRIALS);

    animalData(iFile).allLatencies_use = animalData(iFile).allLatencies(1:commonEndLatencyTrials);

    earlyLat = animalData(iFile).allLatencies_use(1:EARLY_LATE_NTRIALS);
    lateLat  = animalData(iFile).allLatencies_use(end-EARLY_LATE_NTRIALS+1:end);

    animalData(iFile).earlyLat_1000 = mean(earlyLat, 'omitnan');
    animalData(iFile).lateLat_1000  = mean(lateLat,  'omitnan');
    animalData(iFile).deltaLat_1000 = animalData(iFile).lateLat_1000 - animalData(iFile).earlyLat_1000;
end

ctrlLatEarly = [animalData(isCtrl).earlyLat_1000]';
ctrlLatLate  = [animalData(isCtrl).lateLat_1000]';
ctrlLatDelta = [animalData(isCtrl).deltaLat_1000]';

lesLatEarly  = [animalData(isLes).earlyLat_1000]';
lesLatLate   = [animalData(isLes).lateLat_1000]';
lesLatDelta  = [animalData(isLes).deltaLat_1000]';

mCtrlLatEarly = mean(ctrlLatEarly, 'omitnan');
mCtrlLatLate  = mean(ctrlLatLate,  'omitnan');
mLesLatEarly  = mean(lesLatEarly,  'omitnan');
mLesLatLate   = mean(lesLatLate,   'omitnan');

semCtrlLatEarly = std(ctrlLatEarly, 'omitnan') / sqrt(max(1, numel(ctrlLatEarly)));
semCtrlLatLate  = std(ctrlLatLate,  'omitnan') / sqrt(max(1, numel(ctrlLatLate)));
semLesLatEarly  = std(lesLatEarly,  'omitnan') / sqrt(max(1, numel(lesLatEarly)));
semLesLatLate   = std(lesLatLate,   'omitnan') / sqrt(max(1, numel(lesLatLate)));

[pCtrlLat_signrank, ~, statsCtrlLat_signrank] = signrank(ctrlLatEarly, ctrlLatLate);
[pLesLat_signrank,  ~, statsLesLat_signrank]  = signrank(lesLatEarly, lesLatLate);
[pLatDelta_ranksum, ~, statsLatDelta_ranksum] = ranksum(ctrlLatDelta, lesLatDelta);

fprintf('\n=== LATENCY EARLY vs LATE SUMMARY (USING COMMON PLOTTED END) ===\n');
fprintf('Control early mean = %.6f ms | late mean = %.6f ms | Delta(late-early) = %.6f ms\n', ...
    mCtrlLatEarly, mCtrlLatLate, mean(ctrlLatDelta, 'omitnan'));
fprintf('DMS lesion early mean = %.6f ms | late mean = %.6f ms | Delta(late-early) = %.6f ms\n', ...
    mLesLatEarly, mLesLatLate, mean(lesLatDelta, 'omitnan'));

fprintf('\n=== LATENCY STATISTICAL TESTS ===\n');
fprintf('Control Early vs Late latency: Wilcoxon signed-rank p = %.6g', pCtrlLat_signrank);
if isfield(statsCtrlLat_signrank, 'signedrank')
    fprintf(' | signed-rank = %g\n', statsCtrlLat_signrank.signedrank);
else
    fprintf('\n');
end

fprintf('DMS lesion Early vs Late latency: Wilcoxon signed-rank p = %.6g', pLesLat_signrank);
if isfield(statsLesLat_signrank, 'signedrank')
    fprintf(' | signed-rank = %g\n', statsLesLat_signrank.signedrank);
else
    fprintf('\n');
end

fprintf('Control vs DMS lesion latency Delta: Wilcoxon rank-sum p = %.6g', pLatDelta_ranksum);
if isfield(statsLatDelta_ranksum, 'ranksum')
    fprintf(' | rank-sum = %g\n', statsLatDelta_ranksum.ranksum);
else
    fprintf('\n');
end

fprintf('\nControl latency deltas:\n');
disp(ctrlLatDelta(:)');
fprintf('DMS lesion latency deltas:\n');
disp(lesLatDelta(:)');

figLatBar = figure('Color', 'w', 'Position', figPosBar);
axLatBar = axes(figLatBar); hold(axLatBar, 'on');

latBarVals = [mCtrlLatEarly, mCtrlLatLate, mLesLatEarly, mLesLatLate];
latSemVals = [semCtrlLatEarly, semCtrlLatLate, semLesLatEarly, semLesLatLate];

bLat = bar(axLatBar, xPos, latBarVals, barW, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', barAxesTickLW);

bLat.CData(1,:) = colCtrlEarly;
bLat.CData(2,:) = colCtrlLate;
bLat.CData(3,:) = colLesEarly;
bLat.CData(4,:) = colLesLate;
bLat.FaceAlpha = barFaceAlpha;

errorbar(axLatBar, xPos, latBarVals, latSemVals, 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', barAxesTickLW, ...
    'CapSize', 18);

rng(10);
jitCtrlLatEarly = 1 + (rand(size(ctrlLatEarly)) - 0.5) * 2 * jit;
jitCtrlLatLate  = 2 + (rand(size(ctrlLatLate))  - 0.5) * 2 * jit;

for i = 1:numel(ctrlLatEarly)
    plot(axLatBar, [jitCtrlLatEarly(i) jitCtrlLatLate(i)], [ctrlLatEarly(i) ctrlLatLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(axLatBar, jitCtrlLatEarly, ctrlLatEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axLatBar, jitCtrlLatLate, ctrlLatLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

rng(11);
jitLesLatEarly = 4 + (rand(size(lesLatEarly)) - 0.5) * 2 * jit;
jitLesLatLate  = 5 + (rand(size(lesLatLate))  - 0.5) * 2 * jit;

for i = 1:numel(lesLatEarly)
    plot(axLatBar, [jitLesLatEarly(i) jitLesLatLate(i)], [lesLatEarly(i) lesLatLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(axLatBar, jitLesLatEarly, lesLatEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axLatBar, jitLesLatLate, lesLatLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

yMaxLat = max([latBarVals + latSemVals, ctrlLatEarly', ctrlLatLate', lesLatEarly', lesLatLate']);
if ~isfinite(yMaxLat) || isempty(yMaxLat)
    yMaxLat = 1;
end

set(axLatBar, 'XLim', [0.3 5.7], ...
    'XTick', xPos, ...
    'XTickLabel', {'Control Early', 'Control Late', 'DMS lesion Early', 'DMS lesion Late'}, ...
    'YLim', [0 yMaxLat * 1.20]);

ylabel(axLatBar, 'Latency (ms)', 'FontSize', barLabelFontSize);
title(axLatBar, sprintf('Early vs Late latency (%d trials each)', EARLY_LATE_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', barTitleFontSize);

set(axLatBar, 'FontSize', barTickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', barAxesTickLW, ...
    'TickLength', barAxesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

xtickangle(axLatBar, 45);

addSigIfNeeded_(axLatBar, 1, 2, pCtrlLat_signrank, [mCtrlLatEarly mCtrlLatLate], barAxesTickLW, barTickFontSize);
addSigIfNeeded_(axLatBar, 4, 5, pLesLat_signrank,  [mLesLatEarly  mLesLatLate],  barAxesTickLW, barTickFontSize);

latBarResults = struct();
latBarResults.matFiles              = matFiles;
latBarResults.groupLabels           = groupLabels;
latBarResults.nTrials               = EARLY_LATE_NTRIALS;
latBarResults.commonEndLatencyTrials = commonEndLatencyTrials;
latBarResults.ctrlLatEarly          = ctrlLatEarly;
latBarResults.ctrlLatLate           = ctrlLatLate;
latBarResults.ctrlLatDelta          = ctrlLatDelta;
latBarResults.lesLatEarly           = lesLatEarly;
latBarResults.lesLatLate            = lesLatLate;
latBarResults.lesLatDelta           = lesLatDelta;
latBarResults.meanVals              = latBarVals;
latBarResults.semVals               = latSemVals;
latBarResults.pCtrlLat_signrank     = pCtrlLat_signrank;
latBarResults.pLesLat_signrank      = pLesLat_signrank;
latBarResults.pLatDelta_ranksum     = pLatDelta_ranksum;
latBarResults.statsCtrlLat_signrank = statsCtrlLat_signrank;
latBarResults.statsLesLat_signrank  = statsLesLat_signrank;
latBarResults.statsLatDelta_ranksum = statsLatDelta_ranksum;

save(latBarOutMat, 'latBarResults', '-v7.3');

set(figLatBar, 'PaperPositionMode', 'auto');
exportgraphics(figLatBar, latBarOutPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figLatBar, latBarOutFig);
try
    saveas(figLatBar, latBarOutSvg);
catch
    print(figLatBar, latBarOutSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', latBarOutPng, latBarOutFig, latBarOutSvg, latBarOutMat);

figLatDelta = figure('Color', 'w', 'Position', [120, 120, 520, 700]);
axLatDelta = axes(figLatDelta); hold(axLatDelta, 'on');

mCtrlLatDelta = mean(ctrlLatDelta, 'omitnan');
mLesLatDelta  = mean(lesLatDelta,  'omitnan');

semCtrlLatDelta = std(ctrlLatDelta, 'omitnan') / sqrt(max(1, numel(ctrlLatDelta)));
semLesLatDelta  = std(lesLatDelta,  'omitnan') / sqrt(max(1, numel(lesLatDelta)));

bdLat = bar(axLatDelta, xDelta, [mCtrlLatDelta mLesLatDelta], deltaBarW, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', barAxesTickLW);
bdLat.CData(1,:) = [0.45 0.45 0.45];
bdLat.CData(2,:) = lesColor;
bdLat.FaceAlpha = barFaceAlpha;

errorbar(axLatDelta, xDelta, [mCtrlLatDelta mLesLatDelta], [semCtrlLatDelta semLesLatDelta], 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', barAxesTickLW, ...
    'CapSize', 18);

rng(12);
jitCtrlLatDelta = 1 + (rand(size(ctrlLatDelta)) - 0.5) * 2 * jit;
jitLesLatDelta  = 2 + (rand(size(lesLatDelta))  - 0.5) * 2 * jit;

plot(axLatDelta, jitCtrlLatDelta, ctrlLatDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(axLatDelta, jitLesLatDelta, lesLatDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

yline(axLatDelta, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

yAllLatDelta = [ctrlLatDelta(:); lesLatDelta(:); mCtrlLatDelta; mLesLatDelta];
yAllLatDelta = yAllLatDelta(isfinite(yAllLatDelta));
if isempty(yAllLatDelta)
    ylimLatDelta = [-1 1];
else
    padLatDelta = 0.20 * max(abs(yAllLatDelta));
    if padLatDelta == 0
        padLatDelta = 1;
    end
    ylimLatDelta = [min(yAllLatDelta)-padLatDelta, max(yAllLatDelta)+padLatDelta];
end

set(axLatDelta, 'XLim', [0.5 2.5], ...
    'XTick', [1 2], ...
    'XTickLabel', {'Control', 'DMS lesion'}, ...
    'YLim', ylimLatDelta);

ylabel(axLatDelta, '\Delta Latency (Late - Early, ms)', 'FontSize', barLabelFontSize);
title(axLatDelta, sprintf('Change from Early to Late latency (%d trials each)', EARLY_LATE_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', barTitleFontSize);

set(axLatDelta, 'FontSize', barTickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', barAxesTickLW, ...
    'TickLength', barAxesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

addSigIfNeeded_(axLatDelta, 1, 2, pLatDelta_ranksum, [mCtrlLatDelta mLesLatDelta], barAxesTickLW, barTickFontSize);

latDeltaResults = struct();
latDeltaResults.matFiles               = matFiles;
latDeltaResults.groupLabels            = groupLabels;
latDeltaResults.nTrials                = EARLY_LATE_NTRIALS;
latDeltaResults.commonEndLatencyTrials = commonEndLatencyTrials;
latDeltaResults.ctrlLatDelta           = ctrlLatDelta;
latDeltaResults.lesLatDelta            = lesLatDelta;
latDeltaResults.meanCtrlLatDelta       = mCtrlLatDelta;
latDeltaResults.meanLesLatDelta        = mLesLatDelta;
latDeltaResults.semCtrlLatDelta        = semCtrlLatDelta;
latDeltaResults.semLesLatDelta         = semLesLatDelta;
latDeltaResults.pLatDelta_ranksum      = pLatDelta_ranksum;
latDeltaResults.statsLatDelta_ranksum  = statsLatDelta_ranksum;

save(latDeltaOutMat, 'latDeltaResults', '-v7.3');

set(figLatDelta, 'PaperPositionMode', 'auto');
exportgraphics(figLatDelta, latDeltaOutPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figLatDelta, latDeltaOutFig);
try
    saveas(figLatDelta, latDeltaOutSvg);
catch
    print(figLatDelta, latDeltaOutSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', latDeltaOutPng, latDeltaOutFig, latDeltaOutSvg, latDeltaOutMat);

end

%% ========================= HELPERS =========================

function fitRes = initEmptyFitResult_()
    fitRes = struct( ...
        'name', '', ...
        'group', '', ...
        'success', false, ...
        'nSessionsUsed', 0, ...
        'B', NaN, ...
        'A', NaN, ...
        'lambda', NaN, ...
        'k', NaN, ...
        'T10', NaN, ...
        'T50', NaN, ...
        'T90', NaN, ...
        'A_CI', [NaN NaN], ...
        'asymSess', NaN, ...
        'groupPct', [], ...
        'groupW', []);
end

function [pHat, A_CI, asymSess] = fitWeibullWeighted_(x, y, w, weibullFun, LB, UB, opts, nBoot)

    B0 = max(0,  min(100, median(y(1:min(3,end))) ));
    A0 = max(0,  min(100, median(y(max(1,end-2):end)) ));
    if A0 < B0 + 1
        A0 = min(100, B0 + max(5, std(y)));
    end
    lambda0 = max(1, round(median(x)));
    k0      = 2;
    p0 = [B0, A0, lambda0, k0];

    obj = @(p) weightedSSE_withBounds_(p, x, y, w, weibullFun, LB, UB);

    pHat = fminsearch(obj, p0, opts);
    pHat = min(UB, max(LB, pHat));

    Aboot = nan(nBoot,1);
    pSamp = w(:) / sum(w);

    for b = 1:nBoot
        idx = randsample(numel(x), numel(x), true, pSamp);

        xb = x(idx);
        yb = y(idx);
        wb = w(idx);

        objb = @(p) weightedSSE_withBounds_(p, xb, yb, wb, weibullFun, LB, UB);
        ph   = fminsearch(objb, pHat, opts);
        ph   = min(UB, max(LB, ph));

        Aboot(b) = ph(2);
    end

    A_CI = prctile(Aboot, [2.5 97.5]);
    A_lo = A_CI(1);

    asymSess = NaN;
    for i = 1:numel(x)
        if y(i) >= A_lo
            asymSess = x(i);
            break
        end
    end
end

function sse = weightedSSE_withBounds_(p, x, y, w, f, LB, UB)
    pen = 0;
    for i = 1:numel(p)
        if p(i) < LB(i)
            pen = pen + (LB(i)-p(i))^2 * 1e6;
        elseif p(i) > UB(i)
            pen = pen + (p(i)-UB(i))^2 * 1e6;
        end
    end

    yhat = f(p, x);
    r = y - yhat;
    sse = sum(w .* (r.^2)) + pen;

    if p(2) < p(1)
        sse = sse + (p(1)-p(2))^2 * 1e6;
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
        yl = ylim(ax);
        yRange = yl(2) - yl(1);
        y = yMaxData + 0.08*yRange;
        h = 0.03*yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01*yRange, stars, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0]);
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

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct','ratBEHstruct_unit','rat_BEHstruct_unit'};

    for k = 1:numel(cands)
        if isfield(S, cands{k}) && isstruct(S.(cands{k}))
            beh = S.(cands{k});
            return;
        end
    end

    f = fieldnames(S);
    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v) && numel(v) > 1
            beh = v;
            return;
        end
    end

    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v)
            beh = v;
            return;
        end
    end
end

function latVals = extractLatenciesFromSession_(sess)
    latVals = [];

    cueRaw = [];
    if isfield(sess, 'CuedTimes') && ~isempty(sess.CuedTimes)
        cueRaw = sess.CuedTimes;
    elseif isfield(sess, 'cuedTimes') && ~isempty(sess.cuedTimes)
        cueRaw = sess.cuedTimes;
    end

    pokeRaw = [];
    if isfield(sess, 'pokeTimes') && ~isempty(sess.pokeTimes)
        pokeRaw = sess.pokeTimes;
    end

    if isempty(cueRaw) || isempty(pokeRaw)
        return;
    end

    nItems = min(numel(cueRaw), numel(pokeRaw));
    if nItems == 0
        return;
    end

    lat = nan(nItems,1);

    for i = 1:nItems
        cueT  = getFirstScalarAt_(cueRaw, i);
        pokeT = getFirstScalarAt_(pokeRaw, i);

        if ~isfinite(cueT) || ~isfinite(pokeT)
            continue;
        end

        thisLat = pokeT - cueT;
        if isfinite(thisLat) && thisLat >= 0
            lat(i) = thisLat;
        end
    end

    latVals = lat(isfinite(lat));
end

function v = getFirstScalarAt_(x, idx)
    v = NaN;

    if isempty(x) || idx < 1 || idx > numel(x)
        return;
    end

    if iscell(x)
        item = x{idx};
    else
        item = x(idx);
    end

    v = firstNumericValue_(item);
end

function v = firstNumericValue_(item)
    v = NaN;

    if isempty(item)
        return;
    end

    if isnumeric(item)
        item = double(item(:));
        if ~isempty(item)
            v = item(1);
        end
        return;
    end

    if iscell(item)
        if ~isempty(item)
            v = firstNumericValue_(item{1});
        end
        return;
    end

    if isstring(item) || ischar(item)
        vv = str2double(string(item));
        if isfinite(vv)
            v = vv;
        end
        return;
    end
end