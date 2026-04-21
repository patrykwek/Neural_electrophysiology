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
%
% OUTPUTS SAVED IN:
%   /Volumes/WD_BLACK/AB/AB
%
% OUTPUT FILES:
%   accuracy_per_100trials_all_animals.png
%   accuracy_per_100trials_all_animals.svg
%   accuracy_per_100trials_weibull_results.mat
%   accuracy_per_100trials_weibull_results.txt

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