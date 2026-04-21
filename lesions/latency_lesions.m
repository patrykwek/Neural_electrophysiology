function plot_latency_early_late__DMS_FILES()
%% Plot cue->press latency for Early and Late epochs across Control and DMS lesion groups
% and compare them with the same logic/statistics as the diagonal-accuracy summary.
%
% This script:
%   1) Loads behavior files
%   2) Extracts cue times from cuedTimes / CuedTimes and press times from pokeTimes
%   3) Computes latency = firstPressTime - cueTime for each valid trial
%   4) Keeps only finite latencies >= 0
%   5) Uses:
%        - Early = first 1000 valid trials for each animal
%        - Late  = last 1000 valid trials from the COMMON plotted end,
%                  i.e. from the first commonEndTrials trials across animals
%   6) Computes per-animal:
%        - mean Early latency
%        - mean Late latency
%        - delta latency = Late - Early
%   7) Runs animal-level statistics:
%        - Control Early vs Late: Wilcoxon signed-rank
%        - DMS lesion Early vs Late: Wilcoxon signed-rank
%        - Control delta vs DMS lesion delta: Wilcoxon rank-sum
%   8) Plots:
%        - paired barplot of Early vs Late latency
%        - delta plot
%      and adds significance stars where p < 0.05
%
% INPUT FILES:
%   Control:
%     /Volumes/WD_BLACK/AB/AB/B1Bay_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat
%     /Volumes/WD_BLACK/AB/AB/F4Fig2_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/J3Jelly_ratBEHstruct_CUEDNAMES_MIN100.mat
%
%   DMS lesion:
%     /Volumes/WD_BLACK/AB/AB/J5Joy_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/J7Jasmine_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/L3Lychee_ratBEHstruct_CUEDNAMES_MIN100.mat
%     /Volumes/WD_BLACK/AB/AB/T8Truffle_ratBEHstruct_CUEDNAMES_MIN100.mat
%
% OUTPUTS SAVED IN:
%   /Volumes/WD_BLACK/AB/AB
%
% OUTPUT FILES:
%   latency_early_late_barplot.png
%   latency_early_late_barplot.svg
%   latency_early_late_barplot.fig
%   latency_early_late_barplot_results.mat
%   latency_delta_plot.png
%   latency_delta_plot.svg
%   latency_delta_plot.fig
%   latency_delta_plot_results.mat

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

barOutBase = fullfile(outDir, 'latency_early_late_barplot');
barOutPng  = [barOutBase '.png'];
barOutSvg  = [barOutBase '.svg'];
barOutFig  = [barOutBase '.fig'];
barOutMat  = [barOutBase '_results.mat'];

deltaOutBase = fullfile(outDir, 'latency_delta_plot');
deltaOutPng  = [deltaOutBase '.png'];
deltaOutSvg  = [deltaOutBase '.svg'];
deltaOutFig  = [deltaOutBase '.fig'];
deltaOutMat  = [deltaOutBase '_results.mat'];

EPOCH_NTRIALS = 1000;

%% ---------------- STYLE ----------------
figPosBar = [120, 120, 640, 760];
figPosDelta = [120, 120, 520, 700];

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

lesColor      = [1 0 0];
colLine       = [0.45 0.45 0.45];
greyDot       = [0.60 0.60 0.60];
barFaceAlpha  = 0.25;

%% ---------------- LOAD FILES / EXTRACT LATENCIES ----------------
nFiles = numel(matFiles);
animalData = struct([]);

for iFile = 1:nFiles
    matFile = matFiles{iFile};
    assert(exist(matFile, 'file') == 2, 'File not found: %s', matFile);

    S = load(matFile);
    beh = pickBehStruct_(S);
    assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in: %s', matFile);

    [~, baseName, ~] = fileparts(matFile);

    allLatencies = [];

    for sIdx = 1:numel(beh)
        latSess = extractLatenciesFromSession_(beh(sIdx));
        if isempty(latSess)
            continue;
        end
        allLatencies = [allLatencies; latSess(:)]; %#ok<AGROW>
    end

    animalData(iFile).matFile      = matFile;
    animalData(iFile).baseName     = baseName;
    animalData(iFile).group        = groupLabels{iFile};
    animalData(iFile).allLatencies = allLatencies;
    animalData(iFile).nValidTrials = numel(allLatencies);

    fprintf('%s: valid latency trials = %d\n', baseName, numel(allLatencies));
end

allNValid = [animalData.nValidTrials];
commonEndTrials = min(allNValid);

fprintf('\nCommon plotted end across animals = %d valid trials\n', commonEndTrials);

assert(commonEndTrials >= EPOCH_NTRIALS, ...
    'Not enough common trials. Need at least %d, but common end is %d.', ...
    EPOCH_NTRIALS, commonEndTrials);

%% ---------------- EARLY / LATE PER ANIMAL ----------------
for iFile = 1:nFiles
    latUse = animalData(iFile).allLatencies(1:commonEndTrials);

    earlyLat = latUse(1:EPOCH_NTRIALS);
    lateLat  = latUse(end-EPOCH_NTRIALS+1:end);

    animalData(iFile).latUse      = latUse;
    animalData(iFile).earlyLat    = earlyLat;
    animalData(iFile).lateLat     = lateLat;
    animalData(iFile).meanEarly   = mean(earlyLat, 'omitnan');
    animalData(iFile).meanLate    = mean(lateLat,  'omitnan');
    animalData(iFile).deltaLat    = animalData(iFile).meanLate - animalData(iFile).meanEarly;
end

isCtrl = strcmp({animalData.group}, 'Control');
isLes  = strcmp({animalData.group}, 'DMS lesion');

ctrlEarly = [animalData(isCtrl).meanEarly]';
ctrlLate  = [animalData(isCtrl).meanLate]';
ctrlDelta = [animalData(isCtrl).deltaLat]';

lesEarly  = [animalData(isLes).meanEarly]';
lesLate   = [animalData(isLes).meanLate]';
lesDelta  = [animalData(isLes).deltaLat]';

%% ---------------- STATISTICAL TESTS ----------------
[pCtrl_signrank, ~, statsCtrl_signrank] = signrank(ctrlEarly, ctrlLate);
[pLes_signrank,  ~, statsLes_signrank]  = signrank(lesEarly, lesLate);
[pDelta_ranksum, ~, statsDelta_ranksum] = ranksum(ctrlDelta, lesDelta);

mCtrlEarly = mean(ctrlEarly, 'omitnan');
mCtrlLate  = mean(ctrlLate,  'omitnan');
mLesEarly  = mean(lesEarly,  'omitnan');
mLesLate   = mean(lesLate,   'omitnan');

semCtrlEarly = std(ctrlEarly, 'omitnan') / sqrt(max(1, numel(ctrlEarly)));
semCtrlLate  = std(ctrlLate,  'omitnan') / sqrt(max(1, numel(ctrlLate)));
semLesEarly  = std(lesEarly,  'omitnan') / sqrt(max(1, numel(lesEarly)));
semLesLate   = std(lesLate,   'omitnan') / sqrt(max(1, numel(lesLate)));

mCtrlDelta = mean(ctrlDelta, 'omitnan');
mLesDelta  = mean(lesDelta,  'omitnan');

semCtrlDelta = std(ctrlDelta, 'omitnan') / sqrt(max(1, numel(ctrlDelta)));
semLesDelta  = std(lesDelta,  'omitnan') / sqrt(max(1, numel(lesDelta)));

fprintf('\n=== LATENCY SUMMARY (%d trials per epoch) ===\n', EPOCH_NTRIALS);
fprintf('Control Early mean = %.6f ms | Late mean = %.6f ms | Delta = %.6f ms\n', ...
    mCtrlEarly, mCtrlLate, mean(ctrlDelta, 'omitnan'));
fprintf('DMS lesion Early mean = %.6f ms | Late mean = %.6f ms | Delta = %.6f ms\n', ...
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

fprintf('\nControl latency deltas (Late - Early, ms):\n');
disp(ctrlDelta(:)');
fprintf('DMS lesion latency deltas (Late - Early, ms):\n');
disp(lesDelta(:)');

%% ---------------- EARLY vs LATE BARPLOT ----------------
figBar = figure('Color', 'w', 'Position', figPosBar);
axBar = axes(figBar); hold(axBar, 'on');

xPos = [1 2 4 5];
barVals = [mCtrlEarly, mCtrlLate, mLesEarly, mLesLate];
semVals = [semCtrlEarly, semCtrlLate, semLesEarly, semLesLate];

b = bar(axBar, xPos, barVals, 0.60, ...
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

ylBar = [0, max([barVals + semVals, ctrlEarly', ctrlLate', lesEarly', lesLate']) * 1.20];
if ~all(isfinite(ylBar)) || ylBar(2) <= ylBar(1)
    ylBar = [0 1000];
end

set(axBar, 'XLim', [0.3 5.7], ...
    'XTick', xPos, ...
    'XTickLabel', {'Control Early', 'Control Late', 'DMS lesion Early', 'DMS lesion Late'}, ...
    'YLim', ylBar);

ylabel(axBar, 'Latency (ms)', 'FontSize', barLabelFontSize);
title(axBar, sprintf('Early vs Late latency (%d trials each)', EPOCH_NTRIALS), ...
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
barResults.matFiles            = matFiles;
barResults.groupLabels         = groupLabels;
barResults.epochNTrials        = EPOCH_NTRIALS;
barResults.commonEndTrials     = commonEndTrials;
barResults.animalData          = animalData;
barResults.ctrlEarly           = ctrlEarly;
barResults.ctrlLate            = ctrlLate;
barResults.ctrlDelta           = ctrlDelta;
barResults.lesEarly            = lesEarly;
barResults.lesLate             = lesLate;
barResults.lesDelta            = lesDelta;
barResults.meanVals            = barVals;
barResults.semVals             = semVals;
barResults.pCtrl_signrank      = pCtrl_signrank;
barResults.pLes_signrank       = pLes_signrank;
barResults.pDelta_ranksum      = pDelta_ranksum;
barResults.statsCtrl_signrank  = statsCtrl_signrank;
barResults.statsLes_signrank   = statsLes_signrank;
barResults.statsDelta_ranksum  = statsDelta_ranksum;

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
figDelta = figure('Color', 'w', 'Position', figPosDelta);
axDelta = axes(figDelta); hold(axDelta, 'on');

bd = bar(axDelta, [1 2], [mCtrlDelta mLesDelta], 0.60, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', barAxesTickLW);
bd.CData(1,:) = [0.45 0.45 0.45];
bd.CData(2,:) = lesColor;
bd.FaceAlpha = barFaceAlpha;

errorbar(axDelta, [1 2], [mCtrlDelta mLesDelta], [semCtrlDelta semLesDelta], 'k', ...
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

yAll = [ctrlDelta(:); lesDelta(:); mCtrlDelta-semCtrlDelta; mCtrlDelta+semCtrlDelta; mLesDelta-semLesDelta; mLesDelta+semLesDelta; 0];
yAll = yAll(isfinite(yAll));
if isempty(yAll)
    ylDelta = [-100 100];
else
    pad = 0.20 * max(abs(yAll));
    if pad == 0, pad = 50; end
    ylDelta = [min(yAll)-pad, max(yAll)+pad];
end

set(axDelta, 'XLim', [0.5 2.5], ...
    'XTick', [1 2], ...
    'XTickLabel', {'Control', 'DMS lesion'}, ...
    'YLim', ylDelta);

ylabel(axDelta, '\Delta Latency (Late - Early, ms)', 'FontSize', barLabelFontSize);
title(axDelta, sprintf('Change from Early to Late latency (%d trials each)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', barTitleFontSize);

set(axDelta, 'FontSize', barTickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', barAxesTickLW, ...
    'TickLength', barAxesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

addSigIfNeeded_(axDelta, 1, 2, pDelta_ranksum, [mCtrlDelta mLesDelta], barAxesTickLW, barTickFontSize);

deltaResults = struct();
deltaResults.matFiles            = matFiles;
deltaResults.groupLabels         = groupLabels;
deltaResults.epochNTrials        = EPOCH_NTRIALS;
deltaResults.commonEndTrials     = commonEndTrials;
deltaResults.ctrlDelta           = ctrlDelta;
deltaResults.lesDelta            = lesDelta;
deltaResults.meanCtrlDelta       = mCtrlDelta;
deltaResults.meanLesDelta        = mLesDelta;
deltaResults.semCtrlDelta        = semCtrlDelta;
deltaResults.semLesDelta         = semLesDelta;
deltaResults.pDelta_ranksum      = pDelta_ranksum;
deltaResults.statsDelta_ranksum  = statsDelta_ranksum;

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

end

%% ========================= HELPERS =========================

function lat = extractLatenciesFromSession_(sess)
    lat = [];

    cueRaw = [];
    if isfield(sess, 'cuedTimes')
        cueRaw = sess.cuedTimes;
    elseif isfield(sess, 'CuedTimes')
        cueRaw = sess.CuedTimes;
    end

    pokeRaw = [];
    if isfield(sess, 'pokeTimes')
        pokeRaw = sess.pokeTimes;
    end

    if isempty(cueRaw) || isempty(pokeRaw)
        return;
    end

    cueCell  = toCellArray_(cueRaw);
    pokeCell = toCellArray_(pokeRaw);

    n = min(numel(cueCell), numel(pokeCell));
    if n == 0
        return;
    end

    latVec = nan(n,1);

    for i = 1:n
        cueT  = firstNumericValue_(cueCell{i});
        pokeT = firstNumericValue_(pokeCell{i});

        if ~isfinite(cueT) || ~isfinite(pokeT)
            continue;
        end

        thisLat = pokeT - cueT;

        if isfinite(thisLat) && thisLat >= 0
            latVec(i) = thisLat;
        end
    end

    lat = latVec(isfinite(latVec));
end

function x = toCellArray_(v)
    if isempty(v)
        x = {};
        return;
    end

    if iscell(v)
        x = v(:);
        return;
    end

    if isnumeric(v)
        if isvector(v)
            x = num2cell(v(:));
        else
            x = num2cell(v, 2);
            x = x(:);
        end
        return;
    end

    if isstring(v)
        x = cellstr(v(:));
        return;
    end

    if ischar(v)
        if size(v,1) > 1
            x = cellstr(v);
        else
            x = {v};
        end
        x = x(:);
        return;
    end

    try
        x = num2cell(v(:));
    catch
        x = {};
    end
end

function v = firstNumericValue_(x)
    v = NaN;

    if isempty(x)
        return;
    end

    if isnumeric(x)
        x = x(:);
        if ~isempty(x)
            v = double(x(1));
        end
        return;
    end

    if iscell(x)
        if isempty(x)
            return;
        end
        v = firstNumericValue_(x{1});
        return;
    end

    if isstring(x)
        x = char(x);
    end

    if ischar(x)
        nums = sscanf(x, '%f');
        if ~isempty(nums)
            v = double(nums(1));
        else
            vv = str2double(strtrim(x));
            if isfinite(vv)
                v = double(vv);
            end
        end
        return;
    end

    try
        vv = double(x);
        vv = vv(:);
        if ~isempty(vv)
            v = vv(1);
        end
    catch
    end
end

function addSigIfNeeded_(ax, x1, x2, p, means, lw, fs)
    if ~isfinite(p) || p >= 0.05
        return;
    end
    stars = pToStars_(p);

    yl = ylim(ax);
    yMaxData = max(means);
    yRange = yl(2) - yl(1);
    y = yMaxData + 0.08 * yRange;
    h = 0.03 * yRange;

    if y + h + 0.05 * yRange > yl(2)
        ylim(ax, [yl(1), y + h + 0.08 * yRange]);
        yl = ylim(ax);
        yRange = yl(2) - yl(1);
        y = yMaxData + 0.08 * yRange;
        h = 0.03 * yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01 * yRange, stars, ...
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