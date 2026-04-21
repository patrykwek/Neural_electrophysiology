function plot_cue_action_confusion_early_late__DMS_FILES()
%% Plot cue->action confusion matrices for Early and Late epochs
% and summarize diagonal accuracy statistically.
%
% This script:
%   1) Loads behavior files
%   2) Extracts cue labels from cuedNames and chosen-action labels from pokeNames
%   3) Keeps only trials with valid L/C/R labels
%   4) Uses:
%        - Early  = first 1000 valid trials for each animal
%        - Late   = last 1000 valid trials from the COMMON plotted end,
%                   i.e. from the first commonEndTrials trials across animals
%   5) Builds per-animal row-normalized 3x3 confusion matrices:
%        rows = cue (L,C,R), cols = chosen action / poke (L,C,R)
%   6) Averages those matrices within group to make 4 heatmaps:
%        - Control Early
%        - Control Late
%        - DMS lesion Early
%        - DMS lesion Late
%   7) Computes per-animal diagonal accuracy summary:
%        mean(diag(row-normalized confusion matrix)) * 100
%   8) Runs animal-level statistics:
%        - Control Early vs Late: Wilcoxon signed-rank
%        - DMS lesion Early vs Late: Wilcoxon signed-rank
%        - Control delta vs DMS lesion delta: Wilcoxon rank-sum
%   9) Plots:
%        - 4 confusion matrices
%        - Early/Late paired barplot
%        - Delta plot
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
%   cue_action_confusion_early_late.png
%   cue_action_confusion_early_late.svg
%   cue_action_confusion_early_late.fig
%   cue_action_confusion_early_late_results.mat

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

outBase = fullfile(outDir, 'cue_action_confusion_early_late');
outPng  = [outBase '.png'];
outSvg  = [outBase '.svg'];
outFig  = [outBase '.fig'];
outMat  = [outBase '_results.mat'];

EPOCH_NTRIALS = 1000;
labelsLCROrder = {'L','C','R'};

%% ---------------- STYLE ----------------
figPos = [120, 120, 1500, 900];

titleFontSize  = 26;
labelFontSize  = 24;
tickFontSize   = 22;
smallTitleFS   = 22;

axesTickLW   = 4.0;
axesTickLen  = [0.02 0.02];
heatTextFS   = 18;

scatterMS     = 10;
scatterEdgeLW = 1.0;
jit           = 0.06;

colCtrlEarly = [0.75 0.75 0.75];
colCtrlLate  = [0.35 0.35 0.35];
colLesEarly  = [1.00 0.70 0.70];
colLesLate   = [1.00 0.00 0.00];

ctrlColor = [0 0 0];
lesColor  = [1 0 0];
colLine   = [0.45 0.45 0.45];
greyDot   = [0.60 0.60 0.60];
barFaceAlpha = 0.25;

%% ---------------- LOAD ALL FILES / EXTRACT VALID CUE-ACTION TRIALS ----------------
nFiles = numel(matFiles);
animalData = struct([]);

for iFile = 1:nFiles
    matFile = matFiles{iFile};
    assert(exist(matFile,'file') == 2, 'File not found: %s', matFile);

    S = load(matFile);
    beh = pickBehStruct_(S);
    assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in: %s', matFile);

    [~, baseName, ~] = fileparts(matFile);

    cueAll  = cell(0,1);
    pokeAll = cell(0,1);

    for sIdx = 1:numel(beh)
        [cueSess, pokeSess] = extractCuePokeLabelsFromSession_(beh(sIdx));
        if isempty(cueSess), continue; end

        cueAll  = [cueAll;  cueSess(:)];  %#ok<AGROW>
        pokeAll = [pokeAll; pokeSess(:)]; %#ok<AGROW>
    end

    nValidTrials = numel(cueAll);

    animalData(iFile).matFile     = matFile;
    animalData(iFile).baseName    = baseName;
    animalData(iFile).group       = groupLabels{iFile};
    animalData(iFile).cueAll      = cueAll;
    animalData(iFile).pokeAll     = pokeAll;
    animalData(iFile).nValidTrials = nValidTrials;

    fprintf('%s: valid L/C/R cue-action trials = %d\n', baseName, nValidTrials);
end

allNValid = [animalData.nValidTrials];
commonEndTrials = min(allNValid);

fprintf('\nCommon plotted end across animals = %d valid trials\n', commonEndTrials);

assert(commonEndTrials >= EPOCH_NTRIALS, ...
    'Not enough common trials. Need at least %d, but common end is %d.', ...
    EPOCH_NTRIALS, commonEndTrials);

%% ---------------- BUILD EARLY / LATE MATRICES PER ANIMAL ----------------
for iFile = 1:nFiles
    cueAll  = animalData(iFile).cueAll;
    pokeAll = animalData(iFile).pokeAll;

    cueUse  = cueAll(1:commonEndTrials);
    pokeUse = pokeAll(1:commonEndTrials);

    cueEarly  = cueUse(1:EPOCH_NTRIALS);
    pokeEarly = pokeUse(1:EPOCH_NTRIALS);

    cueLate   = cueUse(end-EPOCH_NTRIALS+1:end);
    pokeLate  = pokeUse(end-EPOCH_NTRIALS+1:end);

    confEarly = makeConfusionMatrix_(cueEarly, pokeEarly);
    confLate  = makeConfusionMatrix_(cueLate,  pokeLate);

    diagEarly = 100 * mean(diag(confEarly), 'omitnan');
    diagLate  = 100 * mean(diag(confLate),  'omitnan');
    deltaDiag = diagLate - diagEarly;

    animalData(iFile).cueUse      = cueUse;
    animalData(iFile).pokeUse     = pokeUse;
    animalData(iFile).confEarly   = confEarly;
    animalData(iFile).confLate    = confLate;
    animalData(iFile).diagEarly   = diagEarly;
    animalData(iFile).diagLate    = diagLate;
    animalData(iFile).deltaDiag   = deltaDiag;
end

%% ---------------- GROUP-AVERAGED CONFUSION MATRICES ----------------
isCtrl = strcmp({animalData.group}, 'Control');
isLes  = strcmp({animalData.group}, 'DMS lesion');

ctrlEarlyMats = cat(3, animalData(isCtrl).confEarly);
ctrlLateMats  = cat(3, animalData(isCtrl).confLate);
lesEarlyMats  = cat(3, animalData(isLes).confEarly);
lesLateMats   = cat(3, animalData(isLes).confLate);

ctrlEarlyMean = mean(ctrlEarlyMats, 3, 'omitnan');
ctrlLateMean  = mean(ctrlLateMats,  3, 'omitnan');
lesEarlyMean  = mean(lesEarlyMats,  3, 'omitnan');
lesLateMean   = mean(lesLateMats,   3, 'omitnan');

ctrlDiagEarly = [animalData(isCtrl).diagEarly]';
ctrlDiagLate  = [animalData(isCtrl).diagLate]';
ctrlDiagDelta = [animalData(isCtrl).deltaDiag]';

lesDiagEarly  = [animalData(isLes).diagEarly]';
lesDiagLate   = [animalData(isLes).diagLate]';
lesDiagDelta  = [animalData(isLes).deltaDiag]';

%% ---------------- STATISTICAL TESTS ----------------
[pCtrl_signrank, ~, statsCtrl_signrank] = signrank(ctrlDiagEarly, ctrlDiagLate);
[pLes_signrank,  ~, statsLes_signrank]  = signrank(lesDiagEarly,  lesDiagLate);
[pDelta_ranksum, ~, statsDelta_ranksum] = ranksum(ctrlDiagDelta, lesDiagDelta);

fprintf('\n=== DIAGONAL ACCURACY SUMMARY (%s trials per epoch) ===\n', num2str(EPOCH_NTRIALS));
fprintf('Control Early mean = %.6f | Late mean = %.6f | Delta = %.6f\n', ...
    mean(ctrlDiagEarly, 'omitnan'), mean(ctrlDiagLate, 'omitnan'), mean(ctrlDiagDelta, 'omitnan'));
fprintf('DMS lesion Early mean = %.6f | Late mean = %.6f | Delta = %.6f\n', ...
    mean(lesDiagEarly, 'omitnan'), mean(lesDiagLate, 'omitnan'), mean(lesDiagDelta, 'omitnan'));

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

fprintf('\nControl diagonal deltas:\n');
disp(ctrlDiagDelta(:)');
fprintf('DMS lesion diagonal deltas:\n');
disp(lesDiagDelta(:)');

%% ---------------- FIGURE ----------------
fig = figure('Color', 'w', 'Position', figPos);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

cmap = plasma(256);

% ---- Control Early ----
ax1 = nexttile(tl, 1); hold(ax1, 'on');
imagesc(ax1, ctrlEarlyMean);
axis(ax1, 'image');
set(ax1, 'YDir', 'normal');
colormap(ax1, cmap);
caxis(ax1, [0 1]);
title(ax1, sprintf('Control Early\n(first %d trials)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', smallTitleFS);
xlabel(ax1, 'Chosen action', 'FontSize', labelFontSize);
ylabel(ax1, 'Cue', 'FontSize', labelFontSize);
set(ax1, 'XTick', 1:3, 'XTickLabel', labelsLCROrder, ...
         'YTick', 1:3, 'YTickLabel', labelsLCROrder, ...
         'FontSize', tickFontSize, 'TickDir', 'out', ...
         'LineWidth', axesTickLW, 'TickLength', axesTickLen, ...
         'Box', 'off', 'Layer', 'top');
addHeatText_(ax1, ctrlEarlyMean, heatTextFS);

% ---- Control Late ----
ax2 = nexttile(tl, 2); hold(ax2, 'on');
imagesc(ax2, ctrlLateMean);
axis(ax2, 'image');
set(ax2, 'YDir', 'normal');
colormap(ax2, cmap);
caxis(ax2, [0 1]);
title(ax2, sprintf('Control Late\n(last %d common-end trials)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', smallTitleFS);
xlabel(ax2, 'Chosen action', 'FontSize', labelFontSize);
ylabel(ax2, 'Cue', 'FontSize', labelFontSize);
set(ax2, 'XTick', 1:3, 'XTickLabel', labelsLCROrder, ...
         'YTick', 1:3, 'YTickLabel', labelsLCROrder, ...
         'FontSize', tickFontSize, 'TickDir', 'out', ...
         'LineWidth', axesTickLW, 'TickLength', axesTickLen, ...
         'Box', 'off', 'Layer', 'top');
addHeatText_(ax2, ctrlLateMean, heatTextFS);

% ---- Lesion Early ----
ax3 = nexttile(tl, 4); hold(ax3, 'on');
imagesc(ax3, lesEarlyMean);
axis(ax3, 'image');
set(ax3, 'YDir', 'normal');
colormap(ax3, cmap);
caxis(ax3, [0 1]);
title(ax3, sprintf('DMS lesion Early\n(first %d trials)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', smallTitleFS);
xlabel(ax3, 'Chosen action', 'FontSize', labelFontSize);
ylabel(ax3, 'Cue', 'FontSize', labelFontSize);
set(ax3, 'XTick', 1:3, 'XTickLabel', labelsLCROrder, ...
         'YTick', 1:3, 'YTickLabel', labelsLCROrder, ...
         'FontSize', tickFontSize, 'TickDir', 'out', ...
         'LineWidth', axesTickLW, 'TickLength', axesTickLen, ...
         'Box', 'off', 'Layer', 'top');
addHeatText_(ax3, lesEarlyMean, heatTextFS);

% ---- Lesion Late ----
ax4 = nexttile(tl, 5); hold(ax4, 'on');
imagesc(ax4, lesLateMean);
axis(ax4, 'image');
set(ax4, 'YDir', 'normal');
colormap(ax4, cmap);
caxis(ax4, [0 1]);
title(ax4, sprintf('DMS lesion Late\n(last %d common-end trials)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', smallTitleFS);
xlabel(ax4, 'Chosen action', 'FontSize', labelFontSize);
ylabel(ax4, 'Cue', 'FontSize', labelFontSize);
set(ax4, 'XTick', 1:3, 'XTickLabel', labelsLCROrder, ...
         'YTick', 1:3, 'YTickLabel', labelsLCROrder, ...
         'FontSize', tickFontSize, 'TickDir', 'out', ...
         'LineWidth', axesTickLW, 'TickLength', axesTickLen, ...
         'Box', 'off', 'Layer', 'top');
addHeatText_(ax4, lesLateMean, heatTextFS);

cb = colorbar(ax4, 'eastoutside');
cb.TickDirection = 'out';
cb.LineWidth = axesTickLW;
cb.Label.String = 'P(chosen action | cue)';
cb.FontSize = 18;
cb.Label.FontSize = 18;
cb.Ticks = [0 0.5 1];

% ---- Diagonal accuracy barplot ----
ax5 = nexttile(tl, 3); hold(ax5, 'on');

mCtrlEarly = mean(ctrlDiagEarly, 'omitnan');
mCtrlLate  = mean(ctrlDiagLate,  'omitnan');
mLesEarly  = mean(lesDiagEarly,  'omitnan');
mLesLate   = mean(lesDiagLate,   'omitnan');

semCtrlEarly = std(ctrlDiagEarly, 'omitnan') / sqrt(max(1, numel(ctrlDiagEarly)));
semCtrlLate  = std(ctrlDiagLate,  'omitnan') / sqrt(max(1, numel(ctrlDiagLate)));
semLesEarly  = std(lesDiagEarly,  'omitnan') / sqrt(max(1, numel(lesDiagEarly)));
semLesLate   = std(lesDiagLate,   'omitnan') / sqrt(max(1, numel(lesDiagLate)));

xPos = [1 2 4 5];
barVals = [mCtrlEarly, mCtrlLate, mLesEarly, mLesLate];
semVals = [semCtrlEarly, semCtrlLate, semLesEarly, semLesLate];

b = bar(ax5, xPos, barVals, 0.60, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', axesTickLW);
b.CData(1,:) = colCtrlEarly;
b.CData(2,:) = colCtrlLate;
b.CData(3,:) = colLesEarly;
b.CData(4,:) = colLesLate;
b.FaceAlpha = barFaceAlpha;

errorbar(ax5, xPos, barVals, semVals, 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', axesTickLW, ...
    'CapSize', 18);

rng(0);
jitCtrlEarly = 1 + (rand(size(ctrlDiagEarly)) - 0.5) * 2 * jit;
jitCtrlLate  = 2 + (rand(size(ctrlDiagLate))  - 0.5) * 2 * jit;

for i = 1:numel(ctrlDiagEarly)
    plot(ax5, [jitCtrlEarly(i) jitCtrlLate(i)], [ctrlDiagEarly(i) ctrlDiagLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(ax5, jitCtrlEarly, ctrlDiagEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(ax5, jitCtrlLate, ctrlDiagLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

rng(1);
jitLesEarly = 4 + (rand(size(lesDiagEarly)) - 0.5) * 2 * jit;
jitLesLate  = 5 + (rand(size(lesDiagLate))  - 0.5) * 2 * jit;

for i = 1:numel(lesDiagEarly)
    plot(ax5, [jitLesEarly(i) jitLesLate(i)], [lesDiagEarly(i) lesDiagLate(i)], '-', ...
        'Color', colLine, 'LineWidth', 1.5);
end

plot(ax5, jitLesEarly, lesDiagEarly, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(ax5, jitLesLate, lesDiagLate, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

set(ax5, 'XLim', [0.3 5.7], ...
    'XTick', xPos, ...
    'XTickLabel', {'Ctrl E', 'Ctrl L', 'Les E', 'Les L'}, ...
    'YLim', [0 100]);
ylabel(ax5, 'Diagonal accuracy (%)', 'FontSize', labelFontSize);
title(ax5, 'Diagonal accuracy', 'FontWeight', 'bold', 'FontSize', smallTitleFS);
set(ax5, 'FontSize', tickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');
xtickangle(ax5, 30);

addSigIfNeeded_(ax5, 1, 2, pCtrl_signrank, [mCtrlEarly mCtrlLate], axesTickLW, tickFontSize);
addSigIfNeeded_(ax5, 4, 5, pLes_signrank,  [mLesEarly  mLesLate],  axesTickLW, tickFontSize);

% ---- Delta plot ----
ax6 = nexttile(tl, 6); hold(ax6, 'on');

mCtrlDelta = mean(ctrlDiagDelta, 'omitnan');
mLesDelta  = mean(lesDiagDelta,  'omitnan');

semCtrlDelta = std(ctrlDiagDelta, 'omitnan') / sqrt(max(1, numel(ctrlDiagDelta)));
semLesDelta  = std(lesDiagDelta,  'omitnan') / sqrt(max(1, numel(lesDiagDelta)));

bd = bar(ax6, [1 2], [mCtrlDelta mLesDelta], 0.60, ...
    'FaceColor', 'flat', ...
    'EdgeColor', [0 0 0], ...
    'LineWidth', axesTickLW);
bd.CData(1,:) = [0.45 0.45 0.45];
bd.CData(2,:) = lesColor;
bd.FaceAlpha = barFaceAlpha;

errorbar(ax6, [1 2], [mCtrlDelta mLesDelta], [semCtrlDelta semLesDelta], 'k', ...
    'LineStyle', 'none', ...
    'LineWidth', axesTickLW, ...
    'CapSize', 18);

rng(2);
jitCtrlDelta = 1 + (rand(size(ctrlDiagDelta)) - 0.5) * 2 * jit;
jitLesDelta  = 2 + (rand(size(lesDiagDelta))  - 0.5) * 2 * jit;

plot(ax6, jitCtrlDelta, ctrlDiagDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

plot(ax6, jitLesDelta, lesDiagDelta, 'o', ...
    'MarkerSize', scatterMS, ...
    'MarkerFaceColor', greyDot, ...
    'MarkerEdgeColor', greyDot, ...
    'LineWidth', scatterEdgeLW);

yline(ax6, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

set(ax6, 'XLim', [0.5 2.5], ...
    'XTick', [1 2], ...
    'XTickLabel', {'Control', 'DMS lesion'});
ylabel(ax6, '\Delta diagonal accuracy (Late - Early)', 'FontSize', labelFontSize);
title(ax6, 'Change in diagonal accuracy', 'FontWeight', 'bold', 'FontSize', smallTitleFS);
set(ax6, 'FontSize', tickFontSize, ...
    'TickDir', 'out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box', 'off', ...
    'Layer', 'top');

addSigIfNeeded_(ax6, 1, 2, pDelta_ranksum, [mCtrlDelta mLesDelta], axesTickLW, tickFontSize);

sgtitle(tl, sprintf('Cue-action mapping summary (%d-trial Early/Late epochs)', EPOCH_NTRIALS), ...
    'FontWeight', 'bold', 'FontSize', titleFontSize);

%% ---------------- SAVE ----------------
results = struct();
results.matFiles              = matFiles;
results.groupLabels           = groupLabels;
results.epochNTrials          = EPOCH_NTRIALS;
results.commonEndTrials       = commonEndTrials;
results.animalData            = animalData;
results.ctrlEarlyMean         = ctrlEarlyMean;
results.ctrlLateMean          = ctrlLateMean;
results.lesEarlyMean          = lesEarlyMean;
results.lesLateMean           = lesLateMean;
results.ctrlDiagEarly         = ctrlDiagEarly;
results.ctrlDiagLate          = ctrlDiagLate;
results.ctrlDiagDelta         = ctrlDiagDelta;
results.lesDiagEarly          = lesDiagEarly;
results.lesDiagLate           = lesDiagLate;
results.lesDiagDelta          = lesDiagDelta;
results.pCtrl_signrank        = pCtrl_signrank;
results.pLes_signrank         = pLes_signrank;
results.pDelta_ranksum        = pDelta_ranksum;
results.statsCtrl_signrank    = statsCtrl_signrank;
results.statsLes_signrank     = statsLes_signrank;
results.statsDelta_ranksum    = statsDelta_ranksum;

save(outMat, 'results', '-v7.3');

set(fig, 'PaperPositionMode', 'auto');
exportgraphics(fig, outPng, 'Resolution', 300, 'BackgroundColor', 'white');
savefig(fig, outFig);
try
    saveas(fig, outSvg);
catch
    print(fig, outSvg, '-dsvg');
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg, outMat);

end

%% ========================= HELPERS =========================

function [cueOut, pokeOut] = extractCuePokeLabelsFromSession_(sess)
    cueOut = cell(0,1);
    pokeOut = cell(0,1);

    if ~isfield(sess, 'cuedNames') || ~isfield(sess, 'pokeNames')
        return;
    end

    cueRaw  = sess.cuedNames;
    pokeRaw = sess.pokeNames;

    cueCell  = toCellArray_(cueRaw);
    pokeCell = toCellArray_(pokeRaw);

    n = min(numel(cueCell), numel(pokeCell));
    if n == 0
        return;
    end

    cueCell  = cueCell(1:n);
    pokeCell = pokeCell(1:n);

    keep = false(n,1);
    cueTmp  = cell(n,1);
    pokeTmp = cell(n,1);

    for i = 1:n
        c = firstLetterLCR_(cueCell{i});
        p = firstLetterLCR_(pokeCell{i});

        if ~isempty(c) && ~isempty(p)
            keep(i) = true;
            cueTmp{i}  = c;
            pokeTmp{i} = p;
        end
    end

    cueOut  = cueTmp(keep);
    pokeOut = pokeTmp(keep);
end

function C = makeConfusionMatrix_(cueLabels, pokeLabels)
    % Rows = cue L/C/R, Cols = chosen action L/C/R
    C = nan(3,3);
    cueOrder = {'L','C','R'};
    pokeOrder = {'L','C','R'};

    for r = 1:3
        rowMask = strcmp(cueLabels, cueOrder{r});
        nRow = sum(rowMask);
        if nRow == 0
            continue;
        end
        for c = 1:3
            C(r,c) = sum(rowMask & strcmp(pokeLabels, pokeOrder{c})) / nRow;
        end
    end
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

    if isnumeric(v)
        x = num2cell(v(:));
        return;
    end

    try
        x = num2cell(v(:));
    catch
        x = {};
    end
end

function lab = firstLetterLCR_(x)
    lab = '';

    if isempty(x)
        return;
    end

    if isstring(x)
        x = char(x);
    elseif isnumeric(x)
        if isscalar(x)
            x = char(string(x));
        else
            return;
        end
    elseif ~ischar(x)
        try
            x = char(string(x));
        catch
            return;
        end
    end

    x = upper(strtrim(x));
    if isempty(x)
        return;
    end

    firstChar = x(1);
    if ismember(firstChar, ['L','C','R'])
        lab = firstChar;
    end
end

function addHeatText_(ax, M, fs)
    for r = 1:size(M,1)
        for c = 1:size(M,2)
            if ~isfinite(M(r,c))
                txt = 'NaN';
                col = [0 0 0];
            else
                txt = sprintf('%.2f', M(r,c));
                if M(r,c) >= 0.60
                    col = [1 1 1];
                else
                    col = [0 0 0];
                end
            end
            text(ax, c, r, txt, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', fs, ...
                'FontWeight', 'bold', ...
                'Color', col);
        end
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