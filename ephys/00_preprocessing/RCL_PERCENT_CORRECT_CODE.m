%% ===== % correct by cueType (L, C, R) + SESSION as unit + paired t-tests + stars + save output =====
clear; clc;

matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/4_CUE_CORRECT';
if ~isfolder(outDir), mkdir(outDir); end

% --- Plot style (match your FIG3E-like aesthetics) ---
figPos  = [120, 120, 1050, 520];

% >>> DOUBLED (as requested) <<<
fsTitle = 40;    % was 20
fsAx    = 32;    % was 16
lwAx    = 5.0;   % was 2.5

S = load(matFile);
beh = pickBehStruct_(S);

cats = ["L","C","R"];
nSess = numel(beh);

% Session x CueType matrix of % correct
P = nan(nSess, numel(cats));

for sIdx = 1:nSess
    if ~isfield(beh(sIdx),'trials') || isempty(beh(sIdx).trials), continue; end
    tr = beh(sIdx).trials;

    if ~isfield(tr,'cueType') || ~isfield(tr,'correct'), continue; end

    % robust extraction
    cue = string({tr.cueType});
    cue = upper(strtrim(cue(:)));

    cor = [tr.correct]';
    cor = double(cor(:));

    keep = ismember(cue, cats) & isfinite(cor);
    cue = cue(keep);
    cor = cor(keep);

    for i = 1:numel(cats)
        idx = (cue == cats(i));
        if any(idx)
            P(sIdx,i) = 100 * mean(cor(idx) == 1);
        end
    end
end

% Use only sessions with ALL 3 cue types present (paired comparisons)
useSess = all(isfinite(P), 2);
P_use   = P(useSess, :);
nUse    = size(P_use, 1);

fprintf('\nSessions total: %d\nSessions used (have L,C,R): %d\n', nSess, nUse);

% Summary for plotting (mean +/- SEM across sessions)
mu  = mean(P_use, 1, 'omitnan');
sem = std(P_use, 0, 1, 'omitnan') ./ sqrt(nUse);

%% ---- Paired Student t-tests (paired by session) + Bonferroni ----
pairs = [1 2; 1 3; 2 3];
pairNames = ["L vs C"; "L vs R"; "C vs R"];

pRaw = nan(3,1);
tStat = nan(3,1);
df = nan(3,1);

fprintf('\n=== Paired t-tests (session as unit) ===\n');
for k = 1:3
    a = P_use(:, pairs(k,1));
    b = P_use(:, pairs(k,2));
    [~, pRaw(k), ~, stats] = ttest(a, b); % paired
    tStat(k) = stats.tstat;
    df(k) = stats.df;
    fprintf('%s: t(%d)=%.4f, p_raw=%.6g\n', pairNames(k), df(k), tStat(k), pRaw(k));
end

pBonf = min(pRaw * 3, 1);
fprintf('\n=== Bonferroni-corrected p-values ===\n');
for k = 1:3
    fprintf('%s: p_bonf=%.6g\n', pairNames(k), pBonf(k));
end

%% ---- Plot (histogram-like bars for L/C/R) ----
fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

xCats = categorical(cats);
b = bar(ax, xCats, mu);

% Colors: L orange, C blue, R purple
b.FaceColor = 'flat';
b.CData(1,:) = [0.90 0.55 0.10]; % L
b.CData(2,:) = [0.10 0.55 0.95]; % C
b.CData(3,:) = [0.55 0.25 0.75]; % R

% Error bars (use the ACTUAL bar centers so the mean matches the bar top)
% >>> DOUBLED (as requested): LineWidth 2 -> 4 <<<
errorbar(ax, b.XEndPoints, mu, sem, 'k.', 'LineWidth', 4);

% Session points (jitter) -- all gray
rng(0);
jitter = 0.12;
for i = 1:numel(cats)
    xj = b.XEndPoints(i) + (rand(nUse,1)-0.5)*2*jitter;
    scatter(ax, xj, P_use(:,i), 28, 'filled', ...
        'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerFaceAlpha', 0.55, ...
        'MarkerEdgeAlpha', 0);
end

ylabel(ax, '% correct trials', 'FontSize', fsAx);
xlabel(ax, 'Cue type',        'FontSize', fsAx);

title(ax, sprintf('%% Correct by cue type', nUse), ...
    'FontSize', fsTitle, 'FontWeight','bold');

ylim(ax, [0 100]);
set(ax, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
box(ax,'off');

%% ---- Significance bars + stars (Bonferroni p) ----
sig = pBonf < 0.05;

if any(sig)
    yMax = max(P_use(:));
    yMin = min(P_use(:));
    yRange = yMax - yMin;
    if ~isfinite(yRange) || yRange == 0, yRange = 1; end

    baseY = max(mu + sem) + 0.06*yRange;  % start above bars
    hBar  = 0.02*yRange;                  % bar height
    stepY = 0.06*yRange;                  % vertical spacing between bars

    % Order comparisons by span (shorter first) to reduce overlap
    order = [1 3 2]; % (L-C), (C-R), (L-R)
    drawCount = 0;

    for kk = 1:3
        k = order(kk);
        if ~sig(k), continue; end
        drawCount = drawCount + 1;

        xA = b.XEndPoints(pairs(k,1));
        xB = b.XEndPoints(pairs(k,2));
        yLevel = baseY + (drawCount-1)*stepY;

        % >>> DOUBLED (as requested): LineWidth 2 -> 4 <<<
        plot(ax, [xA xA xB xB], [yLevel yLevel+hBar yLevel+hBar yLevel], 'k-', 'LineWidth', 4);

        % Stars (no string math; simple if/elseif)
        if pBonf(k) < 0.001
            s = '***';
        elseif pBonf(k) < 0.01
            s = '**';
        elseif pBonf(k) < 0.05
            s = '*';
        else
            s = '';
        end

        % Decreased distance between line and stars (your requested tweak)
        text(ax, mean([xA xB]), yLevel+hBar+0.002*yRange, s, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', fsTitle, 'FontWeight','bold');
    end

    % make sure ylim includes stars
    ylim(ax, [0, min(100, (baseY + drawCount*stepY + 0.06*yRange))]);
end

%% ---- Save outputs ----
outPng = fullfile(outDir, 'PercentCorrect_by_CueType_SESSION_AS_UNIT_TTEST_STARS.png');
outFig = fullfile(outDir, 'PercentCorrect_by_CueType_SESSION_AS_UNIT_TTEST_STARS.fig');
outSvg = fullfile(outDir, 'PercentCorrect_by_CueType_SESSION_AS_UNIT_TTEST_STARS.svg');

saveas(fig, outPng);
savefig(fig, outFig);
saveas(fig, outSvg);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg);

%% ---- helper: find behavior struct in loaded MAT ----
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