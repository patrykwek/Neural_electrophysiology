%% ========= PLOT_MDS_AND_DRIFT__WEIBULL_AXISSTYLE__SAVE_FIG_SVG.m =========
% PURPOSE:
%   From an existing session×session similarity matrix SIM (Pearson r),
%   make:
%     (1) Classical MDS plot (ALL sessions) with:
%           Early (1–8) GREEN and Late (29–36) YELLOW
%           using the SAME shades as the drift background shading
%           middle sessions very light gray
%           NO session-number labels next to dots
%     (2) Drift plot: r(session t, session 1) with the SAME x-axis tick style
%         as your Weibull script, BUT with labels at:
%           Session 1, Session 6, then every 5th (11,16,21,...)
%         and wherever there is a number label, that tick is longer.
%   Also prints Spearman correlation (drift vs session) + a permutation test
%   for Early vs Late mean drift.
%
%   ADDED (REQUESTED):
%     (3) Early vs Late bar plot (mean +/- SEM) with grey dot scatter
%         using the same aesthetics as your decoding-stage bar plots,
%         and annotates significance using the permutation test p-value.
%
% SAVES:
%   - MDS: PNG + SVG + FIG
%   - Drift: PNG + SVG + FIG
%   - Drift Early vs Late Bar: PNG + SVG + FIG

clear; clc;

%% ---- USER SETTINGS ----
simMatFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/3_4_FIGURE_MATRIX/MATRIX/Sessions01_36_sessionEnsemble_similarity_GLOBALWARP.mat';
outDir     = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/3_4_FIGURE_MATRIX/MATRIX';

earlyGroup = 1:8;
lateGroup  = 29:36;
boundaryAt = 8.5;

useGeometricDistance = true;
connectChronological = true; %#ok<NASGU>  % kept for compatibility; connecting lines removed by request
nPerm = 10000;

if ~exist(outDir,'dir'), mkdir(outDir); end
assert(exist(simMatFile,'file')==2, 'Missing SIM mat file: %s', simMatFile);

rng(0);

%% ---- STYLE (match Weibull script) ----
figPosMDS   = [120, 120, 1100, 650];
figPosDrift = [120, 120, 1200, 650];

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

% ---- EXACT SAME SHADES AS DRIFT BACKGROUND ----
colEarly_MDS = [0 1 0];  % pure green (same as drift shading)
colLate_MDS  = [1 1 0];  % pure yellow (same as drift shading)

colMid_MDS   = [0.88 0.88 0.88];
colLine_MDS  = [0.65 0.65 0.65]; %#ok<NASGU> % kept; lines removed

colDrift     = [0 0.45 0.75];   % keep drift line blue
dashColor    = [0.2 0.2 0.2];
bgAlpha      = 0.10;

%% ---- ADDED: BAR-PLOT STYLE (match decoding-stage bars) ----
% (Use decoding-stage colors + aesthetics)
colEarlyShade = [1 1 0]; % yellow
colLateShade  = [0 1 0]; % green
scatterMS     = 8;
scatterEdgeLW = 1.0;
jit           = 0.10;
greyDot       = [0.6 0.6 0.6];

% --- ADDED outputs (drift early vs late bar) ---
outBaseBar = fullfile(outDir, sprintf('DRIFT_EARLY_LATE_BAR_%02dSessions', 36)); % placeholder; overwritten after nSess known
outBarPng  = ''; outBarSvg = ''; outBarFig = '';

%% ---- LOAD SIM ----
S = load(simMatFile);
assert(isfield(S,'SIM'), 'SIM not found in %s', simMatFile);
SIM = S.SIM;

nSess = size(SIM,1);
assert(size(SIM,2)==nSess, 'SIM must be square');

earlyGroup = earlyGroup(earlyGroup>=1 & earlyGroup<=nSess);
lateGroup  = lateGroup(lateGroup>=1 & lateGroup<=nSess);
midGroup   = setdiff(1:nSess, [earlyGroup(:); lateGroup(:)]);

% --- finalize bar output names now that nSess is known ---
outBaseBar = fullfile(outDir, sprintf('DRIFT_EARLY_LATE_BAR_%02dSessions', nSess));
outBarPng  = [outBaseBar '.png'];
outBarSvg  = [outBaseBar '.svg'];
outBarFig  = [outBaseBar '.fig'];

%% ============================================================
%  (1) MDS (ALL sessions)
%% ============================================================
R = SIM;

Rnan = isnan(R);
R(Rnan) = 0;
R = max(-1, min(1, R));

if useGeometricDistance
    D = sqrt(max(0, 2*(1 - R)));
else
    D = max(0, 1 - R);
end

D = 0.5*(D + D');
D(1:nSess+1:end) = 0;

if any(Rnan(:))
    medD = median(D(~Rnan), 'omitnan');
    D(Rnan) = medD;
    D = 0.5*(D + D');
    D(1:nSess+1:end) = 0;
end

[Y, eigvals] = classical_mds_(D, 2);
x = Y(:,1); y = Y(:,2);

isEarly = ismember(1:nSess, earlyGroup);
isLate  = ismember(1:nSess, lateGroup);
isMid   = ismember(1:nSess, midGroup);

figMDS = figure('Color','w','Position',figPosMDS);
axMDS = axes(figMDS); hold(axMDS,'on');

% (CHANGED) remove connecting line plot entirely per request

% (CHANGED) double point size on MDS plot only
mdsMS = 2 * dataMS;

plot(axMDS, x(isMid), y(isMid), 'o', ...
    'MarkerSize', mdsMS, ...
    'MarkerFaceColor', colMid_MDS, ...
    'MarkerEdgeColor', colMid_MDS, ...
    'LineWidth', markerLW);

plot(axMDS, x(isEarly), y(isEarly), 'o', ...
    'MarkerSize', mdsMS, ...
    'MarkerFaceColor', colEarly_MDS, ...
    'MarkerEdgeColor', colEarly_MDS, ...
    'LineWidth', markerLW);

plot(axMDS, x(isLate), y(isLate), 'o', ...
    'MarkerSize', mdsMS, ...
    'MarkerFaceColor', colLate_MDS, ...
    'MarkerEdgeColor', colLate_MDS, ...
    'LineWidth', markerLW);

title(axMDS, 'Classical MDS embedding of session similarity', ...
    'FontWeight','bold', 'FontSize', titleFontSize);
xlabel(axMDS, 'MDS 1', 'FontSize', labelFontSize);
ylabel(axMDS, 'MDS 2', 'FontSize', labelFontSize);

set(axMDS, 'FontSize', tickFontSize, 'Box','off');
set(axMDS, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

legend(axMDS, {'Middle (9–28)','Early (1–8)','Late (29–36)'}, 'Location','best');

if ~isempty(eigvals) && all(isfinite(eigvals)) && sum(max(eigvals,0))>0
    ve = 100 * max(eigvals,0) ./ sum(max(eigvals,0));
    txt = sprintf('Var explained: %.1f%%, %.1f%%', ve(1), ve(2));
    text(axMDS, 0.02, 0.98, txt, 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 18);
end

outBaseMDS = fullfile(outDir, sprintf('MDS_ALLSESS_EarlyMidLate_%02dSessions', nSess));
exportgraphics(figMDS, [outBaseMDS '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figMDS, [outBaseMDS '.svg']); catch, print(figMDS, [outBaseMDS '.svg'], '-dsvg'); end
savefig(figMDS, [outBaseMDS '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseMDS);

%% ============================================================
%  (2) Drift curve (Weibull x-axis style, UPDATED labels)
%% ============================================================
drift = SIM(:,1);

figDrift = figure('Color','w','Position',figPosDrift);
ax = axes(figDrift); hold(ax,'on');

finiteD = drift(isfinite(drift));
if isempty(finiteD)
    yl = [-1 1];
else
    pad = 0.05;
    yl = [max(-1, min(finiteD)-pad), min(1, max(finiteD)+pad)];
end

set(ax, 'XLim',[1 nSess], 'YLim', yl);

patch(ax, [1 boundaryAt boundaryAt 1], [yl(1) yl(1) yl(2) yl(2)], ...
    [0 1 0], 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(ax, [boundaryAt nSess nSess boundaryAt], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 1 0], 'FaceAlpha', bgAlpha, 'EdgeColor','none');

plot(ax, 1:nSess, drift, '-', 'LineWidth', dataLW, 'Color', colDrift);
plot(ax, 1:nSess, drift, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', dataMS, ...
    'MarkerFaceColor', colDrift, ...
    'MarkerEdgeColor', colDrift, ...
    'LineWidth', markerLW);

xline(ax, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Pearson r to Session 1 template', 'FontSize', labelFontSize);
title(ax, 'Representational drift: similarity to Session 1', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSess], 'YLim', yl);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

%% ---- Manual tick drawing (Weibull style) with UPDATED labels: 1,6,11,16,... ----
ax.XTick = 1:nSess;

labIdx = unique([1 6 11:5:nSess]);

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

outBaseDrift = fullfile(outDir, sprintf('DRIFT_toSession1_%02dSessions', nSess));
exportgraphics(figDrift, [outBaseDrift '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figDrift, [outBaseDrift '.svg']); catch, print(figDrift, [outBaseDrift '.svg'], '-dsvg'); end
savefig(figDrift, [outBaseDrift '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseDrift);

%% ============================================================
%  STATS
%% ============================================================
fprintf('\n=== DRIFT STATS ===\n');

sessIdx = (1:nSess)';
[rho_s, p_s] = corr(sessIdx, drift, 'Type','Spearman', 'Rows','complete');
fprintf('Spearman corr(drift, session): rho=%.4f, p=%.4g\n', rho_s, p_s);

earlyVals = drift(earlyGroup);
lateVals  = drift(lateGroup);

earlyVals = earlyVals(isfinite(earlyVals));
lateVals  = lateVals(isfinite(lateVals));

obsDiff = mean(lateVals) - mean(earlyVals);

allVals = [earlyVals(:); lateVals(:)];
nE = numel(earlyVals);
nAll = numel(allVals);

nullDiff = nan(nPerm,1);
for b = 1:nPerm
    perm = allVals(randperm(nAll));
    e = perm(1:nE);
    l = perm(nE+1:end);
    nullDiff(b) = mean(l) - mean(e);
end

p_perm = (sum(abs(nullDiff) >= abs(obsDiff)) + 1) / (nPerm + 1);
fprintf('Permutation test (Late–Early mean drift): obs=%.4f, p=%.4g\n', obsDiff, p_perm);

%% ============================================================
%  (3) ADDED: Early vs Late drift bar plot (decoding-stage bar aesthetics)
%% ============================================================
% Mean +/- SEM with jittered grey scatter; significance from permutation p_perm
figBar = figure('Color','w','Position',[120,120,520,700]);
axB = axes(figBar); hold(axB,'on');

mE = mean(earlyVals,'omitnan');
mL = mean(lateVals,'omitnan');
semE = std(earlyVals,'omitnan') / sqrt(max(1,numel(earlyVals)));
semL = std(lateVals,'omitnan')  / sqrt(max(1,numel(lateVals)));

xPos = [1 2];
barW = 0.60;

b = bar(axB, xPos, [mE mL], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axesTickLW);
b.CData(1,:) = colEarlyShade;
b.CData(2,:) = colLateShade;
b.FaceAlpha  = 0.25;

errorbar(axB, xPos, [mE mL], [semE semL], 'k', 'LineStyle','none', 'LineWidth', axesTickLW, 'CapSize', 18);

if ~isempty(earlyVals)
    plot(axB, 1 + (rand(size(earlyVals))-0.5)*2*jit, earlyVals, 'o', 'MarkerSize', scatterMS, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', scatterEdgeLW);
end
if ~isempty(lateVals)
    plot(axB, 2 + (rand(size(lateVals))-0.5)*2*jit,  lateVals,  'o', 'MarkerSize', scatterMS, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', scatterEdgeLW);
end

set(axB, 'XLim',[0.5 2.5], 'XTick',[1 2], 'XTickLabel',{'Early','Late'});
ylabel(axB, 'Pearson r to Session 1 template', 'FontSize', labelFontSize);
title(axB, 'Early vs late similarity to Session 1', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axB, 'FontSize', tickFontSize, 'Box','off');
set(axB, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Add sig using permutation p-value
addSigIfNeeded_(axB, 1, 2, p_perm, [mE mL], axesTickLW, tickFontSize);

exportgraphics(figBar, outBarPng, 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figBar, outBarSvg); catch, print(figBar, outBarSvg, '-dsvg'); end
savefig(figBar, outBarFig);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseBar);

fprintf('\nDONE.\n');

%% ======================= HELPERS =======================

function [Y, eigvals] = classical_mds_(D, k)

    n = size(D,1);
    D = 0.5*(D + D');
    D(1:n+1:end) = 0;

    if exist('cmdscale','file') == 2
        try
            [Y, E] = cmdscale(D, k);
            eigvals = E(1:k);
            if size(Y,2) < k
                Y(:,end+1:k) = 0;
                eigvals(end+1:k) = 0;
            end
            return;
        catch
        end
    end

    J = eye(n) - (1/n)*ones(n);
    B = -0.5 * J * (D.^2) * J;
    B = 0.5*(B+B');

    [V, L] = eig(B);
    [lams, idx] = sort(diag(L), 'descend');
    lams = max(lams, 0);
    V = V(:, idx);

    eigvals = lams(1:k);
    Y = V(:,1:k) * diag(sqrt(eigvals + eps));
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
        yl = ylim(ax); yRange = yl(2)-yl(1);
        y = yMaxData + 0.08*yRange;
        h = 0.03*yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01*yRange, stars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', fs, 'FontWeight','bold', 'Color', [0 0 0]);
end

function s = pToStars_(p)
    if p < 0.001, s='***';
    elseif p < 0.01, s='**';
    else, s='*';
    end
end