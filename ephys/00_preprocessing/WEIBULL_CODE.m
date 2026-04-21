%% ========= PLOT_PERCENT_CORRECT_WEIBULL_FIT_ASYMPTOTE_STYLED__A_THESIS_FINAL.m =========
% Using the SAME behstruct file we've been using:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_FINAL_WITH_TRIALS.mat
%
% Produces the same figure behavior as your reference script:
%   - % correct per session from Hit (session.Hit or [session.trials.Hit])
%   - Weighted Weibull fit (weights = # trials used in that session)
%   - Prints B, A, lambda, k, T10/T50/T90
%   - Bootstraps CI for A; asymptote = first session with y >= A_lo
%   - Background shading: pre-asymptote green, post-asymptote yellow
%   - X ticks: label every 5th; manual ticks for all sessions
%   - Horizontal x tick labels placed lower; x-label moved further down
%
% Saves PNG to:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/2_FIGURE_PSTH/percent_correct_per_session_WEIBULL_asymptote.png

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/3_1_FIGURE_WEIBULL';
outPng  = fullfile(outDir, 'WEIBULL.png');

% (NEW) also save SVG
outSvg  = fullfile(outDir, 'WEIBULL.svg');

nBoot = 1000;
rng(0);

if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- STYLE ----
figPos = [120, 120, 900, 520];

dataLW      = 5;   % connecting line width
dataMS      = 10;  % marker size
markerLW    = 2.5; % marker edge width
fitLW       = 5;
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

%% ---- LOAD ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in MAT.');
nSess = numel(beh);

%% ---- Compute % correct per session from Hit ----
pctCorrect = nan(nSess,1);
nUsed      = zeros(nSess,1);

for sIdx = 1:nSess
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
    nUsed(sIdx) = numel(hit);
    if nUsed(sIdx)==0, continue; end

    pctCorrect(sIdx) = 100 * (sum(hit==1) / nUsed(sIdx));
end

%% ---- Prepare data for fitting ----
xAll = (1:nSess)';
ok   = isfinite(pctCorrect) & (nUsed > 0);

x = xAll(ok);
y = pctCorrect(ok);
w = double(nUsed(ok));

assert(numel(x) >= 4, 'Not enough sessions with Hit data to fit a Weibull.');

%% ---- Weibull model ----
weibullFun = @(p, s) p(2) - (p(2)-p(1)) .* exp( - (max(s,eps)./p(3)).^p(4) );
% p = [B, A, lambda, k]

%% ---- Init params ----
B0 = max(0,  min(100, median(y(1:min(3,end))) ));
A0 = max(0,  min(100, median(y(max(1,end-2):end)) ));
if A0 < B0 + 1
    A0 = min(100, B0 + max(5, std(y)));
end
lambda0 = max(1, round(median(x)));
k0      = 2;
p0 = [B0, A0, lambda0, k0];

LB = [0,   0,   0.1,  0.2];
UB = [100, 100, 1e3,  10];

obj = @(p) weightedSSE_withBounds_(p, x, y, w, weibullFun, LB, UB);

%% ---- Fit ----
opts = optimset('Display','off', 'MaxFunEvals', 5e4, 'MaxIter', 5e4);
pHat = fminsearch(obj, p0, opts);
pHat = min(UB, max(LB, pHat));

B      = pHat(1);
A      = pHat(2);
lambda = pHat(3);
k      = pHat(4);

%% ---- Latency-to-change summaries ----
T10 = lambda * (-log(1-0.10))^(1/k);
T50 = lambda * (-log(1-0.50))^(1/k);
T90 = lambda * (-log(1-0.90))^(1/k);

fprintf('\n==== WEIBULL FIT RESULTS ====\n');
fprintf('B=%.2f  A=%.2f  lambda=%.2f  k=%.2f\n', B, A, lambda, k);
fprintf('T10=%.2f  T50=%.2f  T90=%.2f\n', T10, T50, T90);

%% ---- Bootstrap 95% CI for A ----
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
A_hi = A_CI(2);

fprintf('A 95%% bootstrap CI = [%.2f, %.2f]\n', A_lo, A_hi);

%% ---- Asymptote reached criterion ----
asymSess = NaN;
for i = 1:numel(x)
    if y(i) >= A_lo
        asymSess = x(i);
        break
    end
end

if isfinite(asymSess)
    fprintf('Asymptote reached at Session %d (first y >= A_lo)\n', asymSess);
else
    fprintf('Asymptote not reached by criterion y >= A_lo\n');
end

%% ---- Plot ----
xFine = linspace(1, nSess, 600)';
yFit  = weibullFun(pHat, xFine);

fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

% -------------------- BACKGROUND SHADING --------------------
yl = [0 100];
set(ax, 'XLim',[1 nSess], 'YLim',yl);

if isfinite(asymSess)
    % pre (green): [1, asymSess]
    patch(ax, [1 asymSess asymSess 1], [yl(1) yl(1) yl(2) yl(2)], ...
        [0 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');

    % post (yellow): [asymSess, nSess]
    patch(ax, [asymSess nSess nSess asymSess], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');
else
    % If asymptote not reached, shade all as "pre" (green)
    patch(ax, [1 nSess nSess 1], [yl(1) yl(1) yl(2) yl(2)], ...
        [0 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');
end
% ------------------------------------------------------------

% --- DATA line ---
plot(ax, x, y, '-', 'LineWidth', dataLW, 'Color', [0 0.45 0.75]);

% --- OPEN circles ---
plot(ax, x, y, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', dataMS, ...
    'MarkerFaceColor', [0 0.45 0.75], ...
    'MarkerEdgeColor', [0 0.45 0.75], ...
    'LineWidth', markerLW);

% --- FIT ---
plot(ax, xFine, yFit, '-', 'LineWidth', fitLW, 'Color', [0.85 0.33 0.10]);

% --- CRITERIA LINES ---
% (CHANGED) removed horizontal dashed line at A_lo
% yline(ax, A_lo, '--', 'LineWidth', eventLineLW, 'Color', [0.3 0.3 0.3]);
if isfinite(asymSess)
    xline(ax, asymSess, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);
end

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, '% correct', 'FontSize', labelFontSize);
title(ax, 'Weibull function fit to learning curve', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSess], 'YLim',[0 100]);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

%% ---- TICKS: label every 5th; draw ticks + labels manually ----
ax.XTick = 1:nSess;
labIdx = 5:5:nSess;

% blank built-in labels; we draw them manually
ax.XTickLabel = repmat({''}, 1, nSess);
xtickangle(ax, 0);

% Hide axis' own tick marks
set(ax, 'TickLength', [0 0]);

% Manual ticks (DATA units)
y0 = ax.YLim(1);
yr = diff(ax.YLim);

majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

tickCol = [0 0 0];

for s = 1:nSess
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
set(ax, 'XLim',[1 nSess], 'YLim',[0 100]);

saveas(fig, outPng);

% (NEW) also save SVG
saveas(fig, outSvg);

fprintf('Saved: %s\n', outPng);
fprintf('Saved: %s\n', outSvg);

%% ================= HELPERS =================

function sse = weightedSSE_withBounds_(p, x, y, w, f, LB, UB)
    pen = 0;
    for i=1:numel(p)
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