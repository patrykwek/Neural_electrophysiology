%% ========= PLOT_PERCENT_CORRECT_LEARNING_CURVE_STYLED__A_THESIS_FINAL.m =========
% Using the SAME behstruct file we've been using:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_FINAL_WITH_TRIALS.mat
%
% Produces:
%   - % correct per session from Hit (session.Hit or [session.trials.Hit])
%   - Background shading by stage boundaries:
%       * Sessions 1..38      = yellow (same shade as your current post-asymptote yellow)
%       * Sessions 38..41     = orange-ish
%       * Sessions 41..end    = different orange shade
%   - Vertical dashed lines at session 38 and 41
%   - X ticks: label every 10th session; long tick marks at labeled sessions, short elsewhere
%
% Saves PNG to:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/WEIBULL/percent_correct_per_session_WEIBULL_asymptote.png

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_TESTFINAL_WITH_TRIALS.mat';
outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/BEHAVIOR_EPHYS';
outPng  = fullfile(outDir, 'behavior.png');

if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- STYLE ----
figPos = [120, 120, 900, 520];

dataLW      = 5;   % connecting line width
dataMS      = 10;  % marker size
markerLW    = 2.5; % marker edge width
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

%% ---- STAGE BOUNDARIES (vertical lines + background shading) ----
stage1_end = 38;   % sessions up to 38
stage2_end = 41;   % sessions 38..41
% sessions >41 are stage 3

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

%% ---- Prepare data (no Weibull; just plot observed sessions) ----
xAll = (1:nSess)';
ok   = isfinite(pctCorrect) & (nUsed > 0);

x = xAll(ok);
y = pctCorrect(ok);

%% ---- Plot ----
fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

yl = [0 100];
set(ax, 'XLim',[1 nSess], 'YLim',yl);

% -------------------- BACKGROUND SHADING BY STAGE --------------------
% Colors:
%   - Stage 1 uses the SAME yellow shade you currently used ([1 1 0])
%   - Stage 2 orange-ish
%   - Stage 3 different orange shade
c_stage1 = [1 1 0];          % original yellow
c_stage2 = [1.00 0.70 0.25]; % orange-ish
c_stage3 = [1.00 0.45 0.10]; % deeper/different orange

% Stage 1: [1, stage1_end]
x1a = 1;
x1b = min(stage1_end, nSess);
if x1b >= x1a
    patch(ax, [x1a x1b x1b x1a], [yl(1) yl(1) yl(2) yl(2)], ...
        c_stage1, 'FaceAlpha', 0.10, 'EdgeColor','none');
end

% Stage 2: [stage1_end, stage2_end]
x2a = min(stage1_end, nSess);
x2b = min(stage2_end, nSess);
if x2b >= x2a
    patch(ax, [x2a x2b x2b x2a], [yl(1) yl(1) yl(2) yl(2)], ...
        c_stage2, 'FaceAlpha', 0.10, 'EdgeColor','none');
end

% Stage 3: [stage2_end, nSess]
x3a = min(stage2_end, nSess);
x3b = nSess;
if x3b >= x3a
    patch(ax, [x3a x3b x3b x3a], [yl(1) yl(1) yl(2) yl(2)], ...
        c_stage3, 'FaceAlpha', 0.10, 'EdgeColor','none');
end
% -------------------------------------------------------------------

% --- DATA line ---
plot(ax, x, y, '-', 'LineWidth', dataLW, 'Color', [0 0.45 0.75]);

% --- FILLED circles ---
plot(ax, x, y, 'o', ...
    'LineStyle','none', ...
    'MarkerSize', dataMS, ...
    'MarkerFaceColor', [0 0.45 0.75], ...
    'MarkerEdgeColor', [0 0.45 0.75], ...
    'LineWidth', markerLW);

% --- VERTICAL STAGE LINES ---
if stage1_end <= nSess
    xline(ax, stage1_end, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);
end
if stage2_end <= nSess
    xline(ax, stage2_end, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);
end

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, '% correct', 'FontSize', labelFontSize);
title(ax, 'Learning curve per stage and session', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSess], 'YLim',[0 100]);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

%% ---- TICKS: label every 10th; draw ticks + labels manually ----
ax.XTick = 1:nSess;
labIdx = 10:10:nSess;

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
fprintf('Saved: %s\n', outPng);

%% ================= HELPERS =================

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