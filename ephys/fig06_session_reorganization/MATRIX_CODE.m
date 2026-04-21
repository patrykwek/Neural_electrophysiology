%% ========= FIG3A_SESSION_ENSEMBLE_SIMILARITYMATRIX__A_THESIS_FINAL__GLOBALWARP__FIG2TRIALS_EXACT__CLASSIFIEDONLY.m =========
% SESSION x SESSION similarity matrix (sessions 1–36 on BOTH axes).
%
% Uses the SAME pipeline as your session-ensemble heatmap script to compute
% a per-session warped trial activity vector:
%   - for each session: average (held-out even-trial mean) warped z-SDF per unit,
%     then ensemble-average across units, then z-score across time => ensZ (1 x nT)
%
% Then computes session-by-session similarity as Pearson correlation between
% session ensemble vectors across time (pairwise rows, to tolerate NaNs).
%
% Requested changes (ONLY):
%   - Make each session row/column crisp: NO smoothing/blur between sessions
%   - Use a better paper colormap for correlation (diverging, centered at 0)
%   - Remove mid-cell grey lines; add stronger boundaries between sessions only
%   - Improve colorbar
%   - Ticks like your heatmap: minor ticks for all, major every 5th (bigger/labeled)
%   - Use the SAME color scale convention as plotHeat_ in the provided code:
%       caxis([-1 1]) and colorbar ticks in that range (no other changes)
%   - USE THE SAME COLORMAP AS YOUR DMS HEATMAPS: parula(256)
%   - Save SVG in addition to PNG (requested)
%
% Added (paper-style stats analog; printed to command window):
%   - Compare within-group session-pair correlations using Wilcoxon rank-sum:
%       Early group = sessions 1–8
%       Late  group = sessions 29–36
%
% Visual change requested NOW:
%   - Make ONLY the boundary BETWEEN sessions 8 and 9 (i.e., at 8.5) a
%     wider dashed grey "striped" line (like xline(...,'--') in your Weibull script).
%   - Do NOT add any other frames or highlights.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/3_4_FIGURE_MATRIX';
outDir  = fullfile(baseOut, 'MATRIX');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Heatmap window definition (relative to cue) — display design
preCueMs   = 500;     % show -500 ms to cue
postLickMs = 3000;    % show +3000 ms after lick

% FIG2-like trial-selection window (relative to cue) — drives trial inclusion filters
fixedWin = [-1000, 6000];   % ms relative to cue

% Trial filters
MIN_RT_MS = 100;
minTrialsPerUnit = 10;
requireValid = true;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Display
xUnit = "s";   % "s" or "ms"
forceIntegerSecondsTicks = true;

% ---- FIGURE SIZE ----
figPos = [120, 120, 980, 900];

% ---- COLORBAR formatting ----
cbFontSize   = 22;
cbHeightFrac = 0.50;     % longer colorbar (like your heatmap script)
cbWidthFrac  = 0.040;    % thinner
cbGapFrac    = 0.07;     % gap from axes (fraction of ax width)

% ---- STYLE ----
titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

axesTickLW    = 4.0;
axesTickLen   = [0.02 0.02];

% ============================================================
% SIMILARITY PLOT SETTINGS
% ============================================================
simInterp = 'nearest';     % crisp cells
simCLim   = [-1 1];        % (match plotHeat_ scale)
simCmapN  = 256;           % colormap resolution

% Boundary line styling (between sessions only)
boundaryLW    = 1.2;       % stronger than before
boundaryColor = [0 0 0];   % black
boundaryAlpha = 0.35;      % visible but not overwhelming

% ============================================================
% PAPER-STYLE GROUPS FOR STABILITY STATISTICS (your mapping)
% ============================================================
earlyGroup = 1:8;
lateGroup  = 29:36;

% ============================================================
% Special dashed boundary between sessions 8 and 9 (at 8.5)
% (match Weibull xline style: '--', darker grey, wider)
% ============================================================
dashBoundaryAt = 8.5;
dashColor = [0.2 0.2 0.2];
dashLW    = 5;              % wider
dashStyle = '--';

%% ---- LOAD BEHSTRUCT ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = numel(beh);

% ============================================================
% ONLY USE SESSIONS 1–36
% ============================================================
nSessions = min(nSessions, 36);
beh = beh(1:nSessions);

%% ---- LOAD CELL TYPES ----
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u} : 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skip)

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (shared) ----
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

        pr = press - cue;
        lr = lick  - cue;

        if ~(pr >= fixedWin(1) && pr <= fixedWin(2)), continue; end
        if ~(lr >= fixedWin(1) && lr <= fixedWin(2)), continue; end

        all_pressLat_global(end+1,1) = rt; %#ok<AGROW>
        all_lickLat_global(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat_global) && ~isempty(all_lickLat_global), ...
    'No trials passed global FIG2-like filters (cannot set shared warp targets).');

Tpress = median(all_pressLat_global);
Tlick  = median(all_lickLat_global);

fprintf('\n=== GLOBAL targets: Tpress=%.1f ms, Tlick=%.1f ms ===\n', Tpress, Tlick);

%% ---- TEMPLATE AXIS FOR ENSEMBLE (display design) ----
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

%% ---- BUILD SESSION ENSEMBLE MATRIX (time x session) ----
PETH = nan(nT, nSessions); % nT x nSessions

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Behavior-only keep trials (EXACT FIG2-like behavioral inclusion; NO spikes here) ----
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

        pr = press - cue;
        lr = lick  - cue;

        if ~(pr >= fixedWin(1) && pr <= fixedWin(2)), continue; end
        if ~(lr >= fixedWin(1) && lr <= fixedWin(2)), continue; end

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

    unitCurves = [];  % nUnitsInSession x nT

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)

            % ONLY KEEP CLASSIFIED UNITS (1/2/3)
            ct = [];
            if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
                    ch <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
                    u <= numel(cell_type{sIdx}{ch})
                ct = cell_type{sIdx}{ch}{u};
            end
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                continue;
            end

            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
            spk_abs = double(spk_abs(:));

            % PER-UNIT trial inclusion: require spike in [cue+fixedWin(1), cue+fixedWin(2)]
            activeTrials = false(numel(idxKeep),1);
            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                w0 = cue0 + fixedWin(1);
                w1 = cue0 + fixedWin(2);
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

            warpedZ = zeros(nTrU, nT);

            for iTr = 1:nTrU
                cue0   = cueU(iTr);
                tPress = rtU(iTr);
                tLick  = lkU(iTr);

                winRight_trial = tLick + postLickMs;
                tgrid_trial = winLeft:dt_ms:winRight_trial;
                nT_trial = numel(tgrid_trial);

                idx0 = timeToIdx_(0,      winLeft, dt_ms, nT_trial);
                idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

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
                y  = sm / (dt_ms/1000);  % Hz

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
                    segC_unwarped = [segC, zeros(1, kC_target - numel(segC))];
                else
                    segC_unwarped = segC(1:kC_target);
                end

                yy = [pre, wA(2:end), wB(2:end), segC_unwarped(2:end)];

                if numel(yy) < nT
                    yy = [yy, zeros(1, nT-numel(yy))];
                elseif numel(yy) > nT
                    yy = yy(1:nT);
                end

                warpedZ(iTr,:) = zscoreTrial_(yy);
            end

            idxPlot = 2:2:nTrU;
            if isempty(idxPlot), idxPlot = 1:nTrU; end
            muPlot = mean(warpedZ(idxPlot,:), 1);

            unitCurves(end+1, :) = muPlot; %#ok<AGROW>
        end
    end

    if isempty(unitCurves), continue; end

    ens  = mean(unitCurves, 1);
    ensZ = zscoreTrial_(ens);

    PETH(:, sIdx) = ensZ(:);
end

%% ---- SESSION x SESSION similarity matrix (corr across time) ----
SIM = corr(PETH, 'Rows','pairwise'); % nSessions x nSessions

%% ============================================================
% PAPER-STYLE STABILITY STATISTICS ON YOUR SIM MATRIX
%   Compare within-group session-pair correlations:
%     Early: sessions 1–8
%     Late : sessions 29–36
%   Test: Wilcoxon rank-sum (two-sided), print to command window.
%% ============================================================

% safety: keep groups within [1..nSessions]
earlyGroup = earlyGroup(earlyGroup >= 1 & earlyGroup <= nSessions);
lateGroup  = lateGroup( lateGroup >= 1 & lateGroup  <= nSessions);

% sessions that actually have ensemble vectors (some timepoints non-NaN)
validSess = find(~all(isnan(PETH),1));

earlyUse  = intersect(earlyGroup, validSess);
lateUse   = intersect(lateGroup,  validSess);

assert(numel(earlyUse) >= 2 && numel(lateUse) >= 2, ...
    'Not enough valid sessions in early or late group to compute stability statistics.');

earlyPairs = nchoosek(earlyUse, 2);
latePairs  = nchoosek(lateUse,  2);

rEarly = nan(size(earlyPairs,1),1);
for k = 1:size(earlyPairs,1)
    i = earlyPairs(k,1); j = earlyPairs(k,2);
    rEarly(k) = SIM(i,j);
end
rEarly = rEarly(isfinite(rEarly));

rLate = nan(size(latePairs,1),1);
for k = 1:size(latePairs,1)
    i = latePairs(k,1); j = latePairs(k,2);
    rLate(k) = SIM(i,j);
end
rLate = rLate(isfinite(rLate));

[p_RS,~,stats_RS] = ranksum(rEarly, rLate); % Wilcoxon rank-sum, two-sided

fprintf('\n=== PAPER-STYLE STABILITY (within-group session-pairs) ===\n');
fprintf('Early group requested: 1–8 | using valid sessions: %s\n', mat2str(earlyUse));
fprintf('Late  group requested: 29–36 | using valid sessions: %s\n', mat2str(lateUse));
fprintf('nPairs: early=%d, late=%d\n', numel(rEarly), numel(rLate));
fprintf('Median r: early=%.4f | late=%.4f\n', median(rEarly,'omitnan'), median(rLate,'omitnan'));
fprintf('Mean±SD r: early=%.4f±%.4f | late=%.4f±%.4f\n', ...
    mean(rEarly,'omitnan'), std(rEarly,'omitnan'), mean(rLate,'omitnan'), std(rLate,'omitnan'));
fprintf('Wilcoxon rank-sum (two-sided): p=%.4g, z=%.3f\n', p_RS, stats_RS.zval);

%% ---- Save matrices + stats ----
outmat = fullfile(outDir, sprintf('Sessions01_%02d_sessionEnsemble_similarity_GLOBALWARP.mat', nSessions));
save(outmat, 'PETH', 'SIM', 'tgrid_ms', 'Tpress', 'Tlick', ...
    'earlyGroup','lateGroup','earlyUse','lateUse','rEarly','rLate','p_RS','stats_RS', '-v7.3');
fprintf('Saved: %s\n', outmat);

%% ---- Plot similarity heatmap (sessions x sessions) ----
fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

hImg = imagesc(ax, 1:nSessions, 1:nSessions, SIM);

% Session 1 on TOP, Session N on BOTTOM
set(ax,'YDir','reverse');

% crisp pixels
try
    set(hImg, 'Interpolation', simInterp);
catch
    set(hImg, 'Interpolation', 'nearest');
end

% USE PLASMA COLORMAP (drop-in; provided below)
colormap(ax, plasma(simCmapN));

% ---- USE SAME SCALE CONVENTION AS plotHeat_ ----
caxis(ax, simCLim);   % [-1 1]

% square cells
axis(ax,'image');
xlim(ax,[0.5 nSessions+0.5]);
ylim(ax,[0.5 nSessions+0.5]);

xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Session', 'FontSize', labelFontSize);

title(ax, 'Pearson correlation (r)', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

%% ---- TICKS like your heatmap: major every 5, minor all others ----
major = 1:5:nSessions;
if major(end) ~= nSessions
    major = unique([major nSessions]);
end

xticks(ax, major);
yticks(ax, major);
xticklabels(ax, string(major));
yticklabels(ax, string(major));

ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = setdiff(1:nSessions, major);
ax.YAxis.MinorTickValues = setdiff(1:nSessions, major);

set(ax, 'FontSize', tickFontSize, ...
    'TickDir','out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box','off', ...
    'Layer','top');

%% ---- REMOVE GRID; draw ONLY boundaries BETWEEN sessions ----
ax.XGrid = 'off';
ax.YGrid = 'off';

for b = 1.5:1:(nSessions-0.5)

    % default boundaries
    thisLW    = boundaryLW;
    thisCol   = boundaryColor;
    thisStyle = '-';
    thisAlpha = boundaryAlpha;

    % ONLY change the boundary between sessions 8 and 9 (at 8.5)
    if abs(b - dashBoundaryAt) < 1e-9
        thisLW    = dashLW;
        thisCol   = dashColor;
        thisStyle = dashStyle;
        thisAlpha = 1.0; % keep dashed line crisp
    end

    lh = line(ax, [0.5 nSessions+0.5], [b b], 'Color', thisCol, 'LineWidth', thisLW, 'LineStyle', thisStyle);
    lv = line(ax, [b b], [0.5 nSessions+0.5], 'Color', thisCol, 'LineWidth', thisLW, 'LineStyle', thisStyle);

    % keep alpha behavior for non-dashed boundaries
    if ~(abs(b - dashBoundaryAt) < 1e-9)
        try
            lh.Color(4) = thisAlpha;
            lv.Color(4) = thisAlpha;
        catch
        end
    end
end

%% ---- COLORBAR: match plotHeat_ tick style (only scale; keep your placement) ----
cb = colorbar(ax, 'eastoutside');
cb.TickDirection = 'out';
cb.LineWidth = axesTickLW;
cb.Label.String = 'Pearson r';
cb.FontSize = cbFontSize;
cb.Label.FontSize = cbFontSize;
cb.Ticks = [-1 -0.5 0 0.5 1];

cb.Units = 'normalized';

axpos = ax.Position;
cbpos = cb.Position;

cbpos(3) = cbWidthFrac;
cbpos(4) = cbHeightFrac * axpos(4);
cbpos(2) = axpos(2) + (axpos(4)-cbpos(4))/2;
cbpos(1) = axpos(1) + axpos(3) + cbGapFrac*axpos(3);
cbpos(1) = min(cbpos(1), 0.92);

cb.Position = cbpos;

%% ---- SAVE (PNG + SVG) ----
outbase = fullfile(outDir, sprintf('Sessions01_%02d_similarityMatrix_GLOBALWARP', nSessions));
set(fig, 'PaperPositionMode', 'auto');

exportgraphics(fig, [outbase '.png'], 'Resolution', 300, 'BackgroundColor', 'white');

try
    saveas(fig, [outbase '.svg']);
catch
    print(fig, [outbase '.svg'], '-dsvg');
end

fprintf('Saved: %s.png\n', outbase);
fprintf('Saved: %s.svg\n', outbase);

fprintf('\nDONE. Similarity matrix saved in: %s\n', outDir);

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
        spikes = S.spikes;
    end
end

% --- DROP-IN PLASMA COLORMAP (no toolbox required) ---
function cmap = plasma(m)
%PLASMA  Approximation of matplotlib's "plasma" colormap.
%   CMAP = PLASMA(M) returns an Mx3 colormap.
%   If M is omitted, uses current figure colormap size.

    if nargin < 1 || isempty(m)
        m = size(get(gcf,'colormap'),1);
    end

    % 17-point anchor table (RGB in [0,1]) approximating "plasma"
    % (compact, self-contained; interpolated to m points)
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

    % ensure bounds
    cmap = max(0, min(1, cmap));
end