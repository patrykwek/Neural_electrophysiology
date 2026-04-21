%% ========= FIG4D_LIKE_SESSIONPAIR_UNITS__EARLY_LATE__GLOBALWARP__CLASSIFIEDONLY__PVALUE.m =========
% Fig 4d-like (analog) for your dataset using SESSION-PAIR units (matrix entries),
% matching the logic of code 2 (Fig2d):
%
%   1) Build per-session ensemble vectors (ens_session(:,s))
%   2) Compute session x session similarity matrix: SIM = corr(ens_session, 'Rows','pairwise')
%   3) Treat each session-pair correlation as a unit:
%        - Early–Early (EE): SIM(i,j) for i<j in earlyValid
%        - Late–Late  (LL): SIM(i,j) for i<j in lateValid
%        - Early–Late (EL): SIM(i,j) for i in earlyValid, j in lateValid
%   4) Stats: Wilcoxon rank-sum (two-sided) comparing distributions
%
% LEFT: histogram overlay (NO p-values shown). X axis starts at 0.5.
% RIGHT: bar+jitter style:
%        bar = mean(r) across session-pairs; errorbar = SEM across session-pairs;
%        jittered points = session-pair correlations (gray). Y axis starts at 0.5.
%        Pairwise p comparisons drawn as significance bars + stars ONLY if significant
%        using the SAME star rule as your FIG3E code.
%
% ONLY CHANGE NOW:
%   - Report p values and test data (names etc) in the command window.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/2_FIGURE_NEURAL_SIMILARITY';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Early/Late definition (match your FIG3A script)
stageRanges = [1 8; 29 36];  % [earlyStart earlyEnd; lateStart lateEnd]

% Heatmap window definition (relative to cue) — display design
preCueMs   = 500;
postLickMs = 3000;

% FIG2-like trial-selection window (relative to cue) — drives per-unit inclusion windows
fixedWin = [-1000, 6000];   % ms relative to cue

% Trial filters
MIN_RT_MS = 100;
minTrialsPerUnit = 10;
requireValid = true;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Plot
nBins = 30;

% Style
figPos = [200 200 1180 520];

% ---- ORIGINAL STYLE (kept, then scaled 2x below) ----
titleFontSize = 24;
labelFontSize = 22;
tickFontSize  = 18;
axesTickLW    = 2.5;
axesTickLen   = [0.02 0.02];
histAxesTickLW = 4.0;

% line widths used in plot elements
histLineLW   = 3;
errLineLW    = 2;
sigLineLW    = 2;

% scatter size
ptSize = 16;

% ---- SCALE FACTOR (2x) ----
SCALE = 2;

titleFontSize = titleFontSize * SCALE;
labelFontSize = labelFontSize * SCALE;
tickFontSize  = tickFontSize  * SCALE;

axesTickLW     = axesTickLW     * SCALE;
histAxesTickLW = histAxesTickLW * SCALE;

histLineLW = histLineLW * SCALE;
errLineLW  = errLineLW  * SCALE;
sigLineLW  = sigLineLW  * SCALE;

ptSize = ptSize * SCALE;

% KEEP YOUR COLORS (from your screenshot): EE blue, LL red/orange, EL yellow/orange
colEE = [0.10 0.55 0.95];
colLL = [0.90 0.20 0.10];
colEL = [0.95 0.70 0.10];

tag = sprintf('FIG4Dlike_SESSIONPAIRUNITS_EE_LL_EL_RANKSUM_Early%02d_%02d_Late%02d_%02d', ...
    stageRanges(1,1),stageRanges(1,2),stageRanges(2,1),stageRanges(2,2));

%% ---- LOAD BEHSTRUCT ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = numel(beh);

%% ---- LOAD CELL TYPES ----
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

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

        lr = lick - cue;

        all_pressLat_global(end+1,1) = rt; %#ok<AGROW>
        all_lickLat_global(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat_global) && ~isempty(all_lickLat_global), ...
    'No trials passed global filters (cannot set shared warp targets).');

Tpress_global = median(all_pressLat_global);
Tlick_global  = median(all_lickLat_global);

fprintf('\n=== GLOBAL targets: Tpress=%.1f ms, Tlick=%.1f ms ===\n', Tpress_global, Tlick_global);

%% ---- TEMPLATE AXIS (shared) ----
winLeft  = -preCueMs;
winRight_template = Tlick_global + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,             winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress_global, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick_global,  winLeft, dt_ms, nT);

%% ---- SESSION GROUPS ----
e0 = max(1, min(nSessions, stageRanges(1,1)));
e1 = max(1, min(nSessions, stageRanges(1,2)));
l0 = max(1, min(nSessions, stageRanges(2,1)));
l1 = max(1, min(nSessions, stageRanges(2,2)));

earlySess = e0:e1;
lateSess  = l0:l1;

fprintf('Early sessions: %d-%d | Late sessions: %d-%d\n', e0,e1,l0,l1);

%% ---- BUILD PER-SESSION ENSEMBLE VECTORS (nT x nSessions) ----
ens_session = nan(nT, nSessions);
ens_nUnits  = zeros(1, nSessions);

for sIdx = unique([earlySess, lateSess])

    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Behavior-only keep trials (EXACT rules from your FIG3A code) ----
    keepTrial = false(numel(session.trials),1);
    cueAbs   = nan(numel(session.trials),1);
    pressLat = nan(numel(session.trials),1);
    lickLat  = nan(numel(session.trials),1);

    for k = 1:numel(session.trials)
        tr = session.trials(k);

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

        lr = lick - cue;

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

    unitEvenMeans = []; % nUnitsThisSession x nT

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

            % PER-UNIT trial inclusion: spike in [cue + fixedWin(1), cue + lickLat]
            activeTrials = false(numel(idxKeep),1);
            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                lr0  = lickLatK(iTr);

                w0 = cue0 + fixedWin(1);
                w1 = cue0 + lr0;

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

            % Held-out mean = EVEN trials
            idxEven = 2:2:nTrU;
            if isempty(idxEven), idxEven = 1:nTrU; end
            muEven = mean(warpedZ(idxEven,:), 1);

            unitEvenMeans(end+1,:) = muEven; %#ok<AGROW>
        end
    end

    if ~isempty(unitEvenMeans)
        ens = mean(unitEvenMeans, 1);
        ensZ = zscoreTrial_(ens);
        ens_session(:, sIdx) = ensZ(:);
        ens_nUnits(sIdx) = size(unitEvenMeans,1);
    end
end

%% ---- Keep only sessions with valid ensemble vectors ----
earlyValid = earlySess(~any(isnan(ens_session(:,earlySess)),1));
lateValid  = lateSess(~any(isnan(ens_session(:,lateSess)),1));

assert(~isempty(earlyValid) && ~isempty(lateValid), ...
    'No valid early or late sessions after ensemble extraction.');

%% ---- SESSION-PAIR UNITS FROM SIMILARITY MATRIX ----
validSess = unique([earlyValid(:); lateValid(:)])';
SIM = corr(ens_session(:, validSess), 'Rows','pairwise');

[~, idxEarlySIM] = ismember(earlyValid, validSess);
[~, idxLateSIM]  = ismember(lateValid,  validSess);
assert(all(idxEarlySIM>0) && all(idxLateSIM>0), 'Session index mapping failed.');

pairsEE = nchoosek(idxEarlySIM, 2);
r_EE = nan(size(pairsEE,1),1);
for k = 1:size(pairsEE,1)
    r_EE(k) = SIM(pairsEE(k,1), pairsEE(k,2));
end
r_EE = r_EE(isfinite(r_EE));

pairsLL = nchoosek(idxLateSIM, 2);
r_LL = nan(size(pairsLL,1),1);
for k = 1:size(pairsLL,1)
    r_LL(k) = SIM(pairsLL(k,1), pairsLL(k,2));
end
r_LL = r_LL(isfinite(r_LL));

[gridE, gridL] = ndgrid(idxEarlySIM, idxLateSIM);
r_EL = SIM(sub2ind(size(SIM), gridE(:), gridL(:)));
r_EL = r_EL(isfinite(r_EL));

%% ---- STATS: WILCOXON RANK-SUM (TWO-SIDED) ----
[p2_LL_EE,~,~] = ranksum(r_LL, r_EE);
[p2_LL_EL,~,~] = ranksum(r_LL, r_EL);
[p2_EE_EL,~,~] = ranksum(r_EE, r_EL);

%% ---- ONLY CHANGE: REPORT TESTS + P VALUES IN COMMAND WINDOW ----
fprintf('\n=== Rank-sum tests on session-pair correlations (two-sided) ===\n');
fprintf('Groups:\n');
fprintf('  EE (Early–Early): sessions %d-%d, valid=%s | nPairs=%d\n', e0,e1, mat2str(earlyValid), numel(r_EE));
fprintf('  LL (Late–Late)  : sessions %d-%d, valid=%s | nPairs=%d\n', l0,l1, mat2str(lateValid),  numel(r_LL));
fprintf('  EL (Early–Late) : Early valid=%s vs Late valid=%s | nPairs=%d\n', mat2str(earlyValid), mat2str(lateValid), numel(r_EL));
fprintf('\nTests (ranksum):\n');
fprintf('  LL vs EE : p = %.6g  (nLL=%d, nEE=%d)\n', p2_LL_EE, numel(r_LL), numel(r_EE));
fprintf('  LL vs EL : p = %.6g  (nLL=%d, nEL=%d)\n', p2_LL_EL, numel(r_LL), numel(r_EL));
fprintf('  EE vs EL : p = %.6g  (nEE=%d, nEL=%d)\n', p2_EE_EL, numel(r_EE), numel(r_EL));
fprintf('=============================================================\n\n');

%% ---- SAVE RESULTS (.mat) ----
outMat = fullfile(outDir, [tag '_results.mat']);
save(outMat, ...
    'r_EE','r_LL','r_EL', ...
    'p2_LL_EE','p2_LL_EL','p2_EE_EL', ...
    'SIM','validSess','earlySess','lateSess','earlyValid','lateValid', ...
    'ens_session','ens_nUnits', ...
    'tgrid_ms','Tpress_global','Tlick_global', ...
    'stageRanges','fixedWin','MIN_RT_MS','minTrialsPerUnit','requireValid', ...
    'dt_ms','gaussSigmaMs','preCueMs','postLickMs', ...
    '-v7.3');

%% ---- PLOT: Left histogram + Right bar( mean±SEM ) + jittered points + stars ----
fig = figure('Color','w','Position',figPos);
tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% ---------- LEFT: Histogram overlay (NO p-values; x starts at 0.5) ----------
ax1 = nexttile(tl, 1); hold(ax1,'on');

edges = linspace(0, 1, nBins+1);

hEE = histogram(ax1, r_EE, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);
hLL = histogram(ax1, r_LL, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);
hEL = histogram(ax1, r_EL, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);

hEE.EdgeColor = colEE;
hLL.EdgeColor = colLL;
hEL.EdgeColor = colEL;

xlabel(ax1, 'Correlation r (session-pair units from SIM matrix)', 'FontSize', labelFontSize);
ylabel(ax1, 'Probability', 'FontSize', labelFontSize);
title(ax1, 'Session-pair similarity (EE/LL/EL)', 'FontSize', titleFontSize, 'FontWeight','bold');

legend(ax1, {sprintf('Early–Early (sessions %d-%d)', e0,e1), ...
             sprintf('Late–Late (sessions %d-%d)',  l0,l1), ...
             'Early–Late'}, ...
    'Location','northwest');

xlim(ax1, [0.5 1]);

set(ax1, 'FontSize', tickFontSize, ...
    'TickDir','out', ...
    'LineWidth', histAxesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box','off');

% ---------- RIGHT: Bar + SEM + session-pair points ----------
ax2 = nexttile(tl, 2); hold(ax2,'on');

cats = ["Early–Early","Late–Late","Early–Late"];
xCats = categorical(cats, cats, 'Ordinal', true);

R = {r_EE, r_LL, r_EL};
cols = {colEE, colLL, colEL};

mu  = [mean(r_EE,'omitnan'), mean(r_LL,'omitnan'), mean(r_EL,'omitnan')];
nEE = numel(r_EE); nLL = numel(r_LL); nEL = numel(r_EL);
sem = [std(r_EE,0,'omitnan')/sqrt(max(1,nEE)), ...
       std(r_LL,0,'omitnan')/sqrt(max(1,nLL)), ...
       std(r_EL,0,'omitnan')/sqrt(max(1,nEL))];

b = bar(ax2, xCats, mu);
b.FaceColor = 'flat';
b.CData(1,:) = cols{1};
b.CData(2,:) = cols{2};
b.CData(3,:) = cols{3};

errorbar(ax2, b.XEndPoints, mu, sem, 'k.', 'LineWidth', errLineLW);

rng(0);
jitter = 0.12;

nShow = [min(nEE,3000), min(nLL,3000), min(nEL,3000)];
idxShow = cell(1,3);
idxShow{1} = randperm(nEE, nShow(1));
idxShow{2} = randperm(nLL, nShow(2));
idxShow{3} = randperm(nEL, nShow(3));

for i = 1:3
    rr = R{i}(idxShow{i});
    xj = b.XEndPoints(i) + (rand(numel(rr),1)-0.5)*2*jitter;
    scatter(ax2, xj, rr, ptSize, 'filled', ...
        'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerFaceAlpha', 0.20, ...
        'MarkerEdgeAlpha', 0);
end

ylabel(ax2, 'Correlation r', 'FontSize', labelFontSize);
xlabel(ax2, 'Comparison (Early–Early / Late–Late / Early–Late)', 'FontSize', labelFontSize);

title(ax2, 'Pairwise tests (rank-sum p) with stars', ...
    'FontSize', titleFontSize, 'FontWeight','bold');

ylim(ax2, [0.5 1]);
set(ax2, 'FontSize', tickFontSize, 'LineWidth', axesTickLW, 'TickDir','out');
box(ax2,'off');

%% ---- Significance bars + stars (ONLY if significant; FIG3E rule) ----
pairs = [1 2; 2 3; 1 3];                 % (EE-LL), (LL-EL), (EE-EL)
pVals = [p2_LL_EE; p2_LL_EL; p2_EE_EL];

sig = pVals < 0.05;

if any(sig)
    yMax = max([r_EE(:); r_LL(:); r_EL(:)], [], 'omitnan');
    yMin = min([r_EE(:); r_LL(:); r_EL(:)], [], 'omitnan');
    yRange = yMax - yMin;
    if ~isfinite(yRange) || yRange == 0, yRange = 1; end

    baseY = max(mu + sem) + 0.06*yRange;
    hBar  = 0.02*yRange;
    stepY = 0.06*yRange;

    order = [1 2 3];
    drawCount = 0;

    for kk = 1:3
        k = order(kk);
        if ~sig(k), continue; end
        drawCount = drawCount + 1;

        xA = b.XEndPoints(pairs(k,1));
        xB = b.XEndPoints(pairs(k,2));
        yLevel = baseY + (drawCount-1)*stepY;

        plot(ax2, [xA xA xB xB], [yLevel yLevel+hBar yLevel+hBar yLevel], 'k-', 'LineWidth', sigLineLW);

        if pVals(k) < 0.001
            s = '***';
        elseif pVals(k) < 0.01
            s = '**';
        elseif pVals(k) < 0.05
            s = '*';
        else
            s = '';
        end

        text(ax2, mean([xA xB]), yLevel+hBar+0.002*yRange, s, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', titleFontSize, 'FontWeight','bold');
    end

    yTopNeed = baseY + drawCount*stepY + 0.06*yRange;
    ylim(ax2, [0.5, max(1, yTopNeed)]);
end

%% ---- SAVE FIGURES (PNG + FIG + SVG) ----
outPng = fullfile(outDir, [tag '_hist_plus_barjitter.png']);
outFig = fullfile(outDir, [tag '_hist_plus_barjitter.fig']);
outSvg = fullfile(outDir, [tag '_hist_plus_barjitter.svg']);

saveas(fig, outPng);
savefig(fig, outFig);
try
    print(fig, outSvg, '-dsvg');
catch
    saveas(fig, outSvg);
end

fprintf('Saved PNG: %s\n', outPng);
fprintf('Saved FIG: %s\n', outFig);
fprintf('Saved SVG: %s\n', outSvg);
fprintf('DONE. Outputs in: %s\n', outDir);

%% ================= HELPERS (same as your FIG3A script) =================

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
        spikes = S.spikes; % last resort
    end
end