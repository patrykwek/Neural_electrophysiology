%% ========= FIG3X_SESSION_CUETYPE_DECODING_HEATMAP__A_THESIS_FINAL__GLOBALWARP__CLASSIFIEDONLY.m =========
% Session-wise cue-type decoding heatmap (L / C / R) across warped time.
%
% PURPOSE
%   - Build ONE global warped trial template from ALL sessions (behavior-only targets)
%   - For EACH SESSION:
%       * keep valid trials with cue / press / lick / cueType
%       * use ONLY classified units (cell_type in {1,2,3})
%       * build SDF (Hz) for every trial and every unit
%       * warp cue->press and press->lick to GLOBAL targets
%       * chop warped trial into bins
%       * for each time bin, form a trial x unit feature matrix
%       * split trials STRATIFIED within cue type: 75% train, 25% test
%       * train a within-session decoder on train trials and predict cueType on held-out trials
%       * store held-out decoding ACCURACY for every time bin
%
% OUTPUT
%   - One heatmap:
%       y-axis = session
%       x-axis = warped time bins
%       color  = held-out cue-type decoding accuracy in that session/time bin
%   - One Early vs Late trace with cluster-based permutation test across time:
%       * binwise two-sample t-statistic
%       * cluster-forming threshold at uncorrected p<0.05
%       * cluster mass = sum(abs(t)) within cluster
%       * permutation null = max cluster mass across shuffled group labels
%       * significant clusters marked above the trace
%   - Barplots:
%       * overall mean decoding accuracy (Early vs Late)
%       * cue-segment mean decoding accuracy (Early vs Late)
%       * press-segment mean decoding accuracy (Early vs Late)
%       * lick-segment mean decoding accuracy (Early vs Late)
%
% STYLE
%   - Matches your heatmap / decoding scripts:
%       * thick axes / large fonts
%       * cue / press / lick event lines
%       * plasma colormap from your decoding script
%       * session rows, warped time columns
%
% NOTES
%   - Global warp targets follow the "second code" logic:
%       behavior-only, across all sessions together
%   - Session trial inclusion for decoding uses:
%       valid == true, RT >= MIN_RT_MS, cueType in {L,C,R}, press/lick inside fixedWin
%   - Unit inclusion uses classified-only and active-trial criterion within fixedWin
%   - Decoder is WITHIN SESSION (train 75% of trials, test 25% of trials)
%   - Decoding is done SEPARATELY for each time bin
%
% FIXED
%   - Labels are handled as CATEGORICAL throughout decoding
%   - Removed string/isundefined mismatch that caused warnings
%
% NEW CHANGE
%   - Print overall MIN and MAX decoding accuracy across all sessions/bins in command window
%   - Limit displayed color scale to [0 0.6]
%   - Add cluster-based permutation test across time for Early vs Late trace
%   - Add Early vs Late barplots for overall and per-event/segment accuracy
%   - Save decoding results to .mat for downstream analyses

clear; clc;

%% ---- USER SETTINGS ----
matFile      = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL';
outDir  = fullfile(baseOut, 'DECODING_CUETYPE_HEATMAP');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Warped trial settings
preCueMs   = 500;
postLickMs = 3000;
fixedWin   = [-1000, 6000];   % used for session-trial inclusion + active-trial criterion

% Filters
MIN_RT_MS        = 100;
requireValid     = true;
minTrialsPerCue  = 2;    % minimum TOTAL trials per cue type in a session before split
minTrialsPerUnit = 2;    % active trials required to include a unit
minUnitsPerSess  = 2;    % minimum classified active units to decode a session

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Decoding settings
binSizeMs   = 200;       % requested
trainFrac   = 0.75;      % requested
doZscore    = true;
rngSeed     = 0;

% Group settings for trace / cluster test
earlySess = 1:8;
lateSess  = 29:36;

% Cluster permutation settings
Nperm_cluster = 10000;
clusterAlpha  = 0.05;

% Permutation tests for ALL barplots
Nperm_bar = 10000;

% Display
xUnit = "s";
forceIntegerSecondsTicks = true;

% Figure style (match your scripts)
figPos = [120, 120, 650, 900];
figPosTrace = [120, 120, 650, 300];

cbFontSize    = 22;
titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

eventLineLW = 5.0;
axesTickLW  = 4.0;
axesTickLen = [0.02 0.02];
traceLW     = 4;

scatterMS     = 8;
scatterEdgeLW = 1.0;
jit           = 0.10;
greyDot       = [0.6 0.6 0.6];

colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

colEarlyLine = [0.85 0.75 0.00]; % darker yellow
colLateLine  = [0.00 0.70 0.00]; % green

colEarlyShade = [1 1 0]; % yellow
colLateShade  = [0 1 0]; % green

% Accuracy color limits
climFixed = [0 0.8];

rng(rngSeed);

%% ---- BARPLOT OUTPUTS ----
outBarMeanPng  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__MEANACC.png');
outBarMeanSvg  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__MEANACC.svg');
outBarMeanFig  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__MEANACC.fig');

outCueBarPng   = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__CueSegment.png');
outCueBarSvg   = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__CueSegment.svg');
outCueBarFig   = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__CueSegment.fig');

outPressBarPng = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__PressSegment.png');
outPressBarSvg = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__PressSegment.svg');
outPressBarFig = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__PressSegment.fig');

outLickBarPng  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__LickSegment.png');
outLickBarSvg  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__LickSegment.svg');
outLickBarFig  = fullfile(outDir, 'CUETYPE_DECODING__EARLY_LATE_BAR__LickSegment.fig');

%% ---- LOAD BEHSTRUCT ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = numel(beh);

% Keep same practical session range convention as your recent scripts
nSessions = min(nSessions, 36);
beh = beh(1:nSessions);

earlySess = earlySess(earlySess>=1 & earlySess<=nSessions);
lateSess  = lateSess(lateSess>=1 & lateSess<=nSessions);

fprintf('Loaded beh with %d sessions.\n', nSessions);

%% ---- LOAD CELL TYPES ----
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), ...
    'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (behavior-only; across all sessions together) ----
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
        if lr <= rt, continue; end

        all_pressLat_global(end+1,1) = rt; %#ok<AGROW>
        all_lickLat_global(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat_global) && ~isempty(all_lickLat_global), ...
    'No trials passed global filters for warp targets.');

Tpress = median(all_pressLat_global);
Tlick  = median(all_lickLat_global);

fprintf('\n=== GLOBAL WARP TARGETS ===\n');
fprintf('Tpress = %.1f ms\n', Tpress);
fprintf('Tlick  = %.1f ms\n', Tlick);

%% ---- TEMPLATE AXIS ----
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

% Bin warped template into bins
binBins = max(1, round(binSizeMs / dt_ms));
binStarts = 1:binBins:nT;
nBins = numel(binStarts);
binCenters_idx = zeros(1,nBins);

for b = 1:nBins
    i0 = binStarts(b);
    i1 = min(nT, i0 + binBins - 1);
    binCenters_idx(b) = round(mean([i0 i1]));
end

binCenters_ms = tgrid_ms(binCenters_idx);

% labels by warped-time segment (same technique as comparison code)
segLabels = ones(nBins,1);  % 1 = cue->press
segLabels(binCenters_ms >= Tpress & binCenters_ms < Tlick) = 2; % 2 = press->lick
segLabels(binCenters_ms >= Tlick) = 3;                         % 3 = post-lick

%% ---- AXIS FOR PLOT ----
if xUnit == "s"
    xPlot = binCenters_ms / 1000;
    xlab  = 'Warped time from Cue (s)';
    cueLine   = 0;
    pressLine = Tpress / 1000;
    lickLine  = Tlick / 1000;
else
    xPlot = binCenters_ms;
    xlab  = 'Warped time from Cue (ms)';
    cueLine   = 0;
    pressLine = Tpress;
    lickLine  = Tlick;
end

%% ---- BUILD SESSION x TIMEBIN DECODING MATRIX ----
AccMat = nan(nBins, nSessions);
nUnitsUsed  = nan(nSessions,1);
nTrainUsed  = nan(nSessions,1);
nTestUsed   = nan(nSessions,1);

fprintf('\n=== SESSION-WISE CUE-TYPE DECODING ===\n');

for sIdx = 1:nSessions
    fprintf('\nSession %02d\n', sIdx);

    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session)
        fprintf('  skipped: bad/missing trial struct.\n');
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        fprintf('  skipped: no spikes.\n');
        continue;
    end

    % ---- session trial inclusion ----
    keepTrial = false(numel(trials),1);
    cueAbs    = nan(numel(trials),1);
    pressLat  = nan(numel(trials),1);
    lickLat   = nan(numel(trials),1);
    cueType   = strings(numel(trials),1);

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end

        if ~all(isfield(tr, {'cue','press','lick','cueType'})), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);

        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end

        pr = press - cue;
        lr = lick  - cue;
        if lr <= pr, continue; end

        % keep same fixedWin-style trial bounds as your heatmap script
        if ~(pr >= fixedWin(1) && pr <= fixedWin(2)), continue; end
        if ~(lr >= fixedWin(1) && lr <= fixedWin(2)), continue; end

        ct = upper(strtrim(string(tr.cueType)));
        if ~ismember(ct, ["L","C","R"]), continue; end

        keepTrial(k) = true;
        cueAbs(k)    = cue;
        pressLat(k)  = pr;
        lickLat(k)   = lr;
        cueType(k)   = ct;
    end

    idxKeep = find(keepTrial);
    if isempty(idxKeep)
        fprintf('  skipped: no kept trials.\n');
        continue;
    end

    cueAbsK   = cueAbs(idxKeep);
    pressLatK = pressLat(idxKeep);
    lickLatK  = lickLat(idxKeep);
    cueTypeK  = cueType(idxKeep);

    idxL_all = find(cueTypeK == "L");
    idxC_all = find(cueTypeK == "C");
    idxR_all = find(cueTypeK == "R");

    nL = numel(idxL_all); nC = numel(idxC_all); nR = numel(idxR_all);
    fprintf('  kept trials: total=%d | L=%d C=%d R=%d\n', numel(idxKeep), nL, nC, nR);

    if any([nL nC nR] < minTrialsPerCue)
        fprintf('  skipped: too few trials in at least one cue type.\n');
        continue;
    end

    % ---- stratified 75/25 split within cue type ----
    [trainLocal, testLocal] = stratifiedSplitLCR_(cueTypeK, trainFrac);

    yTrain = categorical(cueTypeK(trainLocal), ["L","C","R"]);
    yTest  = categorical(cueTypeK(testLocal),  ["L","C","R"]);

    if numel(categories(removecats(yTrain))) < 3 || numel(categories(removecats(yTest))) < 3
        fprintf('  skipped: split did not preserve all 3 cue classes in train/test.\n');
        continue;
    end

    nTrainUsed(sIdx) = numel(trainLocal);
    nTestUsed(sIdx)  = numel(testLocal);

    fprintf('  split: train=%d | test=%d\n', numel(trainLocal), numel(testLocal));

    % ---- collect per-unit warped binned activity for ALL kept trials ----
    unitTensor = [];  % [nTrials x nBins x nUnits]
    nUnits = 0;

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)
            ct = safe_get_celltype_(cell_type, sIdx, ch, u);
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                continue;
            end

            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
            spk_abs = double(spk_abs(:));

            % active-trial criterion
            activeTrials = false(numel(idxKeep),1);
            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                if any(spk_abs >= cue0 + fixedWin(1) & spk_abs <= cue0 + fixedWin(2))
                    activeTrials(iTr) = true;
                end
            end

            if nnz(activeTrials) < minTrialsPerUnit
                continue;
            end

            Xunit = nan(numel(idxKeep), nBins);

            for iTr = 1:numel(idxKeep)
                cue0   = cueAbsK(iTr);
                tPress = pressLatK(iTr);
                tLick  = lickLatK(iTr);

                winRight_trial = tLick + postLickMs;
                nT_trial = numel(winLeft:dt_ms:winRight_trial);

                idx0 = timeToIdx_(0,      winLeft, dt_ms, nT_trial);
                idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

                if ~(idx0 < idxP && idxP < idxL)
                    continue;
                end

                % ---- SDF (Hz) ----
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
                y  = sm / (dt_ms/1000);

                % ---- warp cue->press and press->lick to global template ----
                pre = y(1:idx0);

                segA = y(idx0:idxP);
                wA   = warpSegment_floorceilEqn_(segA, max(2, (idxPT - idx0T + 1)));

                segB = y(idxP:idxL);
                wB   = warpSegment_floorceilEqn_(segB, max(2, (idxLT - idxPT + 1)));

                segC = y(idxL:end);
                kC_target = max(1, (nT - idxLT + 1));
                if numel(segC) < kC_target
                    segC_fit = [segC zeros(1, kC_target - numel(segC))];
                else
                    segC_fit = segC(1:kC_target);
                end

                yy = [pre, wA(2:end), wB(2:end), segC_fit(2:end)];

                if numel(yy) < nT
                    yy = [yy zeros(1, nT-numel(yy))];
                elseif numel(yy) > nT
                    yy = yy(1:nT);
                end

                % ---- bin warped trial into bins ----
                xbin = nan(1, nBins);
                for b = 1:nBins
                    i0 = binStarts(b);
                    i1 = min(nT, i0 + binBins - 1);
                    xbin(b) = mean(yy(i0:i1), 'omitnan');
                end

                Xunit(iTr,:) = xbin;
            end

            % require all split trials to be finite for this unit
            if any(~all(isfinite(Xunit(trainLocal,:)),2)) || any(~all(isfinite(Xunit(testLocal,:)),2))
                continue;
            end

            nUnits = nUnits + 1;
            unitTensor(:,:,nUnits) = Xunit; %#ok<AGROW>
        end
    end

    if nUnits < minUnitsPerSess
        fprintf('  skipped: too few included units (%d).\n', nUnits);
        continue;
    end

    nUnitsUsed(sIdx) = nUnits;
    fprintf('  included units: %d\n', nUnits);

    % ---- decode cue type separately for each time bin ----
    accBin = nan(1, nBins);

    for b = 1:nBins
        Xall = squeeze(unitTensor(:, b, :));  % [nTrials x nUnits]
        if isvector(Xall)
            Xall = Xall(:);
        end

        Xtr = Xall(trainLocal, :);
        Xte = Xall(testLocal,  :);

        okTr = all(isfinite(Xtr),2);
        okTe = all(isfinite(Xte),2);

        Xtr = Xtr(okTr,:);
        Xte = Xte(okTe,:);

        yTr = yTrain(okTr);
        yTe = yTest(okTe);

        if size(Xtr,1) < 6 || size(Xte,1) < 3
            continue;
        end
        if numel(categories(removecats(yTr))) < 3 || numel(categories(removecats(yTe))) < 3
            continue;
        end

        if doZscore
            [Xtr, mu, sig] = zscoreSafe_(Xtr);
            Xte = (Xte - mu) ./ sig;
            Xte(~isfinite(Xte)) = 0;
        end

        try
            t = templateLinear('Learner','logistic', 'Lambda',1e-4, 'Regularization','ridge');
            Mdl = fitcecoc(Xtr, yTr, 'Learners', t, 'ClassNames', categorical(["L","C","R"]));
            yHat = predict(Mdl, Xte);

            ok = ~isundefined(yHat) & ~isundefined(yTe);
            if nnz(ok) < 3
                continue;
            end

            accBin(b) = mean(yHat(ok) == yTe(ok));
        catch ME
            warning('Session %02d bin %d decode failed: %s', sIdx, b, ME.message);
            accBin(b) = nan;
        end
    end

    AccMat(:, sIdx) = accBin(:);
end

%% ---- PRINT OVERALL MIN / MAX ----
allVals = AccMat(isfinite(AccMat));
if isempty(allVals)
    fprintf('\n=== OVERALL DECODING RANGE ===\n');
    fprintf('No finite decoding values found across all sessions/bins.\n');
else
    fprintf('\n=== OVERALL DECODING RANGE ===\n');
    fprintf('Min decoding accuracy across all sessions/bins: %.6f\n', min(allVals));
    fprintf('Max decoding accuracy across all sessions/bins: %.6f\n', max(allVals));
end

%% ---- PLOT HEATMAP ----
Pimg = AccMat';

fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

hImg = imagesc(ax, xPlot, 1:nSessions, Pimg);
axis(ax,'tight');
set(ax,'YDir','reverse');
set(hImg,'Interpolation','nearest');

cmap = plasma(256);
colormap(ax, cmap);
caxis(ax, climFixed);

xlabel(ax, xlab, 'FontSize', labelFontSize);
ylabel(ax, 'Session', 'FontSize', labelFontSize);
title(ax, 'Cue-type decoding accuracy per session', 'FontWeight','bold', 'FontSize', titleFontSize);

xline(ax, cueLine,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
xline(ax, pressLine, '--', 'Color', colPress, 'LineWidth', eventLineLW);
xline(ax, lickLine,  '--', 'Color', colLick,  'LineWidth', eventLineLW);

majorY = 1:5:nSessions;
if majorY(end) ~= nSessions
    majorY = unique([majorY nSessions]);
end
yticks(ax, majorY);
yticklabels(ax, string(majorY));
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = setdiff(1:nSessions, majorY);

set(ax, 'FontSize', tickFontSize, ...
    'TickDir','out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box','off', ...
    'Layer','top');

if xUnit == "s" && forceIntegerSecondsTicks
    xl = xlim(ax);
    ts = ceil(xl(1));
    te = floor(xl(2));
    if te >= ts
        xticks(ax, ts:te);
    end
end

%% ---- COLORBAR ----
cb = colorbar(ax, 'eastoutside');
cb.Label.String = 'Accuracy';
cb.FontSize = cbFontSize;
cb.Label.FontSize = cbFontSize;
cb.LineWidth = 2;
cb.TickLength = 0.02;
cb.Ticks = 0:0.2:0.8;

axpos = ax.Position;
cbpos = cb.Position;
cbpos(3) = 0.040;
cbpos(4) = 0.50 * axpos(4);
cbpos(2) = axpos(2) + (axpos(4)-cbpos(4))/2;
cbpos(1) = axpos(1) + axpos(3) + 0.07*axpos(3);
cbpos(1) = min(cbpos(1), 0.92);
cb.Position = cbpos;

%% ---- SAVE HEATMAP ----
outbase = fullfile(outDir, sprintf('Sessions01_%02d_CueTypeDecodingHeatmap_GLOBALWARP_BIN%03dms_train%02d', ...
    nSessions, binSizeMs, round(100*trainFrac)));

set(fig, 'PaperPositionMode', 'auto');

exportgraphics(fig, [outbase '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
savefig(fig, [outbase '.fig']);
try
    saveas(fig, [outbase '.svg']);
catch
    print(fig, [outbase '.svg'], '-dsvg');
end

fprintf('\nSaved:\n  %s.png\n  %s.fig\n  %s.svg\n', outbase, outbase, outbase);

%% ---- EARLY vs LATE TRACE + CLUSTER-BASED PERMUTATION TEST ----
earlyMat = AccMat(:, earlySess)';   % nEarly x nBins
lateMat  = AccMat(:, lateSess)';    % nLate  x nBins

earlyMat = earlyMat(any(isfinite(earlyMat),2), :);
lateMat  = lateMat(any(isfinite(lateMat),2),  :);

mEarly = mean(earlyMat, 1, 'omitnan');
mLate  = mean(lateMat,  1, 'omitnan');

seEarly = std(earlyMat, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(earlyMat),1)));
seLate  = std(lateMat,  0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(lateMat),1)));

obsT = nan(1, nBins);
obsP = nan(1, nBins);

for b = 1:nBins
    e = earlyMat(:,b); e = e(isfinite(e));
    l = lateMat(:,b);  l = l(isfinite(l));

    if numel(e) >= 2 && numel(l) >= 2
        [~, p, ~, stats] = ttest2(e, l, 'Vartype','unequal');
        obsT(b) = stats.tstat;
        obsP(b) = p;
    end
end

clusterMaskObs = isfinite(obsP) & (obsP < clusterAlpha);
obsClusters = findClusters1D_(clusterMaskObs);

obsClusterMass = nan(numel(obsClusters),1);
for iC = 1:numel(obsClusters)
    idx = obsClusters{iC};
    obsClusterMass(iC) = nansum(abs(obsT(idx)));
end

allData = [earlyMat; lateMat];
nEarly = size(earlyMat,1);
nAll   = size(allData,1);

maxClusterMassNull = zeros(Nperm_cluster,1);

rng(rngSeed);
for iPerm = 1:Nperm_cluster
    permIdx = randperm(nAll);
    idxE = permIdx(1:nEarly);
    idxL = permIdx(nEarly+1:end);

    permT = nan(1, nBins);
    permP = nan(1, nBins);

    for b = 1:nBins
        e = allData(idxE,b); e = e(isfinite(e));
        l = allData(idxL,b); l = l(isfinite(l));

        if numel(e) >= 2 && numel(l) >= 2
            [~, p, ~, stats] = ttest2(e, l, 'Vartype','unequal');
            permT(b) = stats.tstat;
            permP(b) = p;
        end
    end

    permMask = isfinite(permP) & (permP < clusterAlpha);
    permClusters = findClusters1D_(permMask);

    if isempty(permClusters)
        maxClusterMassNull(iPerm) = 0;
    else
        permMasses = zeros(numel(permClusters),1);
        for jC = 1:numel(permClusters)
            idx = permClusters{jC};
            permMasses(jC) = nansum(abs(permT(idx)));
        end
        maxClusterMassNull(iPerm) = max(permMasses);
    end
end

obsClusterP = nan(numel(obsClusters),1);
for iC = 1:numel(obsClusters)
    obsClusterP(iC) = (1 + nnz(maxClusterMassNull >= obsClusterMass(iC))) / (Nperm_cluster + 1);
end

fprintf('\n=== CLUSTER-BASED PERMUTATION TEST (EARLY vs LATE) ===\n');
fprintf('nEarly=%d, nLate=%d, clusterAlpha=%.3f, Nperm=%d\n', size(earlyMat,1), size(lateMat,1), clusterAlpha, Nperm_cluster);
if isempty(obsClusters)
    fprintf('No observed supra-threshold clusters.\n');
else
    for iC = 1:numel(obsClusters)
        idx = obsClusters{iC};
        fprintf('Cluster %d: bins %d-%d | time %.3f to %.3f s | mass=%.6f | p=%.6g\n', ...
            iC, idx(1), idx(end), xPlot(idx(1)), xPlot(idx(end)), obsClusterMass(iC), obsClusterP(iC));
    end
end

figTrace = figure('Color','w','Position',figPosTrace);
ax2 = axes(figTrace); hold(ax2,'on');

patch(ax2, [xPlot, fliplr(xPlot)], [mEarly-seEarly, fliplr(mEarly+seEarly)], ...
    colEarlyLine, 'FaceAlpha', 0.15, 'EdgeColor','none');
patch(ax2, [xPlot, fliplr(xPlot)], [mLate-seLate, fliplr(mLate+seLate)], ...
    colLateLine, 'FaceAlpha', 0.15, 'EdgeColor','none');

plot(ax2, xPlot, mEarly, '-', 'Color', colEarlyLine, 'LineWidth', traceLW);
plot(ax2, xPlot, mLate,  '-', 'Color', colLateLine,  'LineWidth', traceLW);

xline(ax2, cueLine,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
xline(ax2, pressLine, '--', 'Color', colPress, 'LineWidth', eventLineLW);
xline(ax2, lickLine,  '--', 'Color', colLick,  'LineWidth', eventLineLW);

xlabel(ax2, xlab, 'FontSize', labelFontSize);
ylabel(ax2, 'Accuracy', 'FontSize', labelFontSize);
title(ax2, 'Cue-type decoding accuracy: Early vs Late', 'FontWeight','bold', 'FontSize', titleFontSize);

legend(ax2, {'Early','Late'}, 'Location','northeast', 'FontSize', tickFontSize);

set(ax2, 'FontSize', tickFontSize, ...
    'TickDir','out', ...
    'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, ...
    'Box','off', ...
    'Layer','top');

if xUnit == "s" && forceIntegerSecondsTicks
    xl = xlim(ax2);
    ts = ceil(xl(1));
    te = floor(xl(2));
    if te >= ts
        xticks(ax2, ts:te);
    end
end

% draw significant clusters above the traces
yl = ylim(ax2);
yRange = yl(2) - yl(1);
yBase = yl(2) + 0.05*yRange;
yTop  = yl(2) + 0.12*yRange;

sigCount = 0;
for iC = 1:numel(obsClusters)
    if obsClusterP(iC) < 0.05
        idx = obsClusters{iC};
        x1 = xPlot(idx(1));
        x2 = xPlot(idx(end));
        sigCount = sigCount + 1;

        plot(ax2, [x1 x2], [yBase yBase], 'k-', 'LineWidth', axesTickLW, 'Clipping','off');
        text(ax2, mean([x1 x2]), yTop, '*', ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', titleFontSize, 'FontWeight','bold', 'Color', [0 0 0], 'Clipping','off');
    end
end

if sigCount > 0
    ylim(ax2, [yl(1), yl(2) + 0.18*yRange]);
end

outbaseTrace = fullfile(outDir, sprintf('Sessions01_%02d_CueTypeDecodingTrace_CLUSTERPERM_GLOBALWARP_BIN%03dms_train%02d', ...
    nSessions, binSizeMs, round(100*trainFrac)));

set(figTrace, 'PaperPositionMode', 'auto');
exportgraphics(figTrace, [outbaseTrace '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
savefig(figTrace, [outbaseTrace '.fig']);
try
    saveas(figTrace, [outbaseTrace '.svg']);
catch
    print(figTrace, [outbaseTrace '.svg'], '-dsvg');
end

fprintf('Saved:\n  %s.png\n  %s.fig\n  %s.svg\n', outbaseTrace, outbaseTrace, outbaseTrace);

%% ---- BARPLOTS OF ACCURACY AND ACCURACY PER EVENT ----
% Session-level overall mean accuracy
sessMeanAcc = mean(AccMat, 1, 'omitnan');

% Session-level segment means using same segment-definition technique as comparison code
cueMask   = (segLabels == 1);
pressMask = (segLabels == 2);
lickMask  = (segLabels == 3);

sessCueAcc   = nan(1, nSessions);
sessPressAcc = nan(1, nSessions);
sessLickAcc  = nan(1, nSessions);

for sIdx = 1:nSessions
    sessCueAcc(sIdx)   = mean(AccMat(cueMask,   sIdx), 'omitnan');
    sessPressAcc(sIdx) = mean(AccMat(pressMask, sIdx), 'omitnan');
    sessLickAcc(sIdx)  = mean(AccMat(lickMask,  sIdx), 'omitnan');
end

earlyMeanAcc = sessMeanAcc(earlySess);
lateMeanAcc  = sessMeanAcc(lateSess);

earlyCueAcc   = sessCueAcc(earlySess);
lateCueAcc    = sessCueAcc(lateSess);

earlyPressAcc = sessPressAcc(earlySess);
latePressAcc  = sessPressAcc(lateSess);

earlyLickAcc  = sessLickAcc(earlySess);
lateLickAcc   = sessLickAcc(lateSess);

% Overall mean accuracy bar
figBarMean = plotGroupBar_perm_(earlyMeanAcc, lateMeanAcc, colEarlyShade, colLateShade, ...
    'Accuracy', 'Mean accuracy (cue-type decoding)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figBarMean, outBarMeanPng); savefig(figBarMean, outBarMeanFig);
try, print(figBarMean, outBarMeanSvg, '-dsvg'); catch ME, warning('SVG save failed (mean accuracy): %s', ME.message); end
close(figBarMean);

% Cue-segment bar
figCueBar = plotGroupBar_perm_(earlyCueAcc, lateCueAcc, colEarlyShade, colLateShade, ...
    'Accuracy', 'Cue-segment accuracy (cue-type decoding)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figCueBar, outCueBarPng); savefig(figCueBar, outCueBarFig);
try, print(figCueBar, outCueBarSvg, '-dsvg'); catch ME, warning('SVG save failed (cue segment): %s', ME.message); end
close(figCueBar);

% Press-segment bar
figPressBar = plotGroupBar_perm_(earlyPressAcc, latePressAcc, colEarlyShade, colLateShade, ...
    'Accuracy', 'Press-segment accuracy (cue-type decoding)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figPressBar, outPressBarPng); savefig(figPressBar, outPressBarFig);
try, print(figPressBar, outPressBarSvg, '-dsvg'); catch ME, warning('SVG save failed (press segment): %s', ME.message); end
close(figPressBar);

% Lick-segment bar
figLickBar = plotGroupBar_perm_(earlyLickAcc, lateLickAcc, colEarlyShade, colLateShade, ...
    'Accuracy', 'Lick-segment accuracy (cue-type decoding)', ...
    labelFontSize, titleFontSize, tickFontSize, axesTickLW, ...
    scatterMS, scatterEdgeLW, jit, greyDot, Nperm_bar);
saveas(figLickBar, outLickBarPng); savefig(figLickBar, outLickBarFig);
try, print(figLickBar, outLickBarSvg, '-dsvg'); catch ME, warning('SVG save failed (lick segment): %s', ME.message); end
close(figLickBar);

%% ---- SAVE DECODING RESULTS MAT ----
outMat = fullfile(outDir, sprintf('Sessions01_%02d_CueTypeDecodingResults_GLOBALWARP_BIN%03dms_train%02d.mat', ...
    nSessions, binSizeMs, round(100*trainFrac)));

save(outMat, ...
    'AccMat', ...
    'binCenters_ms', ...
    'xPlot', ...
    'Tpress', ...
    'Tlick', ...
    'earlySess', ...
    'lateSess', ...
    'nSessions', ...
    'nBins', ...
    'binSizeMs', ...
    'trainFrac', ...
    'segLabels', ...
    'cueLine', ...
    'pressLine', ...
    'lickLine', ...
    'nUnitsUsed', ...
    'nTrainUsed', ...
    'nTestUsed', ...
    'sessMeanAcc', ...
    'sessCueAcc', ...
    'sessPressAcc', ...
    'sessLickAcc', ...
    '-v7.3');

fprintf('Saved decoding results MAT:\n  %s\n', outMat);

%% ---- OPTIONAL SUMMARY PRINT ----
fprintf('\n=== SUMMARY ===\n');
for sIdx = 1:nSessions
    nGoodBins = nnz(isfinite(AccMat(:,sIdx)));
    fprintf('Session %02d: units=%g train=%g test=%g decoded_bins=%d/%d\n', ...
        sIdx, nUnitsUsed(sIdx), nTrainUsed(sIdx), nTestUsed(sIdx), nGoodBins, nBins);
end

%% ================= HELPERS =================

function [Xz, mu, sig] = zscoreSafe_(X)
    mu = mean(X, 1, 'omitnan');
    sig = std(X, 0, 1, 'omitnan');
    sig(sig <= 0 | ~isfinite(sig)) = 1;
    Xz = (X - mu) ./ sig;
    Xz(~isfinite(Xz)) = 0;
end

function [trainIdx, testIdx] = stratifiedSplitLCR_(cueTypeK, trainFrac)
    trainIdx = [];
    testIdx  = [];

    cats = ["L","C","R"];
    for i = 1:numel(cats)
        idx = find(cueTypeK == cats(i));
        idx = idx(randperm(numel(idx)));

        nTrain = max(1, floor(trainFrac * numel(idx)));
        nTrain = min(nTrain, numel(idx)-1);  % preserve at least one test if possible

        tr = idx(1:nTrain);
        te = idx(nTrain+1:end);

        if isempty(te)
            te = tr(end);
            tr = tr(1:end-1);
        end

        trainIdx = [trainIdx; tr(:)]; %#ok<AGROW>
        testIdx  = [testIdx;  te(:)]; %#ok<AGROW>
    end

    trainIdx = sort(trainIdx);
    testIdx  = sort(testIdx);
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k = 1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k}))
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

function ct = safe_get_celltype_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
       ch   <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
       u    <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
    halfWidth = ceil(5 * sigmaMs / dtMs);
    x = (-halfWidth:halfWidth) * dtMs;
    g = exp(-0.5 * (x ./ sigmaMs).^2);
    g = g / sum(g);
end

function idx = timeToIdx_(t_ms, winLeft_ms, dt_ms, nT)
    idx = round((t_ms - winLeft_ms) / dt_ms) + 1;
    idx = max(1, min(nT, idx));
end

function ywarp = warpSegment_floorceilEqn_(y, k)
    y = y(:)';
    n = numel(y);

    if n < 2 || k < 2
        ywarp = zeros(1, k);
        if n >= 1 && k >= 1
            ywarp(1) = y(1);
        end
        return;
    end

    s = (k - 1) / (n - 1);
    ywarp = zeros(1, k);

    for i = 1:k
        x = 1 + (i - 1) / s;
        x0 = floor(x);
        x1 = ceil(x);

        x0 = max(1, min(n, x0));
        x1 = max(1, min(n, x1));

        if x0 == x1
            ywarp(i) = y(x0);
        else
            w1 = x - x0;
            w0 = 1 - w1;
            ywarp(i) = w0 * y(x0) + w1 * y(x1);
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

function clusters = findClusters1D_(mask)
    clusters = {};
    if isempty(mask) || ~any(mask), return; end

    d = diff([false, mask(:)', false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;

    for i = 1:numel(starts)
        clusters{i} = starts(i):stops(i); %#ok<AGROW>
    end
end

function fig = plotGroupBar_perm_(earlyVals, lateVals, colEarly, colLate, ylab, ttl, ...
    labelFS, titleFS, tickFS, axLW, ms, edgeLW, jit, greyDot, Nperm)

    earlyVals = earlyVals(isfinite(earlyVals));
    lateVals  = lateVals(isfinite(lateVals));

    fig = figure('Color','w','Position',[120,120,520,700]);
    ax = axes(fig); hold(ax,'on');

    mE = mean(earlyVals,'omitnan');
    mL = mean(lateVals,'omitnan');
    semE = std(earlyVals,'omitnan') / sqrt(max(1,numel(earlyVals)));
    semL = std(lateVals,'omitnan')  / sqrt(max(1,numel(lateVals)));

    xPos = [1 2];
    barW = 0.60;

    b = bar(ax, xPos, [mE mL], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axLW);
    b.CData(1,:) = colEarly; b.CData(2,:) = colLate;
    b.FaceAlpha = 0.25;

    errorbar(ax, xPos, [mE mL], [semE semL], 'k', 'LineStyle','none', 'LineWidth', axLW, 'CapSize', 18);

    plot(ax, 1 + (rand(size(earlyVals))-0.5)*2*jit, earlyVals, 'o', 'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', edgeLW);
    plot(ax, 2 + (rand(size(lateVals))-0.5)*2*jit,  lateVals,  'o', 'MarkerSize', ms, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', edgeLW);

    set(ax, 'XLim',[0.5 2.5], 'XTick',[1 2], 'XTickLabel',{'Early','Late'});
    ylabel(ax, ylab, 'FontSize', labelFS);
    title(ax, ttl, 'FontWeight','bold', 'FontSize', titleFS);
    set(ax, 'FontSize', tickFS, 'Box','off');
    set(ax, 'TickDir','out', 'LineWidth', axLW, 'Layer','top');

    % permutation test (label-shuffle) on Δ(Late-Early)
    if ~isempty(earlyVals) && ~isempty(lateVals)
        dObs = mean(lateVals,'omitnan') - mean(earlyVals,'omitnan');
        allv = [earlyVals(:); lateVals(:)];
        nE = numel(earlyVals);

        rng(0);
        dPerm = nan(Nperm,1);
        for i = 1:Nperm
            perm = allv(randperm(numel(allv)));
            e = perm(1:nE);
            l = perm(nE+1:end);
            dPerm(i) = mean(l,'omitnan') - mean(e,'omitnan');
        end
        p = (1 + nnz(abs(dPerm) >= abs(dObs))) / (Nperm + 1);

        fprintf('\n%s perm-test: Δ(Late-Early)=%.6g, p=%.6g (Nperm=%d)\n', ttl, dObs, p, Nperm);
        addSigIfNeeded_(ax, 1, 2, p, [mE mL], axLW, tickFS);
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