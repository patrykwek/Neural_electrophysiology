%% ========= FIG3_SPLITHALF_TEMPLATE_RELIABILITY__A_THESIS_FINAL__GLOBALWARP.m =========
% Within-session template reliability (even vs odd trials), using your SAME global-warped SDF pipeline.
%
% UPDATED (requested):
%   - Reliability is now computed PER UNIT first
%   - Keep units that fire in at least 2 trials in a session
%   - For each unit:
%       * compute cue-aligned SDF per trial (10 ms bins; sigma = 25 ms)
%       * globally warp each trial to shared press/lick targets
%       * split active trials into odd vs even
%       * summarize each half across trials using the MEDIAN SDF trace
%       * compute unit reliability as corr(medianOddTrace, medianEvenTrace)
%   - Session reliability is then summarized across units using the MEDIAN
%
% Outputs:
%   (1) Split-half reliability vs session (Weibull-axis tick style + stage shading + dashed boundary)
%   (2) Early vs Late bar (same aesthetics as your decoding-stage bars)
%   (3) Null / chance baseline panels:
%       (A) Histogram of shuffled split-half r (null) with observed median session r as vertical line
%       (B) Paired observed vs shuffled (per-session) with mean±SEM + paired permutation test
%   (4) Performance (% correct) vs split-half r scatter (straight line overlay)
%
% NOTE:
%   - Plot names/variables are kept the same where possible.
%   - Focuses ONLY on split-half reliability figures (no SIM/MDS/drift here).

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/3_2_SPLITHALF_RELIABILITY';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Windowing (same as your SIM script)
preCueMs   = 500;
postLickMs = 3000;
fixedWin   = [-1000, 6000];   % trial inclusion filters (relative to cue)

% Filters
MIN_RT_MS = 100;
minTrialsPerUnit = 2;   % UPDATED: requested minimum is 2 active trials
requireValid = true;

% SDF
dt_ms        = 10;
gaussSigmaMs = 25;

% Sessions
maxSessions = 36;

% Groups (match your thesis convention)
earlyGroup = 1:8;
lateGroup  = 29:36;
boundaryAt = 8.5;

% ---- NULL / SHUFFLE SETTINGS ----
nShuffle = 2000;              % per-session shuffles
rngSeedShuffle = 0;           % deterministic

% ---- PERMUTATION TEST SETTINGS ----
Nperm = 10000;

% ---- STYLE (match your Weibull/drift style) ----
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

% Stage shading + boundary
bgAlpha   = 0.10;
dashColor = [0.2 0.2 0.2];

% Plot color (reuse your drift blue)
colLine   = [0 0.45 0.75];

% Bar aesthetics (match decoding-stage bars)
colEarlyShade = [1 1 0]; % yellow
colLateShade  = [0 1 0]; % green
scatterMS     = 8;
scatterEdgeLW = 1.0;
jit           = 0.10;
greyDot       = [0.6 0.6 0.6];

rng(0);

%% ---- LOAD DATA ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = min(numel(beh), maxSessions);
beh = beh(1:nSessions);

assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u} : 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes

% ---- PRECOMPUTE SDF KERNEL ----
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
    'No trials passed global filters; cannot set shared warp targets.');

Tpress = median(all_pressLat_global);
Tlick  = median(all_lickLat_global);

fprintf('\n=== GLOBAL targets: Tpress=%.1f ms, Tlick=%.1f ms ===\n', Tpress, Tlick);

%% ---- TEMPLATE AXIS ----
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

%% ---- SPLIT-HALF RELIABILITY PER SESSION ----
splitR     = nan(nSessions,1);
nUnitsUsed = nan(nSessions,1);
nTrialsU   = nan(nSessions,1);

% store per-session null medians and pool null values
splitR_null_median = nan(nSessions,1);
splitR_null_all = [];  % pooled null r across sessions (for histogram)

% store per-session null mean as well
splitR_null_mean = nan(nSessions,1);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Behavior-only keep trials (same logic as SIM script) ----
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

    unitR = nan(0,1);
    unitNtr = nan(0,1);
    nullR_byUnit = nan(0, nShuffle);

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

            idxEven = 2:2:nTrU;
            idxOdd  = 1:2:nTrU;

            % need both halves
            if isempty(idxEven) || isempty(idxOdd)
                continue;
            end

            warpedHz = zeros(nTrU, nT);

            for iTr = 1:nTrU
                cue0   = cueU(iTr);
                tPress = rtU(iTr);
                tLick  = lkU(iTr);

                winRight_trial = tLick + postLickMs;
                nT_trial = numel(winLeft:dt_ms:winRight_trial);

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

                warpedHz(iTr,:) = yy;
            end

            % ---- Per-unit observed split-half using MEDIAN across trials ----
            muEven = median(warpedHz(idxEven,:), 1, 'omitnan');
            muOdd  = median(warpedHz(idxOdd,:),  1, 'omitnan');

            rObs = corr(muEven(:), muOdd(:), 'Rows','complete');
            if ~isfinite(rObs)
                continue;
            end

            unitR(end+1,1) = rObs; %#ok<AGROW>
            unitNtr(end+1,1) = nTrU; %#ok<AGROW>

            % ---- Per-unit NULL via circular time-shift of odd median trace ----
            rng(rngSeedShuffle + sIdx); %#ok<RNGR>
            thisNull = nan(1, nShuffle);
            muOddCol = muOdd(:);

            for b = 1:nShuffle
                sh = randi([1, nT-1], 1, 1); % random non-zero circular shift
                oddShift = circshift(muOddCol, sh);
                thisNull(b) = corr(muEven(:), oddShift, 'Rows','complete');
            end
            nullR_byUnit(end+1,:) = thisNull; %#ok<AGROW>
        end
    end

    if isempty(unitR)
        continue;
    end

    % ---- Session observed split-half: MEDIAN across units ----
    splitR(sIdx) = median(unitR, 'omitnan');
    nUnitsUsed(sIdx) = numel(unitR);
    nTrialsU(sIdx) = median(unitNtr, 'omitnan');

    % ---- Session null: for each shuffle, MEDIAN across units ----
    nullR = median(nullR_byUnit, 1, 'omitnan')';
    splitR_null_median(sIdx) = median(nullR, 'omitnan');
    splitR_null_mean(sIdx)   = mean(nullR, 'omitnan');
    splitR_null_all = [splitR_null_all; nullR(:)]; %#ok<AGROW>
end

fprintf('\nComputed split-half for %d/%d sessions.\n', nnz(isfinite(splitR)), nSessions);

%% ---- PRINT SUMMARY STATS FOR REPORTING ----
finiteAll = splitR(isfinite(splitR));
fprintf('\n===== Split-half template reliability summary (all sessions with finite r) =====\n');
fprintf('nSessions (finite r) = %d / %d\n', numel(finiteAll), nSessions);
if ~isempty(finiteAll)
    fprintf('Mean r = %.4f\n', mean(finiteAll, 'omitnan'));
    fprintf('Median r = %.4f\n', median(finiteAll, 'omitnan'));
    fprintf('Std r = %.4f\n', std(finiteAll, 0, 'omitnan'));
    fprintf('Min r = %.4f\n', min(finiteAll));
    fprintf('Max r = %.4f\n', max(finiteAll));
    fprintf('25th percentile = %.4f\n', prctile(finiteAll, 25));
    fprintf('75th percentile = %.4f\n', prctile(finiteAll, 75));
end

%% ---- Spearman vs session index ----
okS = isfinite(splitR) & isfinite((1:nSessions)');
if nnz(okS) >= 3
    [rhoSess, pSess] = corr((1:nSessions)', splitR, 'Type','Spearman', 'Rows','complete');
else
    rhoSess = nan; pSess = nan;
end
fprintf('\n===== Spearman: split-half r vs session index =====\n');
fprintf('n=%d sessions, rho=%.4f, p=%.6g\n', nnz(okS), rhoSess, pSess);

%% ---- Spearman vs performance (% correct) if Hit field exists ----
perfObs = nan(nSessions,1);
for sIdx = 1:nSessions
    sess = beh(sIdx);
    hit = [];
    if isfield(sess,'Hit') && ~isempty(sess.Hit)
        hit = sess.Hit;
    elseif isfield(sess,'trials') && isstruct(sess.trials) && isfield(sess.trials,'Hit')
        try, hit = [sess.trials.Hit]; catch, hit = []; end
    end
    if isempty(hit), continue; end
    hit = double(hit(:));
    hit = hit(isfinite(hit));
    hit = hit(hit==0 | hit==1);
    if isempty(hit), continue; end
    perfObs(sIdx) = 100 * mean(hit==1);
end

okP = isfinite(perfObs) & isfinite(splitR);
if nnz(okP) >= 3
    [rhoPerf, pPerf] = corr(perfObs, splitR, 'Type','Spearman', 'Rows','complete');
else
    rhoPerf = nan; pPerf = nan;
end
fprintf('\n===== Spearman: split-half r vs performance (%% correct) =====\n');
fprintf('n=%d sessions, rho=%.4f, p=%.6g\n', nnz(okP), rhoPerf, pPerf);

%% =========================
% PLOT 0: Null distribution histogram + observed line
%% =========================
figNull = figure('Color','w','Position',[120,120,720,520]);
axN = axes(figNull); hold(axN,'on');

nullVals = splitR_null_all(isfinite(splitR_null_all));
obsVals  = splitR(isfinite(splitR));

if isempty(nullVals)
    nullVals = 0;
end

edges = linspace(-1, 1, 40);
histogram(axN, nullVals, edges, 'Normalization','probability', ...
    'FaceColor',[0.7 0.7 0.7], 'EdgeColor',[0.7 0.7 0.7], 'FaceAlpha',0.55);

xObs = mean(obsVals, 'omitnan');
xline(axN, xObs, '-', 'LineWidth', 5, 'Color', colLine);

% permutation p-value: compare observed mean to null mean distribution via resampling from pooled null
nPermHist = 5000;
rng(0);
Tobs = mean(obsVals, 'omitnan');
Tnull = nan(nPermHist,1);
for p = 1:nPermHist
    rr = nullVals(randi(numel(nullVals), [numel(obsVals) 1]));
    Tnull(p) = mean(rr, 'omitnan');
end
pNull = (1 + nnz(Tnull >= Tobs)) / (nPermHist + 1);

fprintf('\n===== Null (pooled) vs observed mean split-half r =====\n');
fprintf('Observed mean r=%.4f, permutation p=%.6g (N=%d)\n', Tobs, pNull, nPermHist);

xlabel(axN, 'Split-half r', 'FontSize', labelFontSize);
ylabel(axN, 'Probability', 'FontSize', labelFontSize);
title(axN, sprintf('Shuffle null (pooled) vs observed mean (p=%.3g)', pNull), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axN, 'FontSize', tickFontSize, 'Box','off');
set(axN, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');
xlim(axN, [-1 1]);

outBaseNull = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_NULL_HIST_%02dSessions', nSessions));
exportgraphics(figNull, [outBaseNull '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figNull, [outBaseNull '.svg']); catch, print(figNull, [outBaseNull '.svg'], '-dsvg'); end
savefig(figNull, [outBaseNull '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseNull);

%% =========================
% PLOT 0b: Paired observed vs shuffled (per-session) + paired permutation test
%% =========================
figPair = figure('Color','w','Position',[120,120,560,620]);
axP = axes(figPair); hold(axP,'on');

okPair = isfinite(splitR) & isfinite(splitR_null_median);
obsP = splitR(okPair);
nulP = splitR_null_median(okPair);

x1 = 1; x2 = 2;

% paired lines
for i = 1:numel(obsP)
    plot(axP, [x1 x2], [obsP(i) nulP(i)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
end

% points
plot(axP, x1 + (rand(size(obsP))-0.5)*2*0.05, obsP, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', colLine, 'MarkerEdgeColor', colLine);
plot(axP, x2 + (rand(size(nulP))-0.5)*2*0.05, nulP, 'o', ...
    'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.5 0.5 0.5]);

% mean±SEM
m1 = mean(obsP,'omitnan'); s1 = std(obsP,0,'omitnan')/sqrt(max(1,numel(obsP)));
m2 = mean(nulP,'omitnan'); s2 = std(nulP,0,'omitnan')/sqrt(max(1,numel(nulP)));
errorbar(axP, [x1 x2], [m1 m2], [s1 s2], 'k', 'LineStyle','none', 'LineWidth', 3, 'CapSize', 16);

% paired permutation test (sign-flip on differences)
diffs = obsP(:) - nulP(:);
Tobs_pair = mean(diffs, 'omitnan');

Tperm = nan(Nperm,1);
for ip = 1:Nperm
    flip = (rand(size(diffs)) > 0.5) * 2 - 1; % +/-1
    Tperm(ip) = mean(diffs .* flip, 'omitnan');
end
p_pair = (1 + sum(abs(Tperm) >= abs(Tobs_pair))) / (1 + Nperm);

fprintf('\n===== Observed vs null (paired permutation) =====\n');
fprintf('n=%d sessions, mean(diff)=%.4f, p=%.6g (N=%d)\n', numel(diffs), Tobs_pair, p_pair, Nperm);

set(axP, 'XLim',[0.5 2.5], 'XTick',[1 2], 'XTickLabel',{'Observed','Shuffled'});
ylabel(axP, 'Split-half r', 'FontSize', labelFontSize);
title(axP, sprintf('Observed vs shuffle (paired p=%.3g)', p_pair), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axP, 'FontSize', tickFontSize, 'Box','off');
set(axP, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');
ylim(axP, [-1 1]);

outBasePair = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_OBS_VS_SHUFFLED_%02dSessions', nSessions));
exportgraphics(figPair, [outBasePair '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figPair, [outBasePair '.svg']); catch, print(figPair, [outBasePair '.svg'], '-dsvg'); end
savefig(figPair, [outBasePair '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBasePair);

%% =========================
% PLOT 1: Split-half vs session (Weibull tick style) + null line
%% =========================
figRel = figure('Color','w','Position',figPosDrift);
ax = axes(figRel); hold(ax,'on');

finiteR = splitR(isfinite(splitR));
if isempty(finiteR)
    yl = [-1 1];
else
    pad = 0.05;
    yl = [max(-1, min(finiteR)-pad), min(1, max(finiteR)+pad)];
end
set(ax, 'XLim',[1 nSessions], 'YLim', yl);

% background shading (Early green, Late yellow; match your convention)
patch(ax, [1 boundaryAt boundaryAt 1], [yl(1) yl(1) yl(2) yl(2)], ...
    [0 1 0], 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(ax, [boundaryAt nSessions nSessions boundaryAt], [yl(1) yl(1) yl(2) yl(2)], ...
    [1 1 0], 'FaceAlpha', bgAlpha, 'EdgeColor','none');

% null horizontal line: pooled null mean
nullMeanPooled = mean(splitR_null_all, 'omitnan');

if isfinite(nullMeanPooled)
    yl2 = ylim(ax);
    if nullMeanPooled < yl2(1), yl2(1) = nullMeanPooled - 0.02; end
    if nullMeanPooled > yl2(2), yl2(2) = nullMeanPooled + 0.02; end
    ylim(ax, yl2);
end

if isfinite(nullMeanPooled)
    yline(ax, nullMeanPooled, '--', 'LineWidth', eventLineLW, 'Color', [0.35 0.35 0.35]);
end

plot(ax, 1:nSessions, splitR, '-', 'LineWidth', dataLW, 'Color', colLine);
plot(ax, 1:nSessions, splitR, 'o', ...
    'LineStyle','none', 'MarkerSize', dataMS, ...
    'MarkerFaceColor', colLine, 'MarkerEdgeColor', colLine, 'LineWidth', markerLW);

xline(ax, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Split-half r (even vs odd)', 'FontSize', labelFontSize);
title(ax, 'Within-session template reliability', 'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off');
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Manual tick drawing (labels at 1,6,11,16,... ; long ticks when labeled)
ax.XTick = 1:nSessions;
labIdx = unique([1 6 11:5:nSessions]);
ax.XTickLabel = repmat({''}, 1, nSessions);
set(ax, 'TickLength', [0 0]);

y0 = ax.YLim(1);
yr = diff(ax.YLim);
majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

for s = 1:nSessions
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
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', tickFontSize, 'Color', [0 0 0], 'Clipping','off');
end

xCenter = mean(ax.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlab, 'Units','data', 'Position',[xCenter, yXlab, 0]);

ax.Position = [0.12 0.22 0.84 0.70];

outBaseRel = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_%02dSessions', nSessions));
exportgraphics(figRel, [outBaseRel '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figRel, [outBaseRel '.svg']); catch, print(figRel, [outBaseRel '.svg'], '-dsvg'); end
savefig(figRel, [outBaseRel '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseRel);

%% =========================
% Performance (% correct) vs split-half r scatter + straight line overlay
%% =========================
outScatPng = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_VS_PERFORMANCE_SCATTER_%02dSessions.png', nSessions));
outScatSvg = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_VS_PERFORMANCE_SCATTER_%02dSessions.svg', nSessions));
outScatFig = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_VS_PERFORMANCE_SCATTER_%02dSessions.fig', nSessions));

xPerf = perfObs(okP);
yR    = splitR(okP);

figPosSc = figPosDrift;
figPosSc(3) = round(figPosDrift(3) * 0.55);
figSc = figure('Color','w','Position',figPosSc);
axS = axes(figSc); hold(axS,'on');

plot(axS, xPerf, yR, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);

hasLine = numel(xPerf) >= 2;
if hasLine
    pLin  = polyfit(xPerf, yR, 1);
    xLine = linspace(min(xPerf), max(xPerf), 200);
    yLine = polyval(pLin, xLine);
    plot(axS, xLine, yLine, '-', 'LineWidth', axesTickLW, 'Color', colLine);
end

xlabel(axS, '% correct', 'FontSize', labelFontSize);
ylabel(axS, 'Split-half r (even vs odd)', 'FontSize', labelFontSize);
title(axS, sprintf('Performance vs split-half r (Spearman \\rho=%.2f, p=%.3g)', rhoPerf, pPerf), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axS, 'FontSize', tickFontSize, 'Box','off');
set(axS, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% limits with padding
if ~isempty(xPerf)
    xMin = min(xPerf); xMax = max(xPerf);
    xPad = 0.06 * max(eps, (xMax-xMin));
    set(axS, 'XLim', [xMin-xPad, xMax+xPad]);
end
if ~isempty(yR)
    yMin = min(yR); yMax = max(yR);
    yPad = 0.08 * max(eps, (yMax-yMin));
    set(axS, 'YLim', [max(-1, yMin-yPad), min(1, yMax+yPad)]);
end

axS.Position = [0.15 0.18 0.80 0.74];

exportgraphics(figSc, outScatPng, 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figSc, outScatSvg); catch, print(figSc, outScatSvg, '-dsvg'); end
savefig(figSc, outScatFig);

fprintf('Saved: %s\n', outScatPng);
fprintf('Saved: %s\n', outScatSvg);

%% =========================
% PLOT 2: Early vs Late bar + permutation test (two-sample)
%% =========================
earlyGroupUse = earlyGroup(earlyGroup>=1 & earlyGroup<=nSessions);
lateGroupUse  = lateGroup(lateGroup>=1  & lateGroup<=nSessions);

earlyVals = splitR(earlyGroupUse); earlyVals = earlyVals(isfinite(earlyVals));
lateVals  = splitR(lateGroupUse);  lateVals  = lateVals(isfinite(lateVals));

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
ylabel(axB, 'Split-half r', 'FontSize', labelFontSize);
title(axB, 'Within-session reliability', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axB, 'FontSize', tickFontSize, 'Box','off');
set(axB, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

if ~isempty(earlyVals) && ~isempty(lateVals)
    Tobs_el = mean(lateVals,'omitnan') - mean(earlyVals,'omitnan');
    pooled = [earlyVals(:); lateVals(:)];
    nE = numel(earlyVals);
    nL = numel(lateVals);

    Tperm = nan(Nperm,1);
    for ip = 1:Nperm
        permIdx = randperm(nE + nL);
        eIdx = permIdx(1:nE);
        lIdx = permIdx(nE+1:end);
        Tperm(ip) = mean(pooled(lIdx),'omitnan') - mean(pooled(eIdx),'omitnan');
    end

    p_t = (1 + sum(abs(Tperm) >= abs(Tobs_el))) / (1 + Nperm);

    addSigIfNeeded_(axB, 1, 2, p_t, [mE mL], axesTickLW, tickFontSize);
    fprintf('\n===== Early vs Late permutation test (split-half r) =====\n');
    fprintf('Early sessions: %s (n=%d)\n', mat2str(earlyGroupUse), numel(earlyVals));
    fprintf('Late  sessions: %s (n=%d)\n', mat2str(lateGroupUse),  numel(lateVals));
    fprintf('Permutation test (two-sided, N=%d): p=%.6g\n', Nperm, p_t);
end

outBaseBar = fullfile(outDir, sprintf('SPLITHALF_RELIABILITY_EARLY_LATE_BAR_%02dSessions', nSessions));
exportgraphics(figBar, [outBaseBar '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figBar, [outBaseBar '.svg']); catch, print(figBar, [outBaseBar '.svg'], '-dsvg'); end
savefig(figBar, [outBaseBar '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseBar);

fprintf('\nDONE.\n');

%% ================= HELPERS (same as your pipeline) =================

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