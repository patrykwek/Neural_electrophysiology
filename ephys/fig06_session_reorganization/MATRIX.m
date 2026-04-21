%% ========= FIG4_RSA_SESSION_SIM__UNITMATCHED__GLOBALWARP__CLASSIFIEDONLY.m =========
% UPDATED Figure 4 pipeline using RSA-style "temporal organization" similarity:
%
% For each session s:
%   1) Build per-unit warped z-SDF templates (EVEN trials) using your SAME global-warp pipeline.
%   2) Compute time×time correlation across units:
%         C_s(t1,t2) = corr( units(:,t1), units(:,t1) )
%   3) Vectorize upper triangle (excluding diagonal) => v_s
%
% Session×session similarity:
%   SIM(i,j) = corr(v_i, v_j)
%
% UNIT-MATCHED CONTROL (REQUESTED CHANGE):
%   - Keep Nmin computed the same way (from medians), BUT treat it as a CEILING.
%   - EVERY session is included in each graph/analysis (as long as it has >=1 classified unit).
%   - For sessions with > Nmin units: randomly subsample down to Nmin
%   - For sessions with <= Nmin units: use ALL available units (no exclusion)
%   - Repeat nBoot times, recompute C_s, v_s, SIM, and summaries
%
% OUTPUT PANELS (saved separately + .mat):
%   Fig 4A: Session×Session similarity matrix (unit-capped, RSA definition)
%   Fig 4B: EE / LL / EL distributions (hist overlay + bar±SEM + stars)
%   Fig 4C: Drift to Session 1 (unit-capped SIM) + bootstrap CI per session   <-- INCLUDED
%   Fig 4D: MDS embedding (unit-capped SIM) + centroid permutation test
%   Fig 4E: Example/avg time×time matrices (Early avg vs Late avg) with event lines
%
% STYLE CHANGE (requested):
%   - Early = YELLOW, Late = GREEN everywhere (background shading + bars)
%     (i.e., left side yellow, right side green)
%
% ADDITIONAL REQUEST (THIS EDIT ONLY):
%   - Fig 4B histogram x-axis and barplot y-axis are set to [0 0.7]
%
% NEW REQUEST (THIS EDIT ONLY):
%   - Time axes in Fig 4E are labeled in SECONDS (tick labels only; data unchanged).
%   - Event line colors: cue=RED, press=BLUE, lick=GREEN.
%   - Simplify axis labels and titles (keep accurate).
%
% NEW (REQUESTED NOW):
%   - Plot Even–Odd RSA descriptor reliability within session (noise ceiling for v_s)
%   - Plot Local drift (SIM to previous session) with bootstrap CI, styled like drift-to-session-1
%
% NEW (REQUESTED NOW):
%   - Drift vs learning companion FOR LOCAL DRIFT ONLY:
%       x = performance (% correct), y = local drift
%       scatter + fit line + Spearman on plot
%       binned (6–8 bins) mean ± bootstrap CI per bin
%
% NOTE:
%   - Keeps your data path conventions and global warp logic.
%   - Avoids unnecessary cosmetic changes beyond the requested edits.
%   - Sessions with 0 usable classified units remain in the matrices/plots as NaNs.

clear; clc;

%% ---------------- USER SETTINGS ----------------
matFile      = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/FIG4_RSA_UNITMATCHED';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Windowing (same conventions)
preCueMs   = 500;
postLickMs = 3000;
fixedWin   = [-1000, 6000];   % trial inclusion filters (relative to cue)

% Filters
MIN_RT_MS        = 100;
minTrialsPerUnit = 10;
requireValid     = true;

% SDF
dt_ms        = 10;
gaussSigmaMs = 25;

% Sessions
maxSessions = 36;

% Stage groups
earlyGroup = 1:8;
lateGroup  = 29:36;
boundaryAt = 8.5;

% Unit-matched bootstrap
nBoot = 300;         % 100–500 recommended
rngSeed = 0;

% ---------------- STYLE ----------------
figPosMat   = [120, 120, 980, 900];
figPosB     = [200, 200, 1180, 520];
figPosDrift = [120, 120, 1200, 650];
figPosMDS   = [120, 120, 1100, 650];
figPosE     = [120, 120, 1400, 650];

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;
axesTickLW    = 4.0;
axesTickLen   = [0.02 0.02];

dataLW      = 5;
dataMS      = 10;
markerLW    = 2.5;
eventLineLW = 5;

bgAlpha   = 0.10;
dashColor = [0.2 0.2 0.2];
dashLW    = 5;
dashStyle = '--';

% Early/Late colors
colEarlyShade = [1 1 0];   % EARLY = YELLOW
colLateShade  = [0 1 0];   % LATE  = GREEN
colMid_MDS     = [0.88 0.88 0.88];
colDriftLine   = [0 0.45 0.75];

% Similarity heatmap settings
simInterp = 'nearest';
simCLim   = [-0.2 0.6];
simCmapN  = 256;

boundaryLW    = 1.2;
boundaryColor = [0 0 0];
boundaryAlpha = 0.35;

% Fig4B colors (keep your EE/LL/EL scheme)
colEE = [0.10 0.55 0.95];
colLL = [0.90 0.20 0.10];
colEL = [0.95 0.70 0.10];

% Fig4B plot
nBins = 30;
histLineLW   = 6;
histAxesTickLW = 8;
errLineLW    = 4;
sigLineLW    = 4;
ptSize       = 32;

% Bar scatter style (match your decoding-stage)
scatterMS     = 8;
scatterEdgeLW = 1.0;
jit           = 0.10;
greyDot       = [0.6 0.6 0.6];

% Weibull-style x ticks
majorLenFrac      = 0.060;
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;
xLabelDownFrac    = 0.08;

%% ---------------- LOAD DATA ----------------
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT.');
nSessions = min(numel(beh), maxSessions);
beh = beh(1:nSessions);

assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---------------- GLOBAL WARP TARGETS ----------------
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

%% ---------------- TEMPLATE AXIS ----------------
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

%% ============================================================
% STEP 1) Build and cache per-session unit template matrices
%   unitTemplates{s}    = [nUnits_s x nT] (EVEN-trial mean per unit)
%   unitTemplatesOdd{s} = [nUnits_s x nT] (ODD-trial mean per unit)   <-- ADDED
%   nUnits(s) = nUnits_s
%% ============================================================
unitTemplates = cell(nSessions,1);
unitTemplatesOdd = cell(nSessions,1);  % ADDED
nUnits = zeros(nSessions,1);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Behavior-only keep trials (same as your pipeline) ----
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

    unitCurvesEven = []; % nUnitsThisSession x nT
    unitCurvesOdd  = []; % nUnitsThisSession x nT  (ADDED)

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
            if isempty(idxEven), idxEven = 1:nTrU; end
            if isempty(idxOdd),  idxOdd  = 1:nTrU; end

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

            muEven = mean(warpedZ(idxEven,:), 1);
            muOdd  = mean(warpedZ(idxOdd,:),  1);

            unitCurvesEven(end+1,:) = muEven; %#ok<AGROW>
            unitCurvesOdd(end+1,:)  = muOdd;  %#ok<AGROW>
        end
    end

    if ~isempty(unitCurvesEven)
        unitTemplates{sIdx} = unitCurvesEven;
        unitTemplatesOdd{sIdx} = unitCurvesOdd; % ADDED
        nUnits(sIdx) = size(unitCurvesEven,1);
    end
end

fprintf('\nBuilt unit templates for %d/%d sessions.\n', nnz(nUnits>0), nSessions);

%% ============================================================
% STEP 2) Define Nmin for unit-matched control (NOW USED AS A CEILING)
%% ============================================================
earlyUse = intersect(earlyGroup, find(nUnits>0)');
lateUse  = intersect(lateGroup,  find(nUnits>0)');

assert(~isempty(earlyUse) && ~isempty(lateUse), 'Need valid Early and Late sessions with units.');

medE = median(nUnits(earlyUse));
medL = median(nUnits(lateUse));
Nmin = min(medE, medL);
Nmin = floor(Nmin);

fprintf('\n=== UNIT CAPPING (Nmin used as CEILING) ===\n');
fprintf('Median units Early: %.1f | Late: %.1f | Ncap(Nmin)=%d\n', medE, medL, Nmin);

% include ALL sessions with >=1 unit (no >=Nmin floor exclusion)
sessHasN = find(nUnits > 0);
fprintf('Sessions with >=1 classified unit: %d/%d\n', numel(sessHasN), nSessions);

%% ============================================================
% ADDED: Even–Odd RSA descriptor reliability within session (noise ceiling for v_s)
%   For each session s, compute v_even and v_odd from the same capped units,
%   then reliability r_desc(s) = corr(v_even, v_odd)
%% ============================================================
maskUT = triu(true(nT), 1);
r_desc = nan(nSessions,1);

rng(rngSeed); % keep deterministic

for sIdx = sessHasN(:)'
    Ue = unitTemplates{sIdx};
    Uo = unitTemplatesOdd{sIdx};

    if isempty(Ue) || isempty(Uo), continue; end
    if size(Ue,1) ~= size(Uo,1)
        % should not happen if unit order is consistent; skip conservatively
        continue;
    end

    if size(Ue,1) > Nmin
        idx = randperm(size(Ue,1), Nmin);
        UeS = Ue(idx,:);
        UoS = Uo(idx,:);
    else
        UeS = Ue;
        UoS = Uo;
    end

    Ce = corr(UeS, 'Rows','pairwise');
    Co = corr(UoS, 'Rows','pairwise');

    ve = Ce(maskUT);
    vo = Co(maskUT);

    r_desc(sIdx) = corr(ve(:), vo(:), 'Rows','complete');
end

fprintf('\n=== Even–Odd RSA descriptor reliability (within-session) ===\n');
fprintf('Sessions with finite r_desc: %d/%d\n', nnz(isfinite(r_desc)), nSessions);
if any(isfinite(r_desc))
    fprintf('Mean r_desc = %.4f | Median r_desc = %.4f\n', mean(r_desc,'omitnan'), median(r_desc,'omitnan'));
end

% Plot descriptor reliability vs session (styled like drift)
figRdesc = figure('Color','w','Position',figPosDrift);
axR = axes(figRdesc); hold(axR,'on');

finiteR = r_desc(isfinite(r_desc));
if isempty(finiteR)
    ylR = [-1 1];
else
    pad = 0.05;
    ylR = [max(-1, min(finiteR)-pad), min(1, max(finiteR)+pad)];
end
set(axR, 'XLim',[1 nSessions], 'YLim', ylR);

patch(axR, [1 boundaryAt boundaryAt 1], [ylR(1) ylR(1) ylR(2) ylR(2)], ...
    colEarlyShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(axR, [boundaryAt nSessions nSessions boundaryAt], [ylR(1) ylR(1) ylR(2) ylR(2)], ...
    colLateShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');

plot(axR, 1:nSessions, r_desc, '-', 'LineWidth', dataLW, 'Color', colDriftLine);
plot(axR, 1:nSessions, r_desc, 'o', 'LineStyle','none', 'MarkerSize', dataMS, ...
    'MarkerFaceColor', colDriftLine, 'MarkerEdgeColor', colDriftLine, 'LineWidth', markerLW);

xline(axR, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlabR = xlabel(axR, 'Session', 'FontSize', labelFontSize);
ylabel(axR, 'r (v_{even} vs v_{odd})', 'FontSize', labelFontSize);
title(axR, 'RSA descriptor reliability', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axR, 'FontSize', tickFontSize, 'Box','off');
set(axR, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

axR.XTick = 1:nSessions;
labIdx = unique([1 6 11:5:nSessions]);
axR.XTickLabel = repmat({''}, 1, nSessions);
set(axR, 'TickLength', [0 0]);

y0 = axR.YLim(1);
yr = diff(axR.YLim);
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
    line(axR, [s s], [y0, y0 - L], 'Color', [0 0 0], 'LineWidth', lw, 'Clipping','off');
end

yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
for s = labIdx
    text(axR, s, yTickText, num2str(s), ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', tickFontSize, 'Color', [0 0 0], 'Clipping','off');
end

xCenter = mean(axR.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlabR, 'Units','data', 'Position',[xCenter, yXlab, 0]);
axR.Position = [0.12 0.22 0.84 0.70];

outBaseR = fullfile(outDir, sprintf('FIG4_RSA_DESCRIPTOR_RELIABILITY_EVEN_ODD_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figRdesc, [outBaseR '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figRdesc, [outBaseR '.svg']); catch, print(figRdesc, [outBaseR '.svg'], '-dsvg'); end
savefig(figRdesc, [outBaseR '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseR);

%% ============================================================
% STEP 3) Bootstrap: cap high-unit sessions at Nmin and recompute SIM
%% ============================================================
rng(rngSeed);

nV = nnz(maskUT);

SIM_boot = nan(nSessions, nSessions, nBoot);

% Also store bootstrap group means for quick diagnostic
muEE_boot = nan(nBoot,1);
muLL_boot = nan(nBoot,1);
muEL_boot = nan(nBoot,1);

% Precompute session-pair index sets (in original session numbering)
earlyMatch = intersect(earlyGroup, sessHasN');
lateMatch  = intersect(lateGroup,  sessHasN');

pairsEE_sess = nchoosek(earlyMatch, 2);
pairsLL_sess = nchoosek(lateMatch,  2);
[gridE, gridL] = ndgrid(earlyMatch, lateMatch);
pairsEL_sess = [gridE(:), gridL(:)];

fprintf('Included Early sessions used=%s\n', mat2str(earlyMatch));
fprintf('Included Late  sessions used=%s\n', mat2str(lateMatch));

for b = 1:nBoot
    V = nan(nV, nSessions);

    for sIdx = sessHasN(:)'
        U = unitTemplates{sIdx};
        if isempty(U), continue; end

        % cap only if too many units; otherwise use all
        if size(U,1) > Nmin
            idx = randperm(size(U,1), Nmin);
            Us = U(idx,:);
        else
            Us = U;
        end

        C = corr(Us, 'Rows','pairwise');   % nT x nT
        v = C(maskUT);

        V(:, sIdx) = v(:);
    end

    SIMb = corr(V, 'Rows','pairwise');     % nSessions x nSessions
    SIM_boot(:,:,b) = SIMb;

    % group means this bootstrap
    rEE = nan(size(pairsEE_sess,1),1);
    for k = 1:size(pairsEE_sess,1)
        i = pairsEE_sess(k,1); j = pairsEE_sess(k,2);
        rEE(k) = SIMb(i,j);
    end
    rLL = nan(size(pairsLL_sess,1),1);
    for k = 1:size(pairsLL_sess,1)
        i = pairsLL_sess(k,1); j = pairsLL_sess(k,2);
        rLL(k) = SIMb(i,j);
    end
    rEL = nan(size(pairsEL_sess,1),1);
    for k = 1:size(pairsEL_sess,1)
        i = pairsEL_sess(k,1); j = pairsEL_sess(k,2);
        rEL(k) = SIMb(i,j);
    end

    muEE_boot(b) = mean(rEE, 'omitnan');
    muLL_boot(b) = mean(rLL, 'omitnan');
    muEL_boot(b) = mean(rEL, 'omitnan');
end

SIM = mean(SIM_boot, 3, 'omitnan');   % point estimate used for plotting

fprintf('\n=== BOOTSTRAP SUMMARY (means across session-pairs) ===\n');
fprintf('EE mean: %.4f (2.5%%=%.4f, 97.5%%=%.4f)\n', mean(muEE_boot,'omitnan'), prctile(muEE_boot,2.5), prctile(muEE_boot,97.5));
fprintf('LL mean: %.4f (2.5%%=%.4f, 97.5%%=%.4f)\n', mean(muLL_boot,'omitnan'), prctile(muLL_boot,2.5), prctile(muLL_boot,97.5));
fprintf('EL mean: %.4f (2.5%%=%.4f, 97.5%%=%.4f)\n', mean(muEL_boot,'omitnan'), prctile(muEL_boot,2.5), prctile(muEL_boot,97.5));

simFinite = SIM(isfinite(SIM));
if ~isempty(simFinite)
    fprintf('\n=== SIM MATRIX RANGE ===\n');
    fprintf('SIM min = %.4f | SIM max = %.4f\n', min(simFinite), max(simFinite));
end

%% Save core results
outMat = fullfile(outDir, sprintf('FIG4_RSA_UNITMATCHED_SIM_%02dSessions_Nmin%d_nBoot%d.mat', nSessions, Nmin, nBoot));
save(outMat, 'SIM','SIM_boot','muEE_boot','muLL_boot','muEL_boot', ...
    'nUnits','unitTemplates','Nmin','nBoot','sessHasN', ...
    'earlyGroup','lateGroup','earlyMatch','lateMatch', ...
    'tgrid_ms','Tpress','Tlick','dt_ms','preCueMs','postLickMs', '-v7.3');
fprintf('Saved: %s\n', outMat);

%% ============================================================
% FIG 4A — Session × Session similarity matrix (UPDATED, unit-capped)
%% ============================================================
figA = figure('Color','w','Position',figPosMat);
axA = axes(figA); hold(axA,'on');

hImg = imagesc(axA, 1:nSessions, 1:nSessions, SIM);
set(axA,'YDir','reverse');

try
    set(hImg, 'Interpolation', simInterp);
catch
    set(hImg, 'Interpolation', 'nearest');
end

colormap(axA, plasma(simCmapN));
caxis(axA, simCLim);

axis(axA,'image');
xlim(axA,[0.5 nSessions+0.5]);
ylim(axA,[0.5 nSessions+0.5]);

xlabel(axA, 'Session', 'FontSize', labelFontSize);
ylabel(axA, 'Session', 'FontSize', labelFontSize);
title(axA, 'Session similarity', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

% Ticks: major every 5
major = 1:5:nSessions;
if major(end) ~= nSessions, major = unique([major nSessions]); end
xticks(axA, major); yticks(axA, major);
xticklabels(axA, string(major)); yticklabels(axA, string(major));
axA.XAxis.MinorTick = 'on';
axA.YAxis.MinorTick = 'on';
axA.XAxis.MinorTickValues = setdiff(1:nSessions, major);
axA.YAxis.MinorTickValues = setdiff(1:nSessions, major);

set(axA, 'FontSize', tickFontSize, 'TickDir','out', 'LineWidth', axesTickLW, ...
    'TickLength', axesTickLen, 'Box','off', 'Layer','top');

% Boundaries between sessions only, special dashed at 8.5
for b = 1.5:1:(nSessions-0.5)

    thisLW    = boundaryLW;
    thisCol   = boundaryColor;
    thisStyle = '-';
    thisAlpha = boundaryAlpha;

    if abs(b - boundaryAt) < 1e-9
        thisLW    = dashLW;
        thisCol   = dashColor;
        thisStyle = dashStyle;
        thisAlpha = 1.0;
    end

    lh = line(axA, [0.5 nSessions+0.5], [b b], 'Color', thisCol, 'LineWidth', thisLW, 'LineStyle', thisStyle);
    lv = line(axA, [b b], [0.5 nSessions+0.5], 'Color', thisCol, 'LineWidth', thisLW, 'LineStyle', thisStyle);

    if ~(abs(b - boundaryAt) < 1e-9)
        try, lh.Color(4) = thisAlpha; lv.Color(4) = thisAlpha; catch, end
    end
end

cb = colorbar(axA, 'eastoutside');
cb.TickDirection = 'out';
cb.LineWidth = axesTickLW;
cb.Label.String = 'r';
cb.FontSize = 22;
cb.Label.FontSize = 22;
cb.Ticks = -0.2:0.2:0.6;

outBaseA = fullfile(outDir, sprintf('FIG4A_RSA_UNITMATCHED_SIM_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figA, [outBaseA '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figA, [outBaseA '.svg']); catch, print(figA, [outBaseA '.svg'], '-dsvg'); end
savefig(figA, [outBaseA '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseA);

%% ============================================================
% FIG 4B — Within-stage vs across-stage (EE/LL/EL) from unit-capped SIM
%% ============================================================
% Extract session-pair entries from SIM (point estimate)
r_EE = nan(size(pairsEE_sess,1),1);
for k = 1:size(pairsEE_sess,1)
    i = pairsEE_sess(k,1); j = pairsEE_sess(k,2);
    r_EE(k) = SIM(i,j);
end
r_EE = r_EE(isfinite(r_EE));

r_LL = nan(size(pairsLL_sess,1),1);
for k = 1:size(pairsLL_sess,1)
    i = pairsLL_sess(k,1); j = pairsLL_sess(k,2);
    r_LL(k) = SIM(i,j);
end
r_LL = r_LL(isfinite(r_LL));

r_EL = nan(size(pairsEL_sess,1),1);
for k = 1:size(pairsEL_sess,1)
    i = pairsEL_sess(k,1); j = pairsEL_sess(k,2);
    r_EL(k) = SIM(i,j);
end
r_EL = r_EL(isfinite(r_EL));

% ------------------- CHANGED: permutation tests instead of ranksum -------------------
nPermBar = 5000;
[p_LL_EE, ~] = permTestDiffMeans_(r_LL, r_EE, nPermBar);
[p_LL_EL, ~] = permTestDiffMeans_(r_LL, r_EL, nPermBar);
[p_EE_EL, ~] = permTestDiffMeans_(r_EE, r_EL, nPermBar);
% -----------------------------------------------------------------------------------

fprintf('\n=== Permutation tests on SIM session-pair entries (unit-capped SIM point estimate) ===\n');
fprintf('LL vs EE : p = %.6g\n', p_LL_EE);
fprintf('LL vs EL : p = %.6g\n', p_LL_EL);
fprintf('EE vs EL : p = %.6g\n', p_EE_EL);

figB = figure('Color','w','Position',figPosB);
tl = tiledlayout(figB, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% Left: histogram overlay
ax1 = nexttile(tl, 1); hold(ax1,'on');

% axis range 0 to 0.7
edges = linspace(0, 0.7, nBins+1);

hEE = histogram(ax1, r_EE, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);
hLL = histogram(ax1, r_LL, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);
hEL = histogram(ax1, r_EL, edges, 'Normalization','probability', 'DisplayStyle','stairs', 'LineWidth',histLineLW);

hEE.EdgeColor = colEE;
hLL.EdgeColor = colLL;
hEL.EdgeColor = colEL;

xlabel(ax1, 'r', 'FontSize', labelFontSize);
ylabel(ax1, 'Probability', 'FontSize', labelFontSize);
title(ax1, 'Pairwise similarity', 'FontSize', titleFontSize, 'FontWeight','bold');

legend(ax1, {sprintf('Early–Early (%d-%d)', earlyGroup(1), earlyGroup(end)), ...
             sprintf('Late–Late (%d-%d)',  lateGroup(1),  lateGroup(end)), ...
             'Early–Late'}, 'Location','northwest');

xlim(ax1, [0 0.7]);
set(ax1, 'FontSize', tickFontSize, 'TickDir','out', 'LineWidth', histAxesTickLW, ...
    'TickLength', axesTickLen, 'Box','off');

% Right: bar ± SEM + points + stars
ax2 = nexttile(tl, 2); hold(ax2,'on');

cats = ["Early–Early","Late–Late","Early–Late"];
xCats = categorical(cats, cats, 'Ordinal', true);

mu  = [mean(r_EE,'omitnan'), mean(r_LL,'omitnan'), mean(r_EL,'omitnan')];
nEE = numel(r_EE); nLL = numel(r_LL); nEL = numel(r_EL);
sem = [std(r_EE,0,'omitnan')/sqrt(max(1,nEE)), ...
       std(r_LL,0,'omitnan')/sqrt(max(1,nLL)), ...
       std(r_EL,0,'omitnan')/sqrt(max(1,nEL))];

bbar = bar(ax2, xCats, mu);
bbar.FaceColor = 'flat';
bbar.CData(1,:) = colEE;
bbar.CData(2,:) = colLL;
bbar.CData(3,:) = colEL;

errorbar(ax2, bbar.XEndPoints, mu, sem, 'k.', 'LineWidth', errLineLW);

rng(0);
jitter = 0.12;

R = {r_EE, r_LL, r_EL};
nShow = [min(nEE,3000), min(nLL,3000), min(nEL,3000)];
idxShow = cell(1,3);
idxShow{1} = randperm(nEE, nShow(1));
idxShow{2} = randperm(nLL, nShow(2));
idxShow{3} = randperm(nEL, nShow(3));

for i = 1:3
    rr = R{i}(idxShow{i});
    xj = bbar.XEndPoints(i) + (rand(numel(rr),1)-0.5)*2*jitter;
    scatter(ax2, xj, rr, ptSize, 'filled', ...
        'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerFaceAlpha', 0.20, ...
        'MarkerEdgeAlpha', 0);
end

ylabel(ax2, 'r', 'FontSize', labelFontSize);
xlabel(ax2, 'Group', 'FontSize', labelFontSize);
title(ax2, 'Group mean', 'FontSize', titleFontSize, 'FontWeight','bold');

ylim(ax2, [0 0.7]);

set(ax2, 'FontSize', tickFontSize, 'LineWidth', axesTickLW, 'TickDir','out');
box(ax2,'off');

% Significance bars + stars (only if significant; your rule)
pairs = [1 2; 2 3; 1 3];
pVals = [p_LL_EE; p_LL_EL; p_EE_EL];
sig = pVals < 0.05;

if any(sig)
    yMax = max([r_EE(:); r_LL(:); r_EL(:)], [], 'omitnan');
    yMin = min([r_EE(:); r_LL(:); r_EL(:)], [], 'omitnan');
    yRange = yMax - yMin; if ~isfinite(yRange) || yRange==0, yRange=1; end

    baseY = max(mu + sem) + 0.06*yRange;
    hBar  = 0.02*yRange;
    stepY = 0.06*yRange;

    drawCount = 0;
    for kk = 1:3
        if ~sig(kk), continue; end
        drawCount = drawCount + 1;

        xA = bbar.XEndPoints(pairs(kk,1));
        xB = bbar.XEndPoints(pairs(kk,2));
        yLevel = baseY + (drawCount-1)*stepY;

        plot(ax2, [xA xA xB xB], [yLevel yLevel+hBar yLevel+hBar yLevel], 'k-', 'LineWidth', sigLineLW);

        stars = pToStars_(pVals(kk));
        text(ax2, mean([xA xB]), yLevel+hBar+0.002*yRange, stars, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', titleFontSize, 'FontWeight','bold');
    end
end

outBaseB = fullfile(outDir, sprintf('FIG4B_RSA_UNITMATCHED_EE_LL_EL_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figB, [outBaseB '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figB, [outBaseB '.svg']); catch, print(figB, [outBaseB '.svg'], '-dsvg'); end
savefig(figB, [outBaseB '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseB);

%% ============================================================
% FIG 4C — Drift relative to Session 1 (unit-capped SIM) + per-session bootstrap CI
%% ============================================================
drift = SIM(:,1);

% --- CI PER SESSION (2.5%–97.5%) FROM BOOTSTRAP ---
drift_boot = squeeze(SIM_boot(:,1,:));     % nSessions x nBoot
drift_lo   = prctile(drift_boot, 2.5,  2);
drift_hi   = prctile(drift_boot, 97.5, 2);

figC = figure('Color','w','Position',figPosDrift);
axC = axes(figC); hold(axC,'on');

finiteBand = [drift_lo(:); drift_hi(:)];
finiteBand = finiteBand(isfinite(finiteBand));
if isempty(finiteBand)
    yl = [-1 1];
else
    pad = 0.05;
    yl = [max(-1, min(finiteBand)-pad), min(1, max(finiteBand)+pad)];
end
set(axC, 'XLim',[1 nSessions], 'YLim', yl);

% Background shading (Early yellow, Late green)
patch(axC, [1 boundaryAt boundaryAt 1], [yl(1) yl(1) yl(2) yl(2)], ...
    colEarlyShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(axC, [boundaryAt nSessions nSessions boundaryAt], [yl(1) yl(1) yl(2) yl(2)], ...
    colLateShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');

% CI band
xBand = 1:nSessions;
okBand = isfinite(drift_lo) & isfinite(drift_hi);
if any(okBand)
    xx = [xBand(okBand), fliplr(xBand(okBand))];
    yy = [drift_lo(okBand)', fliplr(drift_hi(okBand)')];
    hBand = patch(axC, xx, yy, colDriftLine, 'EdgeColor','none');
    try, hBand.FaceAlpha = 0.18; catch, end
end

plot(axC, 1:nSessions, drift, '-', 'LineWidth', dataLW, 'Color', colDriftLine);
plot(axC, 1:nSessions, drift, 'o', 'LineStyle','none', 'MarkerSize', dataMS, ...
    'MarkerFaceColor', colDriftLine, 'MarkerEdgeColor', colDriftLine, 'LineWidth', markerLW);

xline(axC, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlab = xlabel(axC, 'Session', 'FontSize', labelFontSize);
ylabel(axC, 'r to Session 1', 'FontSize', labelFontSize);
title(axC, 'Drift', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axC, 'FontSize', tickFontSize, 'Box','off');
set(axC, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% Weibull-style ticks: labels at 1,6,11,16,...
axC.XTick = 1:nSessions;
labIdx = unique([1 6 11:5:nSessions]);
axC.XTickLabel = repmat({''}, 1, nSessions);
set(axC, 'TickLength', [0 0]);

y0 = axC.YLim(1);
yr = diff(axC.YLim);
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
    line(axC, [s s], [y0, y0 - L], 'Color', [0 0 0], 'LineWidth', lw, 'Clipping','off');
end

yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
for s = labIdx
    text(axC, s, yTickText, num2str(s), ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', tickFontSize, 'Color', [0 0 0], 'Clipping','off');
end

xCenter = mean(axC.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlab, 'Units','data', 'Position',[xCenter, yXlab, 0]);
axC.Position = [0.12 0.22 0.84 0.70];

outBaseC = fullfile(outDir, sprintf('FIG4C_RSA_UNITMATCHED_DRIFT_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figC, [outBaseC '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figC, [outBaseC '.svg']); catch, print(figC, [outBaseC '.svg'], '-dsvg'); end
savefig(figC, [outBaseC '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseC);

%% ============================================================
% ADDED: Local drift (similarity to previous session) + bootstrap CI
%% ============================================================
localDrift = nan(nSessions,1);
localDrift(2:end) = diag(SIM, -1);

localDrift_boot = nan(nSessions, nBoot);
for b = 1:nBoot
    SIMb = SIM_boot(:,:,b);
    localDrift_boot(2:end, b) = diag(SIMb, -1);
end

localDrift_lo = prctile(localDrift_boot, 2.5,  2);
localDrift_hi = prctile(localDrift_boot, 97.5, 2);

figC2 = figure('Color','w','Position',figPosDrift);
axC2 = axes(figC2); hold(axC2,'on');

finiteBand2 = [localDrift_lo(:); localDrift_hi(:)];
finiteBand2 = finiteBand2(isfinite(finiteBand2));
if isempty(finiteBand2)
    yl2 = [-1 1];
else
    pad = 0.05;
    yl2 = [max(-1, min(finiteBand2)-pad), min(1, max(finiteBand2)+pad)];
end
set(axC2, 'XLim',[1 nSessions], 'YLim', yl2);

patch(axC2, [1 boundaryAt boundaryAt 1], [yl2(1) yl2(1) yl2(2) yl2(2)], ...
    colEarlyShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');
patch(axC2, [boundaryAt nSessions nSessions boundaryAt], [yl2(1) yl2(1) yl2(2) yl2(2)], ...
    colLateShade, 'FaceAlpha', bgAlpha, 'EdgeColor','none');

xBand = 1:nSessions;
okBand2 = isfinite(localDrift_lo) & isfinite(localDrift_hi);
if any(okBand2)
    xx = [xBand(okBand2), fliplr(xBand(okBand2))];
    yy = [localDrift_lo(okBand2)', fliplr(localDrift_hi(okBand2)')];
    hBand2 = patch(axC2, xx, yy, colDriftLine, 'EdgeColor','none');
    try, hBand2.FaceAlpha = 0.18; catch, end
end

plot(axC2, 1:nSessions, localDrift, '-', 'LineWidth', dataLW, 'Color', colDriftLine);
plot(axC2, 1:nSessions, localDrift, 'o', 'LineStyle','none', 'MarkerSize', dataMS, ...
    'MarkerFaceColor', colDriftLine, 'MarkerEdgeColor', colDriftLine, 'LineWidth', markerLW);

xline(axC2, boundaryAt, '--', 'LineWidth', eventLineLW, 'Color', dashColor);

hxlab2 = xlabel(axC2, 'Session', 'FontSize', labelFontSize);
ylabel(axC2, 'r to previous session', 'FontSize', labelFontSize);
title(axC2, 'Local drift', 'FontWeight','bold', 'FontSize', titleFontSize);

set(axC2, 'FontSize', tickFontSize, 'Box','off');
set(axC2, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

axC2.XTick = 1:nSessions;
axC2.XTickLabel = repmat({''}, 1, nSessions);
set(axC2, 'TickLength', [0 0]);

y0 = axC2.YLim(1);
yr = diff(axC2.YLim);
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
    line(axC2, [s s], [y0, y0 - L], 'Color', [0 0 0], 'LineWidth', lw, 'Clipping','off');
end

yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
for s = labIdx
    text(axC2, s, yTickText, num2str(s), ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', tickFontSize, 'Color', [0 0 0], 'Clipping','off');
end

xCenter = mean(axC2.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlab2, 'Units','data', 'Position',[xCenter, yXlab, 0]);
axC2.Position = [0.12 0.22 0.84 0.70];

outBaseC2 = fullfile(outDir, sprintf('FIG4C2_RSA_UNITMATCHED_LOCAL_DRIFT_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figC2, [outBaseC2 '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figC2, [outBaseC2 '.svg']); catch, print(figC2, [outBaseC2 '.svg'], '-dsvg'); end
savefig(figC2, [outBaseC2 '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseC2);

%% ============================================================
% ADDED: Local drift vs learning (performance % correct)
%   - scatter + linear fit
%   - Spearman printed + included in title
%   - bin performance into 7 bins; plot mean ± bootstrap CI per bin
%% ============================================================
% performance extraction (% correct)
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

okLP = isfinite(perfObs) & isfinite(localDrift);
xPerf = perfObs(okLP);
yLD   = localDrift(okLP);

if nnz(okLP) >= 3
    [rhoLD, pLD] = corr(perfObs, localDrift, 'Type','Spearman', 'Rows','complete');
else
    rhoLD = nan; pLD = nan;
end
fprintf('\n=== Spearman: local drift vs performance (%% correct) ===\n');
fprintf('n=%d sessions, rho=%.4f, p=%.6g\n', nnz(okLP), rhoLD, pLD);

% figure
figLP = figure('Color','w','Position',figPosDrift);
figLP.Position(3) = round(figPosDrift(3) * 0.62);
axLP = axes(figLP); hold(axLP,'on');

% scatter (ONLY grey scatter requested)
plot(axLP, xPerf, yLD, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);

% fit line (simple linear)
if numel(xPerf) >= 2 && all(isfinite(xPerf)) && all(isfinite(yLD))
    pLin  = polyfit(xPerf, yLD, 1);
    xLine = linspace(min(xPerf), max(xPerf), 200);
    yLine = polyval(pLin, xLine);
    plot(axLP, xLine, yLine, '-', 'LineWidth', axesTickLW, 'Color', colDriftLine);
end

xlabel(axLP, '% correct', 'FontSize', labelFontSize);
ylabel(axLP, 'Local drift (r to previous session)', 'FontSize', labelFontSize);
title(axLP, sprintf('Local drift vs learning (Spearman \\rho=%.2f, p=%.3g)', rhoLD, pLD), ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(axLP, 'FontSize', tickFontSize, 'Box','off');
set(axLP, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% limits with padding
if ~isempty(xPerf)
    xMin = min(xPerf); xMax = max(xPerf);
    xPad = 0.06 * max(eps, (xMax-xMin));
    set(axLP, 'XLim', [xMin-xPad, xMax+xPad]);
end
if ~isempty(yLD)
    yMin = min(yLD); yMax = max(yLD);
    yPad = 0.08 * max(eps, (yMax-yMin));
    set(axLP, 'YLim', [max(-1, yMin-yPad), min(1, yMax+yPad)]);
end

axLP.Position = [0.15 0.18 0.80 0.74];

outBaseLP = fullfile(outDir, sprintf('FIG4_LOCAL_DRIFT_VS_PERFORMANCE_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figLP, [outBaseLP '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figLP, [outBaseLP '.svg']); catch, print(figLP, [outBaseLP '.svg'], '-dsvg'); end
savefig(figLP, [outBaseLP '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseLP);

%% ============================================================
% FIG 4D — MDS embedding (unit-capped SIM)
% + permutation test on centroid separation (Early vs Late only)
%% ============================================================
R = SIM;
Rnan = isnan(R);
R(Rnan) = 0;
R = max(-1, min(1, R));

D = sqrt(max(0, 2*(1 - R)));
D = 0.5*(D + D');
D(1:nSessions+1:end) = 0;

if any(Rnan(:))
    medD = median(D(~Rnan), 'omitnan');
    D(Rnan) = medD;
    D = 0.5*(D + D');
    D(1:nSessions+1:end) = 0;
end

[Y, eigvals] = classical_mds_(D, 2);
x = Y(:,1); y = Y(:,2);

isEarly = ismember(1:nSessions, earlyGroup);
isLate  = ismember(1:nSessions, lateGroup);
isMid   = ~(isEarly | isLate);

% centroid permutation test (restrict to sessions with units)
Eidx = earlyMatch(:);
Lidx = lateMatch(:);
XY = [x(:), y(:)];
nPerm_MDS = 5000;
[centDist_obs, p_cent] = centroidPermTest2D_(XY, Eidx, Lidx, nPerm_MDS);

fprintf('\n=== MDS centroid permutation test (Early vs Late only) ===\n');
fprintf('Centroid dist = %.6f, p = %.6g (nPerm=%d)\n', centDist_obs, p_cent, nPerm_MDS);

figD = figure('Color','w','Position',figPosMDS);
axD = axes(figD); hold(axD,'on');

mdsMS = 2 * dataMS;

sessionVals = (1:nSessions)';
scatter(axD, x, y, mdsMS^2, sessionVals, 'filled', 'MarkerEdgeColor', 'none');

baseSummer = flipud(summer(256));
keepTo = round(0.78 * size(baseSummer,1));
mdsCmap = baseSummer(round(linspace(1, keepTo, nSessions)), :);
colormap(axD, mdsCmap);
caxis(axD, [1 nSessions]);

title(axD, 'MDS', 'FontWeight','bold', 'FontSize', titleFontSize);
xlabel(axD, 'Dim 1', 'FontSize', labelFontSize);
ylabel(axD, 'Dim 2', 'FontSize', labelFontSize);

set(axD, 'FontSize', tickFontSize, 'Box','off');
set(axD, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

cbD = colorbar(axD, 'eastoutside');
cbD.TickDirection = 'out';
cbD.LineWidth = axesTickLW;
cbD.Label.String = 'Session';
cbD.FontSize = 22;
cbD.Label.FontSize = 22;
cbD.Ticks = [1 nSessions];
cbD.TickLabels = {num2str(1), num2str(nSessions)};

if ~isempty(eigvals) && all(isfinite(eigvals)) && sum(max(eigvals,0))>0
    ve = 100 * max(eigvals,0) ./ sum(max(eigvals,0));
    txt = sprintf('Var: %.1f%%, %.1f%%', ve(1), ve(2));
    text(axD, 0.02, 0.98, txt, 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 18);
end

txtP = sprintf('Centroid perm E vs L: p=%.3g', p_cent);
text(axD, 0.02, 0.90, txtP, 'Units','normalized', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 18);

outBaseD = fullfile(outDir, sprintf('FIG4D_RSA_UNITMATCHED_MDS_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figD, [outBaseD '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figD, [outBaseD '.svg']); catch, print(figD, [outBaseD '.svg'], '-dsvg'); end
savefig(figD, [outBaseD '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseD);

%% ============================================================
% FIG 4E — Example/avg time×time matrices (Early avg vs Late avg), unit-capped
% Use one fixed subsample per session for a clean “average C” visualization.
%% ============================================================
rng(rngSeed);

Cearly_sum = zeros(nT,nT);
Clate_sum  = zeros(nT,nT);
nE = 0; nL = 0;

for sIdx = earlyMatch(:)'
    U = unitTemplates{sIdx};
    if isempty(U), continue; end

    if size(U,1) > Nmin
        idx = randperm(size(U,1), Nmin);
        C = corr(U(idx,:), 'Rows','pairwise');
    else
        C = corr(U, 'Rows','pairwise');
    end

    if ~any(isfinite(C(:))), continue; end

    Cearly_sum = Cearly_sum + C;
    nE = nE + 1;
end

for sIdx = lateMatch(:)'
    U = unitTemplates{sIdx};
    if isempty(U), continue; end

    if size(U,1) > Nmin
        idx = randperm(size(U,1), Nmin);
        C = corr(U(idx,:), 'Rows','pairwise');
    else
        C = corr(U, 'Rows','pairwise');
    end

    if ~any(isfinite(C(:))), continue; end

    Clate_sum = Clate_sum + C;
    nL = nL + 1;
end

Cearly = Cearly_sum / max(1,nE);
Clate  = Clate_sum  / max(1,nL);

figE = figure('Color','w','Position',figPosE);
tlE = tiledlayout(figE,1,2,'TileSpacing','compact','Padding','compact');

tCue   = 0;
tPress = Tpress;
tLick  = Tlick;

iCue   = timeToIdx_(tCue,   winLeft, dt_ms, nT);
iPress = timeToIdx_(tPress, winLeft, dt_ms, nT);
iLick  = timeToIdx_(tLick,  winLeft, dt_ms, nT);

tickStep_s  = 1; % 1 second
t0_s = winLeft/1000;
t1_s = winRight_template/1000;
tickVals_s = ceil(t0_s/tickStep_s)*tickStep_s : tickStep_s : floor(t1_s/tickStep_s)*tickStep_s;
tickPos = arrayfun(@(ts) timeToIdx_(ts*1000, winLeft, dt_ms, nT), tickVals_s);
tickLabs = arrayfun(@(ts) sprintf('%.0f', ts), tickVals_s, 'UniformOutput', false);

axE1 = nexttile(tlE,1); hold(axE1,'on');
imagesc(axE1, Cearly);
axis(axE1,'image'); set(axE1,'YDir','normal');
colormap(axE1, plasma(simCmapN));
caxis(axE1, simCLim);

title(axE1, sprintf('Early (sessions %d–%d)', earlyGroup(1), earlyGroup(end)), ...
    'FontWeight','bold','FontSize',titleFontSize);

xlabel(axE1,'Time (s)','FontSize',labelFontSize);
ylabel(axE1,'Time (s)','FontSize',labelFontSize);

set(axE1,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,'Box','off');
set(axE1,'XTick',tickPos,'XTickLabel',tickLabs,'YTick',tickPos,'YTickLabel',tickLabs);

drawEventLines_(axE1, iCue, iPress, iLick);

axE2 = nexttile(tlE,2); hold(axE2,'on');
imagesc(axE2, Clate);
axis(axE2,'image'); set(axE2,'YDir','normal');
colormap(axE2, plasma(simCmapN));
caxis(axE2, simCLim);

title(axE2, sprintf('Late (sessions %d–%d)', lateGroup(1), lateGroup(end)), ...
    'FontWeight','bold','FontSize',titleFontSize);

xlabel(axE2,'Time (s)','FontSize',labelFontSize);
ylabel(axE2,'Time (s)','FontSize',labelFontSize);

set(axE2,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,'Box','off');
set(axE2,'XTick',tickPos,'XTickLabel',tickLabs,'YTick',tickPos,'YTickLabel',tickLabs);

drawEventLines_(axE2, iCue, iPress, iLick);

cbE = colorbar(axE2,'eastoutside');
cbE.TickDirection = 'out';
cbE.LineWidth = axesTickLW;
cbE.Label.String = 'Correlation';
cbE.FontSize = 22;
cbE.Label.FontSize = 22;
cbE.Ticks = [0 0.25 0.5 0.75 1];

outBaseE = fullfile(outDir, sprintf('FIG4E_RSA_UNITMATCHED_CtXt_EarlyLate_%02dSessions_Nmin%d', nSessions, Nmin));
exportgraphics(figE, [outBaseE '.png'], 'Resolution', 300, 'BackgroundColor','white');
try, saveas(figE, [outBaseE '.svg']); catch, print(figE, [outBaseE '.svg'], '-dsvg'); end
savefig(figE, [outBaseE '.fig']);
fprintf('Saved: %s.[png/svg/fig]\n', outBaseE);

fprintf('\nDONE.\n');

%% ============================ HELPERS ============================

function drawEventLines_(ax, iCue, iPress, iLick)
    cueCol   = [1 0 0];
    pressCol = [0 0.45 0.85];
    lickCol  = [0 0.6 0];

    lw = 3.5;

    xline(ax, iCue,   '-', 'LineWidth', lw, 'Color', cueCol);
    yline(ax, iCue,   '-', 'LineWidth', lw, 'Color', cueCol);

    xline(ax, iPress, '-', 'LineWidth', lw, 'Color', pressCol);
    yline(ax, iPress, '-', 'LineWidth', lw, 'Color', pressCol);

    xline(ax, iLick,  '-', 'LineWidth', lw, 'Color', lickCol);
    yline(ax, iLick,  '-', 'LineWidth', lw, 'Color', lickCol);
end

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

function s = pToStars_(p)
    if p < 0.001, s='***';
    elseif p < 0.01, s='**';
    else, s='*';
    end
end

function [distObs, p] = centroidPermTest2D_(XY, idxE, idxL, nPerm)
    idxE = idxE(:);
    idxL = idxL(:);

    allIdx = [idxE; idxL];
    allIdx = allIdx(:);

    ok = all(isfinite(XY(allIdx,:)), 2);
    allIdx = allIdx(ok);

    idxE = intersect(idxE, allIdx);
    idxL = intersect(idxL, allIdx);

    nE = numel(idxE);
    nL = numel(idxL);
    assert(nE>0 && nL>0, 'Need non-empty Early and Late sets for centroid test.');

    muE = mean(XY(idxE,:), 1);
    muL = mean(XY(idxL,:), 1);
    distObs = norm(muE - muL);

    rng(0);
    Tnull = nan(nPerm,1);

    pool = [idxE; idxL];
    nP = numel(pool);

    for pp = 1:nPerm
        perm = pool(randperm(nP));
        pE = perm(1:nE);
        pL = perm(nE+1:end);

        muEp = mean(XY(pE,:), 1);
        muLp = mean(XY(pL,:), 1);
        Tnull(pp) = norm(muEp - muLp);
    end

    p = (1 + nnz(Tnull >= distObs)) / (nPerm + 1);
end

% ------------------- ADDED HELPER (minimal) -------------------
function [p, diffObs] = permTestDiffMeans_(a, b, nPerm)
    a = a(:); b = b(:);
    a = a(isfinite(a));
    b = b(isfinite(b));
    if isempty(a) || isempty(b)
        p = nan; diffObs = nan; return;
    end

    diffObs = mean(a,'omitnan') - mean(b,'omitnan');

    pool = [a; b];
    nA = numel(a);
    nP = numel(pool);

    rng(0);
    T = nan(nPerm,1);
    for i = 1:nPerm
        idx = randperm(nP);
        aP = pool(idx(1:nA));
        bP = pool(idx(nA+1:end));
        T(i) = mean(aP,'omitnan') - mean(bP,'omitnan');
    end

    p = (1 + nnz(abs(T) >= abs(diffObs))) / (nPerm + 1); % two-sided
end
% --------------------------------------------------------------