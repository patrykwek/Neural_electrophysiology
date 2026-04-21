%% ========= FIG3C_STAGE_ZSCORE_TRACES__A_THESIS_FINAL__GLOBALWARP__FIG2TRIALS_EXACT__CLASSIFIEDONLY.m =========
% Population z-score traces (Mizes Fig4E-style):
%   - Build SDF (Hz) per trial
%   - Warp cue->press and press->lick to GLOBAL targets (shared across stages)
%   - Z-score EACH TRIAL across time (WHOLE-TRIAL z; not baseline-z)
%   - (Optional) light smoothing AFTER warp (in warped time) to reduce edge sharpening
%   - Average across trials -> per-unit trace
%   - Average across units -> population trace (+/- SEM)
%
% Changes implemented now:
%   1) Force y-axis limits to [-1, +1] on BOTH produced plots (one per stage).
%      (Does NOT change how z-score is calculated.)
%   2) Match plot styling to FIG3A_STAGE_HEATMAPS:
%      - Match font sizes, axis linewidths, tick lengths, and event line width
%   3) Make x-axis tick labels "straight" (no rotation), like FIG3A.
%   4) Make the figure 1/3 of its previous height (keep y-lims [-1,1]).
%   No other unnecessary changes.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/4_2_ZSCORE';
outDir  = fullfile(baseOut, '4_2_ZSCORE');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Windowing (match FIG3A design)
preCueMs   = 500;
postLickMs = 3000;

% Keep fixedWin ONLY for defining pre-cue offset in per-unit spike window
% (fixedWin(2) not used for trial inclusion)
fixedWin = [-1000, 6000];   % ms relative to cue (fixedWin(2) not used for inclusion)

% Trial filters
MIN_RT_MS    = 100;
minTrialsPerUnit = 10;
requireValid = true;

% SDF settings
dt_ms        = 10;

% ============================================================
% MATCH FIG3A SMOOTHING METHOD ONLY:
%   - Gaussian kernel, unit area, sigma = 25 ms (same as FIG3A)
% ===========================================================
gaussSigmaMs = 50;   % (left unchanged per your original unless you set to 25)
% ============================================================

% Post-warp smoothing to reduce interpolation edge sharpening (set true/false)
doPostWarpSmooth = true;

% Display
xUnit = "s";
forceIntegerSecondsTicks = true;

% ============================================================
% PLOT STYLE (match FIG3A_STAGE_HEATMAPS)
% ============================================================
% Previous was [120,120,650,900]; make height 1/3 -> 300
figPos = [120, 120, 650, 300];

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

axesTickLW    = 4.0;
axesTickLen   = [0.02 0.02];

traceLW = 4;          % keep your trace thickness
eventLineLW = 5.0;    % match FIG3A event lines

% Event colors (match FIG3A palette)
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

% Force z-score display range on y-axis
yLimFixed = [-1 1];

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
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u} : 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skip)

%% ---- PRECOMPUTE SDF KERNEL (unit area) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (shared across BOTH stage plots) ----
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

        lr = lick  - cue;

        all_pressLat_global(end+1,1) = rt; %#ok<AGROW>
        all_lickLat_global(end+1,1)  = lr; %#ok<AGROW>
    end
end

assert(~isempty(all_pressLat_global) && ~isempty(all_lickLat_global), ...
    'No trials passed global filters (cannot set shared warp targets).');

Tpress_global = median(all_pressLat_global);
Tlick_global  = median(all_lickLat_global);

fprintf('\n=== GLOBAL targets (shared): Tpress=%.1f ms, Tlick=%.1f ms ===\n', ...
    Tpress_global, Tlick_global);

%% ---- STAGE DEFINITIONS (match your FIG3A heatmaps) ----
stageRanges = [1 8; 29 36];

for st = 1:size(stageRanges,1)
    s0 = stageRanges(st,1);
    s1 = stageRanges(st,2);

    % clamp to available sessions
    s0 = max(1, min(nSessions, s0));
    s1 = max(1, min(nSessions, s1));
    if s1 < s0
        fprintf('Stage %d range invalid after clamp; skipping.\n', st);
        continue;
    end

    sessIdxList = s0:s1;
    fprintf('\n=== STAGE %d: sessions %d-%d ===\n', st, s0, s1);

    %% 1) USE GLOBAL TARGETS
    Tpress = Tpress_global;
    Tlick  = Tlick_global;

    %% 2) TEMPLATE AXIS (same as FIG3A)
    winLeft  = -preCueMs;
    winRight_template = Tlick + postLickMs;

    tgrid_ms = winLeft:dt_ms:winRight_template;
    nT = numel(tgrid_ms);

    idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
    idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
    idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

    %% 3) COLLECT PER-UNIT MEAN(Z_trial)(t) TRACES
    Zmat = [];   % nT x nUnitsIncluded
    nUnits = 0;

    for sIdx = sessIdxList
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        % ---- Behavior-only keepTrial: valid, RT>=MIN_RT (NO fixedWin press/lick requirement) ----
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

            lr = lick  - cue;

            keepTrial(k) = true;
            cueAbs(k)    = cue;
            pressLat(k)  = rt;  % RT
            lickLat(k)   = lr;  % lick latency
        end

        idxKeep = find(keepTrial);
        if isempty(idxKeep), continue; end

        cueAbsK   = cueAbs(idxKeep);
        pressLatK = pressLat(idxKeep);
        lickLatK  = lickLat(idxKeep);

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

                % PER-UNIT activeTrials: spike in [cue + fixedWin(1), cue + lickLat]
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

                zTrials = zeros(nTrU, nT);

                for iTr = 1:nTrU
                    cue0   = cueU(iTr);
                    tPress = rtU(iTr);
                    tLick  = lkU(iTr);

                    winRight_trial = tLick + postLickMs;
                    nT_trial = numel(winLeft:dt_ms:winRight_trial);

                    idx0 = timeToIdx_(0,      winLeft, dt_ms, nT_trial);
                    idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                    idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

                    % SDF (Hz)
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

                    % Warp cue->press, press->lick; keep pre-cue and post-lick unwarped
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
                        segC_fit = [segC, zeros(1, kC_target - numel(segC))];
                    else
                        segC_fit = segC(1:kC_target);
                    end

                    yy = [pre, wA(2:end), wB(2:end), segC_fit(2:end)];

                    if numel(yy) < nT
                        yy = [yy, zeros(1, nT-numel(yy))];
                    elseif numel(yy) > nT
                        yy = yy(1:nT);
                    end

                    if doPostWarpSmooth
                        yy = conv(yy, g, 'same');
                    end

                    zTrials(iTr,:) = zscoreTrial_(yy);
                end

                muZ = mean(zTrials, 1);

                nUnits = nUnits + 1;
                Zmat(:, nUnits) = muZ(:);
            end
        end
    end

    if nUnits == 0
        fprintf('No units in stage %d; skipping.\n', st);
        continue;
    end

    %% 4) POPULATION TRACE + SEM
    zMean = mean(Zmat, 2);
    zSEM  = std(Zmat, 0, 2) / sqrt(nUnits);

    %% 5) PLOT
    if xUnit == "s"
        xPlot = tgrid_ms / 1000;
        xlab  = 'Warped time from cue (s)';
        cueLine   = 0;
        pressLine = Tpress/1000;
        lickLine  = Tlick/1000;
    else
        xPlot = tgrid_ms;
        xlab  = 'Warped time from cue (ms)';
        cueLine   = 0;
        pressLine = Tpress;
        lickLine  = Tlick;
    end

    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    patch(ax, [xPlot, fliplr(xPlot)], ...
        [(zMean - zSEM)', fliplr((zMean + zSEM)')], ...
        [0 0 0], 'FaceAlpha', 0.15, 'EdgeColor','none');

    plot(ax, xPlot, zMean, 'k-', 'LineWidth', traceLW);

    xline(ax, cueLine,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
    xline(ax, pressLine, '--', 'Color', colPress, 'LineWidth', eventLineLW);
    xline(ax, lickLine,  '--', 'Color', colLick,  'LineWidth', eventLineLW);

    xlabel(ax, xlab, 'FontSize', labelFontSize);
    ylabel(ax, 'Z-score', 'FontSize', labelFontSize);
    title(ax, sprintf('Sessions %d-%d: nUnits=%d', s0, s1, nUnits), ...
        'FontWeight','bold', 'FontSize', titleFontSize);

    set(ax, 'FontSize', tickFontSize, ...
        'TickDir','out', ...
        'LineWidth', axesTickLW, ...
        'TickLength', axesTickLen, ...
        'Box','off', ...
        'Layer','top', ...
        'XAxisLocation','bottom', ...
        'YAxisLocation','left');

    % Force y-axis limits
    ylim(ax, yLimFixed);

    % Make x tick labels "straight" (no rotation)
    xtickangle(ax, 0);

    if xUnit == "s" && forceIntegerSecondsTicks
        xl = xlim(ax);
        ts = ceil(xl(1)); te = floor(xl(2));
        if te >= ts, xticks(ax, ts:1:te); end
    end

    outpng = fullfile(outDir, sprintf('Stage%02d_sessions%02d_%02d_GLOBALWARP_WHOLETRIALZ.png', st, s0, s1));
    saveas(fig, outpng);
    close(fig);

    fprintf('Saved: %s\n', outpng);
end

fprintf('\nDONE. Saved in: %s\n', outDir);

%% ================= HELPERS =================
function z = zscoreTrial_(x)
% Whole-trial z-score across timepoints
    x = x(:)';
    mu = mean(x);
    sd = std(x);
    if ~isfinite(sd) || sd <= 0, sd = 1; end
    z = (x - mu) / sd;
end

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
% Gaussian kernel with sum=1 (unit area in bins)
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
% Floor/ceil weighted interpolation in INDEX space
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