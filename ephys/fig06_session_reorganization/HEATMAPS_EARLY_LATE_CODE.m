%% ========= FIG3A_STAGE_HEATMAPS__A_THESIS_FINAL__GLOBALWARP__FIG2TRIALS_EXACT__CLASSIFIEDONLY.m =========
% Makes FIG3A-style stage heatmaps (early vs late) for your CURRENT ratbehstruct.
%
% Key properties:
%   1) TRIAL SELECTION: valid, RT>=MIN_RT (NO fixedWin press/lick requirement)
%   2) PER-UNIT trial inclusion additionally requires >=1 spike in trial-defined window:
%        [cue + fixedWin(1), cue + lickLat]  (cue-1s to cue->lick)
%      (fixedWin(2) is NOT used for trial inclusion)
%   3) ONLY CLASSIFIED units are used: cell_type in {1,2,3} (MSN/FSI/TAN)
%   4) GLOBAL warp targets are computed ONCE from ALL sessions (behavior-only, no spike requirement)
%      and used for BOTH stage plots so timing is comparable.
%   5) For each unit: compute warped, per-trial z-scored SDF; use half-trial sort (odd) and
%      half-trial plot (even) to reduce circular sorting bias.
%
% CHANGE REQUEST IMPLEMENTED:
%   - Removed the -1000..6000 ms fixedWin requirement for press/lick trial inclusion.
%   - No other unnecessary changes.

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/FIGURES/4_1_HEATMAPS';
outDir  = fullfile(baseOut, '4_HEATMAPS');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Heatmap window definition (relative to cue) — display design
preCueMs   = 500;     % show -500 ms to cue
postLickMs = 3000;    % show +3000 ms after lick

% Keep fixedWin ONLY for defining pre-cue offset in per-unit spike window
% (press/lick inclusion in this script does NOT require fixedWin(2))
fixedWin = [-1000, 6000];   % ms relative to cue (fixedWin(2) not used for inclusion)

% Trial filters
MIN_RT_MS = 100;
minTrialsPerUnit = 10;
requireValid = true;

% SDF settings
dt_ms        = 10;

% ============================================================
% MATCH FIG3A_SESSION_ENSEMBLE_HEATMAP smoothing method ONLY:
%   - Gaussian kernel, unit area, sigma = 25 ms
%   - NO display-only smoothing (no conv2), NO bilinear interpolation
% ============================================================
gaussSigmaMs = 25;
% ============================================================

% FIXED DISPLAY RANGE (same as your script)
climFixed = [-2 2];

% Display
xUnit = "s";   % "s" or "ms"
forceIntegerSecondsTicks = true;

% ---- FIGURE SIZE (MATCH FIG3A_SESSION_ENSEMBLE_HEATMAP) ----
figPos = [120, 120, 650, 900];

% ---- COLORBAR formatting ----
cbFontSize   = 22;

% ---- STYLE ----
titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

eventLineLW   = 5.0;
axesTickLW    = 4.0;
axesTickLen   = [0.02 0.02];

% Behavioral event colors
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

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

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (shared across ALL stage plots) ----
% Behavior-only inclusion (NO per-unit spike requirement)
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

        % (Removed fixedWin press/lick inclusion requirement)

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

%% ---- STAGE DEFINITIONS (EARLY vs LATE) ----
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

    %% 1) USE GLOBAL WARP TARGETS (shared across stages)
    Tpress = Tpress_global;
    Tlick  = Tlick_global;
    fprintf('Stage %d using GLOBAL targets: Tpress=%.1f ms, Tlick=%.1f ms\n', st, Tpress, Tlick);

    %% 2) TEMPLATE AXIS FOR HEATMAP (display design)
    winLeft  = -preCueMs;
    winRight_template = Tlick + postLickMs;

    tgrid_ms = winLeft:dt_ms:winRight_template;
    nT = numel(tgrid_ms);

    idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
    idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
    idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

    %% 3) BUILD HEATMAP MATRIX (nT x nUnits) + sorting indices
    PETH = [];      % nT x nUnits (held-out mean)
    peakIdx = [];   % peak index from held-in half
    nUnits = 0;

    for sIdx = sessIdxList
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        % ---- Behavior-only keep trials (NO fixedWin press/lick requirement) ----
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

            % (Removed fixedWin press/lick inclusion requirement)

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

                    % trial-specific right edge (to include post-lick tail)
                    winRight_trial = tLick + postLickMs;
                    tgrid_trial = winLeft:dt_ms:winRight_trial;
                    nT_trial = numel(tgrid_trial);

                    idx0 = timeToIdx_(0,      winLeft, dt_ms, nT_trial);
                    idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                    idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

                    % SDF (Hz) on the trial window
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

                    % Warp only cue->press and press->lick; keep pre-cue and post-lick unwarped
                    pre = y(1:idx0);

                    segA = y(idx0:idxP);
                    kA   = max(2, (idxPT - idx0T + 1));
                    wA   = warpSegment_floorceilEqn_(segA, kA);

                    segB = y(idxP:idxL);
                    kB   = max(2, (idxLT - idxPT + 1));
                    wB   = warpSegment_floorceilEqn_(segB, kB);

                    segC = y(idxL:end);
                    kC_target = max(1, (nT - idxLT + 1)); % include idxLT..end
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

                    % per-trial z-score across time
                    warpedZ(iTr,:) = zscoreTrial_(yy);
                end

                % Half-trial sorting (odd) applied to held-out plot (even)
                idxSort = 1:2:nTrU;
                idxPlot = 2:2:nTrU;
                if isempty(idxPlot), idxPlot = idxSort; end
                if isempty(idxSort), idxSort = idxPlot; end

                muSort = mean(warpedZ(idxSort,:), 1);
                [~, pk] = max(muSort);

                muPlot = mean(warpedZ(idxPlot,:), 1);

                nUnits = nUnits + 1;
                PETH(:, nUnits) = muPlot(:);
                peakIdx(1, nUnits) = pk; %#ok<AGROW>
            end
        end
    end

    if nUnits == 0
        fprintf('No active units in stage %d; skipping.\n', st);
        continue;
    end

    fprintf('Stage %d: kept %d active CLASSIFIED units.\n', st, nUnits);

    %% 4) Sort by peak, then reverse order
    [~, ord] = sort(peakIdx, 'ascend');
    ord = ord(end:-1:1);
    P = PETH(:, ord);

    %% 5) X axis + event lines
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

    %% 6) Plot heatmap
    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    Pimg = P'; % units x time

    hImg = imagesc(ax, xPlot, 1:nUnits, Pimg); axis(ax,'tight');
    set(ax,'YDir','normal');

    try
        set(hImg, 'Interpolation', 'nearest');
    catch
        set(hImg, 'Interpolation', 'nearest');
    end

    colormap(ax, flipud(gray(256)));

    xlabel(ax, xlab, 'FontSize', labelFontSize);
    ylabel(ax, 'Unit', 'FontSize', labelFontSize);

    title(ax, sprintf('Sessions %d-%d: nUnits=%d | GLOBAL Tpress=%.0fms | GLOBAL Tlick=%.0fms', ...
        s0, s1, nUnits, Tpress, Tlick), 'FontWeight','bold', 'FontSize', titleFontSize);

    xline(ax, cueLine,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
    xline(ax, pressLine, '--', 'Color', colPress, 'LineWidth', eventLineLW);
    xline(ax, lickLine,  '--', 'Color', colLick,  'LineWidth', eventLineLW);

    caxis(ax, climFixed);

    set(ax, 'FontSize', tickFontSize, ...
        'TickDir','out', ...
        'LineWidth', axesTickLW, ...
        'TickLength', axesTickLen, ...
        'Box','off', ...
        'Layer','top', ...
        'XAxisLocation','bottom', ...
        'YAxisLocation','left');

    if xUnit == "s" && forceIntegerSecondsTicks
        xl = xlim(ax);
        ts = ceil(xl(1)); te = floor(xl(2));
        if te >= ts, xticks(ax, ts:1:te); end
    end

    %% ---- COLORBAR ----
    cb = colorbar(ax, 'eastoutside');
    cb.Label.String = 'Z score';
    cb.FontSize = cbFontSize;
    cb.Label.FontSize = cbFontSize;

    cb.LineWidth  = 2;
    cb.TickLength = 0.02;
    cb.Units = 'normalized';

    axpos = ax.Position;
    cbpos = cb.Position;

    cbpos(3) = 0.040;
    cbpos(4) = 0.50 * axpos(4);
    cbpos(2) = axpos(2) + (axpos(4) - cbpos(4))/2;
    cbpos(1) = axpos(1) + axpos(3) + 0.07*axpos(3);
    cbpos(1) = min(cbpos(1), 0.92);

    cb.Position = cbpos;

    %% ---- SAVE (PNG + SVG) ----
    outbase = fullfile(outDir, sprintf('Stage%02d_sessions%02d_%02d_heatmap_GLOBALWARP', st, s0, s1));
    outpng  = [outbase '.png'];
    outsvg  = [outbase '.svg'];

    saveas(fig, outpng);

    try
        saveas(fig, outsvg);
    catch
        print(fig, outsvg, '-dsvg');
    end
    close(fig);

    fprintf('Saved: %s\n', outpng);
    fprintf('Saved: %s\n', outsvg);
end

fprintf('\nDONE. Heatmaps saved in: %s\n', outDir);

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

function k = gaussian1dPix_(sigmaPix)
    if ~isfinite(sigmaPix) || sigmaPix <= 0
        k = 1; return;
    end
    hw = max(1, ceil(4*sigmaPix));
    x = -hw:hw;
    k = exp(-0.5*(x./sigmaPix).^2);
    k = k / sum(k);
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