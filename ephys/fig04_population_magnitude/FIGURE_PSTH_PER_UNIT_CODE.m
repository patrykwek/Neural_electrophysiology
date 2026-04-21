%% ========= FIG2_like_raster_and_warpedPSTH_ALLUNITS_HZ__FINAL_WITH_TRIALS__CELLTYPE.m =========
% Same analysis/plotting logic as your FIG2-like script, PLUS:
%   - Loads cell_type arranged session->channel->unit
%   - Appends cell-type label (MSN/FSI/TAN) to plot title and filename
%   - ALSO saves every figure as SVG (in addition to PNG)
%
% cell_type coding:
%   1 = MSN
%   2 = FSI
%   3 = TAN
%
% Behavior+spikes dataset:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_FINAL_WITH_TRIALS.mat
% Cell-type dataset:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/celltype_all_sessions_full.mat
% Output figures:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/2_FIGURE_PSTH

clear; clc;

%% ---- USER SETTINGS ----
matFile      = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
outDir       = '/Volumes/WD_BLACK/A_THESIS_FINAL/2_FIGURE_PSTH';

MIN_RT_MS = 100;
fixedWin  = [-1000, 5000];   % ms (relative to cue), plotting + SDF grid

% ---- SDF / PSTH SETTINGS ----
dt_ms        = 10;       % 10 ms grid
gaussSigmaMs = 50;       % sigma=50ms
bootstrapN   = 500;

% IMPORTANT:
% If you z-score, negatives are expected by design.
doZscore      = false;        % Fig 2-like rate: keep false
baselineWinMs = [-1000, 0];   % used only if doZscore==true

% ---- DISPLAY ----
xUnit = "s";
forceIntegerSecondsTicks = true;

% ---- STYLE ----
spikeMarkerSize = 15;
tickHalfHeight  = 1;
tickLW          = 7.5;
psthLW          = 5;
ciAlpha         = 0.20;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

colSpk   = [0 0 0];
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

% CI + event line styling
colCI       = [0.3 0.3 0.3];   % grey CI fill
eventLineLW = 5.0;             % thicker event xlines
axesTickLW  = 4.0;             % more visible ticks/axes
axesTickLen = [0.04 0.04];     % longer ticks

% ---- FIGURE SIZE ----
figPos = [120, 120, 680, 880];

if ~exist(outDir,'dir'), mkdir(outDir); end

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile); %#ok<NASGU>
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found in MAT file.');
nSessions_beh = numel(beh);

%% ---- LOAD CELL TYPE ----
assert(exist(cellTypeFile,'file')==2, 'Cell-type file not found: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type'), 'Cell-type MAT file lacks variable "cell_type".');
cell_type = CT.cell_type;
assert(iscell(cell_type), '"cell_type" must be a cell array indexed by session.');

nSessions_ct = numel(cell_type);
if nSessions_ct ~= nSessions_beh
    warning('Session count mismatch: beh=%d, cell_type=%d. Using min().', nSessions_beh, nSessions_ct);
end
nSessions = min(nSessions_beh, nSessions_ct);

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms); % area=1 in "bins"

%% ---- MAIN ----
plot_count = 0;

for sIdx = 1:nSessions
    session = beh(sIdx);

    if ~isGoodTrialsStruct_(session)
        fprintf('[sess %d] Missing/invalid session.trials; skipping.\n', sIdx);
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        fprintf('[sess %d] No spikes; skipping.\n', sIdx);
        continue;
    end

    sessDir = fullfile(outDir, sprintf('Session_%02d', sIdx));
    if ~exist(sessDir,'dir'), mkdir(sessDir); end

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)
            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
            spk_abs = double(spk_abs(:));

            % ---- Cell-type label (MSN/FSI/TAN/UNK) ----
            [ct_code, ct_label] = getCellTypeLabel_(cell_type, sIdx, ch, u); %#ok<NASGU>

            % ============================================================
            % TRIAL SELECTION (PER-UNIT, spikes within plotted window)
            % ============================================================
            keepMask = false(1, numel(trials));

            for k = 1:numel(trials)
                tr = trials(k);

                if ~isfield(tr,'valid') || ~tr.valid, continue; end
                if ~all(isfield(tr, {'cue','press','lick'})), continue; end

                cue   = double(tr.cue);
                press = double(tr.press);
                lick  = double(tr.lick);

                if any(~isfinite([cue, press, lick])), continue; end

                rt = press - cue;
                if rt < MIN_RT_MS, continue; end

                pr = press - cue;
                lr = lick  - cue;
                if ~(pr >= fixedWin(1) && pr <= fixedWin(2)), continue; end
                if ~(lr >= fixedWin(1) && lr <= fixedWin(2)), continue; end

                % per-unit inclusion requires a spike in the ACTUAL plotted window:
                w0 = cue + fixedWin(1);
                w1 = cue + fixedWin(2);
                if ~any(spk_abs >= w0 & spk_abs <= w1), continue; end

                keepMask(k) = true;
            end

            keepIdx = find(keepMask);
            if isempty(keepIdx), continue; end
            % ============================================================

            % ---- Extract ABS event times for kept trials ----
            cue_abs   = nan(numel(keepIdx),1);
            press_abs = nan(numel(keepIdx),1);
            lick_abs  = nan(numel(keepIdx),1);

            for iTr = 1:numel(keepIdx)
                k = keepIdx(iTr);
                cue_abs(iTr)   = double(trials(k).cue);
                press_abs(iTr) = double(trials(k).press);
                lick_abs(iTr)  = double(trials(k).lick);
            end

            rt_ms      = press_abs - cue_abs;
            licklat_ms = lick_abs  - cue_abs;

            % ---- Sort by RT ascending ----
            [rt_ms, ord] = sort(rt_ms, 'ascend');
            cue_abs      = cue_abs(ord);
            press_abs    = press_abs(ord);
            lick_abs     = lick_abs(ord);
            licklat_ms   = licklat_ms(ord);

            % ---- Common grid ----
            tgrid_ms = fixedWin(1):dt_ms:fixedWin(2);
            nT  = numel(tgrid_ms);
            nTr = numel(cue_abs);

            % ---- Build per-trial firing-rate SDF in Hz ----
            SDF_Hz = zeros(nTr, nT);

            for iTr = 1:nTr
                cue0 = cue_abs(iTr);

                spk_rel = spk_abs - cue0; % ms relative to cue
                spk_rel = spk_rel(spk_rel >= fixedWin(1) & spk_rel <= fixedWin(2));

                counts = zeros(1, nT); % spikes per dt_ms bin
                if ~isempty(spk_rel)
                    idx = round((spk_rel - fixedWin(1))/dt_ms) + 1;
                    idx = idx(idx >= 1 & idx <= nT);
                    for jj = 1:numel(idx)
                        counts(idx(jj)) = counts(idx(jj)) + 1;
                    end
                end

                smoothed_counts = conv(counts, g, 'same'); % spikes/bin (bin=dt_ms)
                hz = smoothed_counts / (dt_ms/1000);       % spikes/sec = Hz
                SDF_Hz(iTr,:) = hz;
            end

            % ---- Targets (medians) ----
            T2 = median(rt_ms);
            T3 = median(licklat_ms);

            % Convert target event times to indices on the common grid
            idx0  = timeToIdx_(0,  fixedWin(1), dt_ms, nT);
            idxT2 = timeToIdx_(T2, fixedWin(1), dt_ms, nT);
            idxT3 = timeToIdx_(T3, fixedWin(1), dt_ms, nT);

            % ---- Warp each trial ----
            warpedHz = zeros(nTr, nT);

            for iTr = 1:nTr
                t2 = rt_ms(iTr);
                t3 = licklat_ms(iTr);

                idx2 = timeToIdx_(t2, fixedWin(1), dt_ms, nT);
                idx3 = timeToIdx_(t3, fixedWin(1), dt_ms, nT);

                y = SDF_Hz(iTr,:);

                % Pre-cue: [start .. 0] unchanged
                pre = y(1:idx0);

                % Segment A: [0 .. press] -> [0 .. medianPress]
                segA = y(idx0:idx2);
                kA   = max(2, (idxT2 - idx0 + 1));
                wA   = warpSegment_floorceilEqn_(segA, kA);

                % Segment B: [press .. lick] -> [medianPress .. medianLick]
                segB = y(idx2:idx3);
                kB   = max(2, (idxT3 - idxT2 + 1));
                wB   = warpSegment_floorceilEqn_(segB, kB);

                % Segment C: [lick .. end] -> fill remaining length
                segC = y(idx3:nT);
                kC   = max(2, (nT - idxT3 + 1));
                wC   = warpSegment_floorceilEqn_(segC, kC);

                yy = [pre, wA(2:end), wB(2:end), wC(2:end)]; % drop overlap endpoints
                if numel(yy) < nT
                    yy = [yy, zeros(1, nT-numel(yy))];
                elseif numel(yy) > nT
                    yy = yy(1:nT);
                end

                warpedHz(iTr,:) = yy;
            end

            % ---- Optional z-score ----
            if doZscore
                baseIdx = (tgrid_ms >= baselineWinMs(1)) & (tgrid_ms <= baselineWinMs(2));
                if ~any(baseIdx), baseIdx = (tgrid_ms < 0); end

                baseVals = warpedHz(:, baseIdx);
                muBase = mean(baseVals(:));
                sdBase = std(baseVals(:));
                if sdBase <= 0, sdBase = 1; end

                warpedHz = (warpedHz - muBase) / sdBase;
            end

            % ---- Mean + bootstrap CI ----
            mu = mean(warpedHz, 1);

            bootMu = zeros(bootstrapN, nT);
            for b = 1:bootstrapN
                samp = randi(nTr, [nTr, 1]);
                bootMu(b,:) = mean(warpedHz(samp,:), 1);
            end
            lo = prctile(bootMu, 2.5, 1);
            hi = prctile(bootMu, 97.5, 1);

            % ---- Plot units ----
            if xUnit == "s"
                xPlot     = tgrid_ms / 1000;
                pressPlot = rt_ms / 1000;
                lickPlot  = licklat_ms / 1000;
                T2p       = T2 / 1000;
                T3p       = T3 / 1000;
                xlab      = 'Time from cue (s)';
                xlimPlot  = fixedWin / 1000;
            else
                xPlot     = tgrid_ms;
                pressPlot = rt_ms;
                lickPlot  = licklat_ms;
                T2p       = T2;
                T3p       = T3;
                xlab      = 'Time from cue (ms)';
                xlimPlot  = fixedWin;
            end

            % ============================================================
            % FIGURE: Raster + Warped PSTH
            % ============================================================
            fig = figure('Visible','off','Color','w','Position',figPos);
            tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

            % ----- RASTER -----
            ax1 = nexttile(tl, 1); hold(ax1,'on');

            for r = 1:nTr
                cue0 = cue_abs(r);
                spk_rel = spk_abs - cue0; % ms
                spk_rel = spk_rel(spk_rel >= fixedWin(1) & spk_rel <= fixedWin(2));
                if xUnit == "s", spk_rel = spk_rel / 1000; end

                % spikes first
                if ~isempty(spk_rel)
                    plot(ax1, spk_rel, r*ones(size(spk_rel)), '.', ...
                        'Color', colSpk, 'MarkerSize', spikeMarkerSize);
                end

                % event ticks after spikes (on top)
                hCue   = drawTick_(ax1, 0,            r, tickHalfHeight, tickLW, colCue);
                hPress = drawTick_(ax1, pressPlot(r), r, tickHalfHeight, tickLW, colPress);
                hLick  = drawTick_(ax1, lickPlot(r),  r, tickHalfHeight, tickLW, colLick);

                if ~isempty(hCue),   uistack(hCue,   'top'); end
                if ~isempty(hPress), uistack(hPress, 'top'); end
                if ~isempty(hLick),  uistack(hLick,  'top'); end
            end

            xlim(ax1, xlimPlot);
            ylim(ax1, [0, nTr+1]);
            set(ax1, 'YDir', 'reverse', 'FontSize', tickFontSize);
            ylabel(ax1, 'Trials', 'FontSize', labelFontSize);

            if xUnit == "s" && forceIntegerSecondsTicks
                xl = xlim(ax1);
                ts = ceil(xl(1)); te = floor(xl(2));
                if te >= ts, xticks(ax1, ts:1:te); end
            end

            % remove top x-axis numbers but keep ticks
            xticklabels(ax1, []);

            set(ax1, 'TickDir','out', 'LineWidth', axesTickLW, ...
                'TickLength', axesTickLen, 'Layer','top');

            title(ax1, sprintf('Session %d, Ch %d, Unit %d (%s) (n=%d trials)', ...
                sIdx, ch, u, ct_label, nTr), 'FontWeight','bold', 'FontSize', titleFontSize);

            % ----- PSTH -----
            ax2 = nexttile(tl, 2); hold(ax2,'on');

            patch(ax2, [xPlot, fliplr(xPlot)], [lo, fliplr(hi)], colCI, ...
                'FaceAlpha', ciAlpha, 'EdgeColor','none');

            plot(ax2, xPlot, mu, '-', 'Color', colSpk, 'LineWidth', psthLW);

            xline(ax2, 0,   '--', 'Color', colCue,   'LineWidth', eventLineLW);
            xline(ax2, T2p, '--', 'Color', colPress, 'LineWidth', eventLineLW);
            xline(ax2, T3p, '--', 'Color', colLick,  'LineWidth', eventLineLW);

            xlim(ax2, xlimPlot);
            set(ax2, 'FontSize', tickFontSize);

            if xUnit == "s" && forceIntegerSecondsTicks
                xl = xlim(ax2);
                ts = ceil(xl(1)); te = floor(xl(2));
                if te >= ts, xticks(ax2, ts:1:te); end
            end

            xlabel(ax2, xlab, 'FontSize', labelFontSize);
            if doZscore
                ylabel(ax2, 'Warped firing (z)', 'FontSize', labelFontSize);
            else
                ylabel(ax2, 'Firing rate (Hz)', 'FontSize', labelFontSize);
            end

            set(ax2, 'TickDir','out', 'LineWidth', axesTickLW, ...
                'TickLength', axesTickLen, 'Layer','top');

            % ---- Filename includes cell type ----
            ct_label_safe = regexprep(ct_label, '\s+', '');
            outpng = fullfile(sessDir, sprintf('sess%02d_ch%02d_unit%02d_%s_fig2like_HZ.png', ...
                sIdx, ch, u, ct_label_safe));

            saveas(fig, outpng);

            % ---- ALSO SAVE SVG (ONLY CHANGE REQUESTED) ----
            outsvg = strrep(outpng, '.png', '.svg');
            try
                print(fig, outsvg, '-dsvg');
            catch ME
                warning('SVG save failed for %s: %s', outsvg, ME.message);
            end
            % -----------------------------------------------

            close(fig);

            plot_count = plot_count + 1;
        end
    end
end

fprintf('\nDONE.\nSaved %d figures in: %s\n', plot_count, outDir);

%% ================= HELPERS =================

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
% Discrete Gaussian kernel with SUM = 1 (unit area in bins).
% Then smoothed_counts is spikes/bin, and dividing by dt (sec) yields Hz.
    halfWidth = ceil(5*sigmaMs/dtMs);
    x = (-halfWidth:halfWidth) * dtMs;
    g = exp(-0.5*(x./sigmaMs).^2);
    g = g / sum(g);
end

function idx = timeToIdx_(t_ms, winLeft_ms, dt_ms, nT)
% Convert time (ms) relative to cue into 1-based index on grid winLeft:dt:winRight
    idx = round((t_ms - winLeft_ms)/dt_ms) + 1;
    idx = max(1, min(nT, idx));
end

function ywarp = warpSegment_floorceilEqn_(y, k)
% Paper-style floor/ceil weighted interpolation in INDEX space.
% y: 1xN, k: desired length
    y = y(:)';
    n = numel(y);

    if n < 2 || k < 2
        ywarp = zeros(1, k);
        if n >= 1 && k >= 1, ywarp(1) = y(1); end
        return;
    end

    % index-domain scaling: map [1..n] -> [1..k]
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

function h = drawTick_(ax, x, r, halfHeight, lw, col)
% Return handle so we can uistack() and guarantee top rendering.
    h = [];
    if ~isfinite(x), return; end
    h = plot(ax, [x x], [r-halfHeight r+halfHeight], '-', ...
        'Color', col, 'LineWidth', lw);
end

function beh = pickBehStruct_(S)
% Prefer exact field names; fallback to first struct found.
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k = 1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k}))
            beh = S.(cands{k});
            return;
        end
    end
    f = fieldnames(S);
    for i=1:numel(f)
        v = S.(f{i});
        if isstruct(v) && numel(v)>1
            beh = v;
            return;
        end
    end
    for i=1:numel(f)
        v = S.(f{i});
        if isstruct(v)
            beh = v;
            return;
        end
    end
end

function tf = isGoodTrialsStruct_(session)
% Require session.trials struct with required fields.
    tf = false;
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    if ~isstruct(session.trials), return; end
    req = {'valid','t_start','t_end','cue','press','lick'};
    tf = all(isfield(session.trials, req));
end

function spikes = getSpikesForSession_(session, S, sessionIdx)
% Returns spikes_by_ch as cell array: spikes_by_ch{ch}{unit} = spikeTimesAbs
    spikes = [];

    if isfield(session,'spikes') && ~isempty(session.spikes)
        spikes = session.spikes;
        return
    end
    if isfield(S,'spikes_session') && ~isempty(S.spikes_session) && ...
            numel(S.spikes_session) >= sessionIdx && ~isempty(S.spikes_session{sessionIdx})
        spikes = S.spikes_session{sessionIdx};
        return
    end
    if isfield(S,'spikes_persession') && ~isempty(S.spikes_persession) && ...
            numel(S.spikes_persession) >= sessionIdx && ~isempty(S.spikes_persession{sessionIdx})
        spikes = S.spikes_persession{sessionIdx};
        return
    end
    if isfield(S,'spikes') && iscell(S.spikes) && ~isempty(S.spikes)
        spikes = S.spikes; % last resort
    end
end

function [ct_code, ct_label] = getCellTypeLabel_(cell_type, sIdx, ch, u)
% Robustly fetch cell type code and label from cell_type{sess}{ch}{unit}.
    ct_code  = NaN;
    ct_label = 'UNK';

    if isempty(cell_type) || ~iscell(cell_type), return; end
    if sIdx < 1 || sIdx > numel(cell_type), return; end

    sCell = cell_type{sIdx};
    if isempty(sCell) || ~iscell(sCell), return; end
    if ch < 1 || ch > numel(sCell), return; end

    chCell = sCell{ch};
    if isempty(chCell), return; end

    if iscell(chCell)
        if u < 1 || u > numel(chCell), return; end
        uVal = chCell{u};
    else
        try
            uVal = chCell(u);
        catch
            return;
        end
    end

    if isempty(uVal), return; end

    if isnumeric(uVal)
        ct_code = double(uVal(1));
    elseif iscell(uVal) && ~isempty(uVal) && isnumeric(uVal{1})
        ct_code = double(uVal{1}(1));
    else
        return;
    end

    if ct_code == 1
        ct_label = 'MSN';
    elseif ct_code == 2
        ct_label = 'FSI';
    elseif ct_code == 3
        ct_label = 'TAN';
    else
        ct_label = sprintf('UNK%d', round(ct_code));
    end
end