%% ========= FIG3_GLM_EVENT_KERNELS__TRIALSTART_END_WINDOWS__SESSIONLEVEL__SIMPLIFIED.m =========
% MINIMAL CHANGE (per your request):
% - Fit kernels using ALL valid trials (same as before)
% - Compute baseline level ONLY from a stricter subset of trials with a clean pre-cue baseline,
%   to reduce contamination from previous-trial reward activity.
% - Then subtract that baseline level from cue/press/lick kernels (per unit).
%
% How baseline-clean trials are enforced (simple rule):
%   Require: (cue time) - (previous trial lick time) >= MIN_GAP_FROM_PREV_LICK_MS
%   for trials k>1. Trial 1 is allowed by default (or you can exclude it).
%
% NOTE: This does NOT change which trials go into the GLM fit. It only changes
%       how baseline is estimated for baseline subtraction.

clear; clc;

%% ---------------- USER SETTINGS ----------------
matFile      = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL';
outDir  = fullfile(baseOut, 'GLM');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Trial filters
MIN_RT_MS    = 100;
requireValid = true;

% Unit inclusion
minTrialsPerUnit        = 1;
requireSpikeInTrialWin  = true;

% Binning / smoothing
dt_ms        = 10;
gaussSigmaMs = 25; %#ok<NASGU>  % kept for provenance only (not used)

% --- Pre-cue baseline ---
BASELINE_PRE_CUE_MS = 500;  % use [-500, 0) ms relative to cue as baseline; required to exist within trial

% --- NEW (baseline-only) "clean gap" rule ---
% Simple countermeasure for contamination by prior reward/lick:
MIN_GAP_FROM_PREV_LICK_MS = 1000;   % (edit) require >= 1 s since previous trial lick to use baseline

% Kernel lag windows
lags_cue_ms   = 0:dt_ms:750;
lags_press_ms = -1500:dt_ms:750;
lags_lick_ms  = -750:dt_ms:1500;

%% ---------------- BASIS SETTINGS ----------------
nBasisCue   = 6;
nBasisPress = 6;
nBasisLick  = 6;

%% ---------------- PROGRESS SETTINGS ----------------
printEverySessions = 1;
printEveryUnits    = 50;

tGlobal = tic;
fprintf('\n=== GLM START (MAX-SIMPLIFIED; SESSION DESIGN REUSE) ===\n');
fprintf('Output dir: %s\n', outDir);

%% ---------------- LOAD ----------------
assert(exist(matFile,'file')==2, 'Missing: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

assert(exist(cellTypeFile,'file')==2, 'Missing: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type.');
cell_type = CT.cell_type;

fprintf('Loaded %d sessions.\n', nSessions);

%% ---------------- Convert lags to bin offsets ----------------
lags_cue_bins   = round(lags_cue_ms   / dt_ms);
lags_press_bins = round(lags_press_ms / dt_ms);
lags_lick_bins  = round(lags_lick_ms  / dt_ms);

%% ---------------- PRECOMPUTE BASIS -> KERNEL MAPS ----------------
[Bcue_lag,   cueLag0Idx]   = make_cosine_basis_for_lags_(lags_cue_bins,   nBasisCue);
[Bpress_lag, pressLag0Idx] = make_cosine_basis_for_lags_(lags_press_bins, nBasisPress);
[Blick_lag,  lickLag0Idx]  = make_cosine_basis_for_lags_(lags_lick_bins,  nBasisLick);

% Precompute indices for pre-cue baseline region on cue lag grid
idxBaselineCue = (lags_cue_ms >= -BASELINE_PRE_CUE_MS) & (lags_cue_ms < 0);

%% ============================================================
% Fit GLM per unit; summarize to session-level (avoid pseudorep)
%% ============================================================
sessKernels_cue   = nan(nSessions, numel(lags_cue_ms));
sessKernels_press = nan(nSessions, numel(lags_press_ms));
sessKernels_lick  = nan(nSessions, numel(lags_lick_ms));

sessR2_full     = nan(nSessions,1);
sess_nUnitsUsed = zeros(nSessions,1);

% session-level mean pre-cue baseline level (in kernel units)
sessBaseline_preCue = nan(nSessions,1);

% Progress counters
nSessionsAttempted = 0;
nSessionsSkipped   = 0;
nSessionsFitOK     = 0;

nUnitsCandidateTotal = 0;
nUnitsUsedTotal      = 0;

for sIdx = 1:nSessions
    tSess = tic;
    nSessionsAttempted = nSessionsAttempted + 1;

    if mod(sIdx, printEverySessions) == 0 || sIdx == 1
        fprintf('\n--- Session %d/%d (%.1f%%) ---\n', sIdx, nSessions, 100*sIdx/nSessions);
    end

    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session)
        nSessionsSkipped = nSessionsSkipped + 1;
        fprintf('  [skip] bad/missing trials struct\n');
        continue;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        nSessionsSkipped = nSessionsSkipped + 1;
        fprintf('  [skip] no spikes_by_ch\n');
        continue;
    end

    % ---- Select trials using t_start/t_end windows (GLM TRIALS; unchanged) ----
    [idxKeep, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK] = ...
        selectTrials_TrialBounds_(trials, requireValid, MIN_RT_MS, BASELINE_PRE_CUE_MS);

    if isempty(idxKeep)
        nSessionsSkipped = nSessionsSkipped + 1;
        fprintf('  [skip] no valid trials after filters\n');
        continue;
    end

    % ---- NEW: determine which of the KEPT trials are "baseline-clean" ----
    % Only used for baseline estimation (not for GLM fit)
    isBaselineClean = baseline_clean_mask_from_prev_lick_(trials, idxKeep, MIN_GAP_FROM_PREV_LICK_MS);

    % ---- Precompute session timebase + design matrix ONCE ----
    nTr = numel(idxKeep);
    trialLensBins = zeros(nTr,1);
    for tr = 1:nTr
        trialLensBins(tr) = numel(tStartK(tr):dt_ms:tEndK(tr));
    end
    totalRows = sum(trialLensBins);
    trialOffsets = [0; cumsum(trialLensBins(1:end-1))];  % 0-based offsets

    % Build sparse basis-stamped design for all kept trials
    [XcueB, XpressB, XlickB] = build_session_design_basis_( ...
        cueAbsK, pressAbsK, lickAbsK, tStartK, dt_ms, trialLensBins, trialOffsets, ...
        Bcue_lag, cueLag0Idx, Bpress_lag, pressLag0Idx, Blick_lag, lickLag0Idx);

    % Assemble X once and cache QR decomposition
    X = [ones(totalRows,1), XcueB, XpressB, XlickB];
    D = decomposition(X, 'qr');  % reuse for all units in session

    % Indices into beta
    idxIntercept = 1;
    idxCue   = (idxIntercept+1) : (idxIntercept + size(XcueB,2));
    idxPress = (idxCue(end)+1)  : (idxCue(end) + size(XpressB,2));
    idxLick  = (idxPress(end)+1): (idxPress(end)+ size(XlickB,2));

    % Accumulate per-unit outputs (then average to session)
    Kcue_units_cell   = {};
    Kpress_units_cell = {};
    Klick_units_cell  = {};
    R2_units_cell     = {};
    Baseline_units_cell = {}; % per-unit pre-cue baseline level from cue-kernel (baseline-clean trials only)

    sessCandidate = 0;
    sessUsed      = 0;
    sessSkippedCt = 0;
    sessSkippedSpike = 0;
    sessSkippedTrials = 0;

    unitCounterThisSession = 0;

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)

            ct = safe_get_celltype_(cell_type, sIdx, ch, u);
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                sessSkippedCt = sessSkippedCt + 1;
                continue;
            end

            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs)
                sessSkippedSpike = sessSkippedSpike + 1;
                continue;
            end
            spk_abs = double(spk_abs(:));
            if isempty(spk_abs)
                sessSkippedSpike = sessSkippedSpike + 1;
                continue;
            end

            sessCandidate = sessCandidate + 1;
            nUnitsCandidateTotal = nUnitsCandidateTotal + 1;
            unitCounterThisSession = unitCounterThisSession + 1;

            if mod(unitCounterThisSession, printEveryUnits) == 0
                fprintf('  ch %d/%d | candidate units so far: %d (used %d)\n', ...
                    ch, numel(spikes_by_ch), sessCandidate, sessUsed);
            end

            % ---- Enforce minTrialsPerUnit (but DO NOT subset trials in the fit) ----
            if requireSpikeInTrialWin
                nActive = count_trials_with_spike_(spk_abs, tStartK, tEndK);
                if nActive < minTrialsPerUnit
                    sessSkippedTrials = sessSkippedTrials + 1;
                    continue;
                end
            end

            % ---- Build y over ALL kept trials on the session timebase ----
            y = build_y_session_(spk_abs, tStartK, tEndK, dt_ms, trialLensBins, trialOffsets, totalRows);

            % ---- OLS via cached decomposition ----
            beta = D \ y;

            % ---- R^2 ----
            yhat = X * beta;
            R2_full = R2_(y, yhat);

            % ---- Basis weights ----
            w_cue   = beta(idxCue);
            w_press = beta(idxPress);
            w_lick  = beta(idxLick);

            % ---- Reconstruct full kernels on original lag grids ----
            k_cue   = (Bcue_lag   * w_cue(:)).';
            k_press = (Bpress_lag * w_press(:)).';
            k_lick  = (Blick_lag  * w_lick(:)).';

            % ---- BASELINE subtraction using ONLY baseline-clean trials (MINIMAL CHANGE) ----
            % If no baseline-clean trials exist, fallback to original behavior (use all kept trials implicitly).
            if any(idxBaselineCue)
                if any(isBaselineClean)
                    % compute baseline from mean pre-cue SPIKE RATE across baseline-clean trials
                    % by directly measuring y in the pre-cue window (NOT from the kernel)
                    baseLevel = estimate_baseline_from_y_precue_( ...
                        y, cueAbsK, tStartK, dt_ms, trialLensBins, trialOffsets, ...
                        BASELINE_PRE_CUE_MS, isBaselineClean);
                else
                    % fallback (previous behavior): baseline from cue-kernel pre-cue region
                    baseLevel = mean(k_cue(idxBaselineCue), 'omitnan');
                end
            else
                baseLevel = 0;
            end

            k_cue   = k_cue   - baseLevel;
            k_press = k_press - baseLevel;
            k_lick  = k_lick  - baseLevel;

            Kcue_units_cell{end+1,1}      = k_cue;      %#ok<AGROW>
            Kpress_units_cell{end+1,1}    = k_press;    %#ok<AGROW>
            Klick_units_cell{end+1,1}     = k_lick;     %#ok<AGROW>
            R2_units_cell{end+1,1}        = R2_full;    %#ok<AGROW>
            Baseline_units_cell{end+1,1}  = baseLevel;  %#ok<AGROW>

            sessUsed = sessUsed + 1;
            nUnitsUsedTotal = nUnitsUsedTotal + 1;
        end
    end

    if isempty(R2_units_cell)
        nSessionsSkipped = nSessionsSkipped + 1;
        fprintf('  [skip] no units survived filters in this session. (ct-skip=%d, spk-skip=%d, trial-skip=%d)\n', ...
            sessSkippedCt, sessSkippedSpike, sessSkippedTrials);
        continue;
    end

    Kcue_units   = cell2mat(Kcue_units_cell);
    Kpress_units = cell2mat(Kpress_units_cell);
    Klick_units  = cell2mat(Klick_units_cell);
    R2_units     = cell2mat(R2_units_cell);
    Base_units   = cell2mat(Baseline_units_cell);

    sessKernels_cue(sIdx,:)   = mean(Kcue_units,   1, 'omitnan');
    sessKernels_press(sIdx,:) = mean(Kpress_units, 1, 'omitnan');
    sessKernels_lick(sIdx,:)  = mean(Klick_units,  1, 'omitnan');

    sessR2_full(sIdx)         = mean(R2_units, 'omitnan');
    sess_nUnitsUsed(sIdx)     = size(Kcue_units,1);
    sessBaseline_preCue(sIdx) = mean(Base_units, 'omitnan');

    nSessionsFitOK = nSessionsFitOK + 1;

    fprintf('  [ok] trials kept=%d | baseline-clean=%d | candidate units=%d | used units=%d | mean R2=%.3f | elapsed %.1fs\n', ...
        numel(idxKeep), nnz(isBaselineClean), sessCandidate, sessUsed, sessR2_full(sIdx), toc(tSess));
end

fprintf('\nDONE fitting (MAX-SIMPLIFIED) GLMs.\n');
fprintf('Sessions fit OK: %d / %d (skipped=%d)\n', nSessionsFitOK, nSessions, nSessionsSkipped);
fprintf('Units: candidate=%d | used=%d\n', nUnitsCandidateTotal, nUnitsUsedTotal);
fprintf('Total elapsed: %.1f min\n\n', toc(tGlobal)/60);

%% ============================================================
% SAVE outputs needed for later window-validity plots
%% ============================================================
saveFile = fullfile(outDir, 'GLM_SESSIONLEVEL_OUTPUTS_for_window_validity_SIMPLIFIED.mat');

GLM_OUTPUTS = struct();
GLM_OUTPUTS.sessKernels_cue   = sessKernels_cue;
GLM_OUTPUTS.sessKernels_press = sessKernels_press;
GLM_OUTPUTS.sessKernels_lick  = sessKernels_lick;

GLM_OUTPUTS.lags_cue_ms   = lags_cue_ms;
GLM_OUTPUTS.lags_press_ms = lags_press_ms;
GLM_OUTPUTS.lags_lick_ms  = lags_lick_ms;

GLM_OUTPUTS.sessR2_full         = sessR2_full;
GLM_OUTPUTS.sess_nUnitsUsed     = sess_nUnitsUsed;
GLM_OUTPUTS.sessBaseline_preCue = sessBaseline_preCue;

% provenance
GLM_OUTPUTS.dt_ms        = dt_ms;
GLM_OUTPUTS.gaussSigmaMs = gaussSigmaMs;
GLM_OUTPUTS.requireValid = requireValid;
GLM_OUTPUTS.MIN_RT_MS    = MIN_RT_MS;
GLM_OUTPUTS.minTrialsPerUnit       = minTrialsPerUnit;
GLM_OUTPUTS.requireSpikeInTrialWin = requireSpikeInTrialWin;

% baseline provenance
GLM_OUTPUTS.BASELINE_PRE_CUE_MS = BASELINE_PRE_CUE_MS;
GLM_OUTPUTS.MIN_GAP_FROM_PREV_LICK_MS = MIN_GAP_FROM_PREV_LICK_MS;

% basis provenance
GLM_OUTPUTS.nBasisCue   = nBasisCue;
GLM_OUTPUTS.nBasisPress = nBasisPress;
GLM_OUTPUTS.nBasisLick  = nBasisLick;

save(saveFile, 'GLM_OUTPUTS', '-v7.3');
fprintf('Saved GLM outputs: %s\n', saveFile);

fprintf('\nALL DONE. Outputs in: %s\n', outDir);

%% =============================== HELPERS ===============================

function [idxKeep, cueAbsK, pressAbsK, lickAbsK, tStartK, tEndK] = selectTrials_TrialBounds_(trials, requireValid, MIN_RT_MS, BASELINE_PRE_CUE_MS)
    keepTrial = false(numel(trials),1);
    cueAbs    = nan(numel(trials),1);
    pressAbs  = nan(numel(trials),1);
    lickAbs   = nan(numel(trials),1);
    tStart    = nan(numel(trials),1);
    tEnd      = nan(numel(trials),1);

    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end

        req = {'cue','press','lick','t_start','t_end'};
        if ~all(isfield(tr, req)), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        ts    = double(tr.t_start);
        te    = double(tr.t_end);

        if any(~isfinite([cue press lick ts te])), continue; end
        if te <= ts, continue; end

        if ~(cue>=ts && cue<=te), continue; end
        if ~(press>=ts && press<=te), continue; end
        if ~(lick>=ts && lick<=te), continue; end

        % require enough pre-cue time for baseline window
        if cue - BASELINE_PRE_CUE_MS < ts, continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end

        keepTrial(k) = true;
        cueAbs(k)    = cue;
        pressAbs(k)  = press;
        lickAbs(k)   = lick;
        tStart(k)    = ts;
        tEnd(k)      = te;
    end

    idxKeep   = find(keepTrial);
    cueAbsK   = cueAbs(idxKeep);
    pressAbsK = pressAbs(idxKeep);
    lickAbsK  = lickAbs(idxKeep);
    tStartK   = tStart(idxKeep);
    tEndK     = tEnd(idxKeep);
end

function isClean = baseline_clean_mask_from_prev_lick_(trials, idxKeep, MIN_GAP_FROM_PREV_LICK_MS)
    % isClean has length = numel(idxKeep), referring to those kept trials.
    % Rule: cue(current) - lick(previous) >= MIN_GAP_FROM_PREV_LICK_MS (for k>1).
    % If k==1 (first trial of session), mark true.
    isClean = false(numel(idxKeep),1);

    for i = 1:numel(idxKeep)
        k = idxKeep(i);
        if k <= 1
            isClean(i) = true;
            continue;
        end

        trCur = trials(k);
        trPrev = trials(k-1);

        if ~isfield(trCur,'cue') || ~isfield(trPrev,'lick')
            isClean(i) = true;
            continue;
        end

        cueCur = double(trCur.cue);
        lickPrev = double(trPrev.lick);

        if ~isfinite(cueCur) || ~isfinite(lickPrev)
            isClean(i) = true;
            continue;
        end

        isClean(i) = (cueCur - lickPrev) >= MIN_GAP_FROM_PREV_LICK_MS;
    end
end

function baseLevel = estimate_baseline_from_y_precue_(y, cueAbsK, tStartK, dt_ms, trialLensBins, trialOffsets, BASELINE_PRE_CUE_MS, isBaselineClean)
    % Estimate baseline as mean firing rate (Hz) in [-BASELINE_PRE_CUE_MS, 0) relative to cue,
    % pooled across baseline-clean trials.
    %
    % Uses the already-built session y (binned spike rate), so it's fast and simple.

    y = y(:);
    nTr = numel(cueAbsK);

    vals = [];

    for tr = 1:nTr
        if ~isBaselineClean(tr), continue; end

        ts = tStartK(tr);
        nT = trialLensBins(tr);
        off = trialOffsets(tr);

        % pre-cue window absolute times
        t0 = cueAbsK(tr) - BASELINE_PRE_CUE_MS;
        t1 = cueAbsK(tr);

        % convert to bin indices within this trial (1..nT)
        i0 = round((t0 - ts)/dt_ms) + 1;
        i1 = round((t1 - ts)/dt_ms) + 1;

        i0 = max(1, min(nT, i0));
        i1 = max(1, min(nT, i1));

        if i1 <= i0
            continue;
        end

        yy = y(off + (i0:(i1-1))); % [-BASELINE, 0)
        vals = [vals; yy(:)]; %#ok<AGROW>
    end

    if isempty(vals)
        baseLevel = 0;
    else
        baseLevel = mean(vals, 'omitnan');
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

function nActive = count_trials_with_spike_(spk_abs, tStartK, tEndK)
    spk = spk_abs(:);
    if ~issorted(spk), spk = sort(spk); end

    nTr = numel(tStartK);
    active = false(nTr,1);

    j = 1; nSpk = numel(spk);
    for iTr = 1:nTr
        ts = tStartK(iTr); te = tEndK(iTr);
        while j <= nSpk && spk(j) < ts
            j = j + 1;
        end
        if j <= nSpk && spk(j) <= te
            active(iTr) = true;
        end
    end
    nActive = nnz(active);
end

function y = build_y_session_(spk_abs, tStartK, tEndK, dt_ms, trialLensBins, trialOffsets, totalRows)
    y = zeros(totalRows,1);

    spk = spk_abs(:);
    if isempty(spk), return; end
    if ~issorted(spk), spk = sort(spk); end

    nTr = numel(tStartK);
    for tr = 1:nTr
        ts = tStartK(tr);
        te = tEndK(tr);
        nT = trialLensBins(tr);
        off = trialOffsets(tr);

        s = spk(spk >= ts & spk <= te);
        if isempty(s), continue; end

        jj = round((s - ts)/dt_ms) + 1;
        jj = jj(jj>=1 & jj<=nT);
        if isempty(jj), continue; end

        counts = accumarray(jj, 1, [nT 1], @sum, 0);
        y(off + (1:nT)) = y(off + (1:nT)) + counts / (dt_ms/1000);
    end
end

function [XcueB, XpressB, XlickB] = build_session_design_basis_( ...
    cueAbsK, pressAbsK, lickAbsK, tStartK, dt_ms, trialLensBins, trialOffsets, ...
    Bcue_lag, cueLag0Idx, Bpress_lag, pressLag0Idx, Blick_lag, lickLag0Idx)

    nTr = numel(cueAbsK);
    totalRows = sum(trialLensBins);

    nBCue   = size(Bcue_lag,2);
    nBPress = size(Bpress_lag,2);
    nBLick  = size(Blick_lag,2);

    XcueB   = spalloc(totalRows, nBCue,   nTr * nBCue   * 200);
    XpressB = spalloc(totalRows, nBPress, nTr * nBPress * 200);
    XlickB  = spalloc(totalRows, nBLick,  nTr * nBLick  * 200);

    for tr = 1:nTr
        ts = tStartK(tr);
        nT = trialLensBins(tr);
        off = trialOffsets(tr);

        icue   = timeToIdx_abs_(cueAbsK(tr),   ts, dt_ms, nT);
        ipress = timeToIdx_abs_(pressAbsK(tr), ts, dt_ms, nT);
        ilick  = timeToIdx_abs_(lickAbsK(tr),  ts, dt_ms, nT);

        rr = off + (1:nT);

        XcueB(rr,:)   = XcueB(rr,:)   + stamp_basis_local_(nT, icue,   Bcue_lag,   cueLag0Idx);
        XpressB(rr,:) = XpressB(rr,:) + stamp_basis_local_(nT, ipress, Bpress_lag, pressLag0Idx);
        XlickB(rr,:)  = XlickB(rr,:)  + stamp_basis_local_(nT, ilick,  Blick_lag,  lickLag0Idx);
    end
end

function Xloc = stamp_basis_local_(nT, eventIdx, B_lag, lag0Idx)
    nLags = size(B_lag,1);
    nB = size(B_lag,2);

    lagRow = (1:nLags);
    timeIdx = eventIdx + (lagRow - lag0Idx);

    keep = (timeIdx >= 1) & (timeIdx <= nT);
    timeIdx = timeIdx(keep);
    Bkeep = B_lag(keep,:);

    [r,c,v] = find(Bkeep);
    rows = timeIdx(r);
    cols = c;
    vals = v;
    Xloc = sparse(rows, cols, vals, nT, nB);
end

function idx = timeToIdx_abs_(t_abs_ms, tStart_abs_ms, dt_ms, nT)
    idx = round((t_abs_ms - tStart_abs_ms)/dt_ms) + 1;
    idx = max(1, min(nT, idx));
end

function r2 = R2_(y, yhat)
    y = y(:); yhat = yhat(:);
    ssr = sum((y - yhat).^2);
    sst = sum((y - mean(y)).^2) + eps;
    r2 = 1 - ssr/sst;
end

function [B_lag, lag0Idx] = make_cosine_basis_for_lags_(lags_bins, nBasis)
    lags_bins = lags_bins(:);
    [~, lag0Idx] = min(abs(lags_bins));

    x = lags_bins;
    xmin = min(x); xmax = max(x);
    if xmax == xmin
        B_lag = ones(numel(x), nBasis);
        B_lag = B_lag ./ max(vecnorm(B_lag,2,1), eps);
        return;
    end

    xn = (x - xmin) / (xmax - xmin);
    centers = linspace(0, 1, nBasis);
    width = 1/(nBasis-1 + eps);

    B_lag = zeros(numel(x), nBasis);
    for b = 1:nBasis
        z = (xn - centers(b)) / (width + eps);
        Bb = zeros(numel(x),1);
        m = abs(z) <= 1;
        Bb(m) = 0.5*(1 + cos(pi*z(m)));
        B_lag(:,b) = Bb;
    end
    B_lag = B_lag ./ max(vecnorm(B_lag,2,1), eps);
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
    req = {'valid','cue','press','lick','t_start','t_end'};
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