%% ========= FIG__TRIAL_ACTIVITY_ZSCORE__EARLY_LATE__CORRECT_INCORRECT__GLOBALWARP.m =========
% Makes 4 ensemble plots of z-scored per-trial activity:
%   1) EARLY (sessions 1–8)   - CORRECT trials only
%   2) EARLY (sessions 1–8)   - INCORRECT trials only
%   3) LATE  (sessions 28–36) - CORRECT trials only
%   4) LATE  (sessions 28–36) - INCORRECT trials only
%
% Assumptions from user:
%   - beh(s).trials is a CELL ARRAY (one cell per trial)
%   - each trial cell contains a TABLE with variable/column named "correct"
%     (rows = trials; correct=1 correct, 0 incorrect)
%
% Pipeline (kept consistent with your prior global-warp scripts):
%   - Behavioral trial filters: valid (optional), RT>=MIN_RT_MS, press/lick inside fixedWin
%   - Per-unit trial inclusion: >=1 spike in [cue+fixedWin(1), cue+fixedWin(2)]
%   - Classified-only units: cell_type in {1,2,3}
%   - GLOBAL warp targets computed ONCE from ALL sessions (behavior-only)
%   - Per-trial SDF -> warp cue->press and press->lick -> zscore across time
%   - For each unit: held-out mean uses even trials (idxPlot = 2:2:nTrU, fallback all)
%   - Ensemble: average across units, then zscore across time
%
% Outputs:
%   - 4 PNGs saved in outDir
%   - .mat file with all four ensemble vectors and metadata

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL';
outDir  = fullfile(baseOut, 'CORRECT_INCORRECT_LATE_EARLY');
if ~exist(outDir,'dir'), mkdir(outDir); end

% Heatmap/ensemble window definition (relative to cue)
preCueMs   = 500;     % show -500 ms to cue
postLickMs = 3000;    % show +3000 ms after lick

% FIG2-like trial-selection window (relative to cue)
fixedWin = [-1000, 6000];   % ms relative to cue

% Trial filters
MIN_RT_MS = 100;
minTrialsPerUnit = 10;
requireValid = true;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Group session ranges
earlySessRange = [1 8];
lateSessRange  = [28 36];

% Plot
figPos = [120, 120, 1050, 720];
lw = 3.0;
fsTitle = 22;
fsAx = 18;

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

% clamp session ranges to available sessions
earlySessRange(1) = max(1, min(nSessions, earlySessRange(1)));
earlySessRange(2) = max(1, min(nSessions, earlySessRange(2)));
lateSessRange(1)  = max(1, min(nSessions, lateSessRange(1)));
lateSessRange(2)  = max(1, min(nSessions, lateSessRange(2)));

earlySessIdxList = earlySessRange(1):earlySessRange(2);
lateSessIdxList  = lateSessRange(1):lateSessRange(2);

%% ---- LOAD CELL TYPES ----
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS (shared across all groups) ----
all_pressLat_global = [];
all_lickLat_global  = [];

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStructFlexible_(session), continue; end
    trials = session.trials;

    nTr = numel(trials);
    for k = 1:nTr
        tr = getTrialStruct_(trials, k);
        if isempty(tr), continue; end

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

%% ---- TEMPLATE AXIS FOR ENSEMBLE ----
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

%% ---- COMPUTE 4 GROUP ENSEMBLES ----
fprintf('\n=== GROUPS ===\n');
fprintf('EARLY sessions: %d-%d\n', earlySessIdxList(1), earlySessIdxList(end));
fprintf('LATE  sessions: %d-%d\n', lateSessIdxList(1),  lateSessIdxList(end));

ens_E_correct   = computeGroupEnsemble_(beh, S, cell_type, earlySessIdxList, true, ...
    Tpress, Tlick, winLeft, postLickMs, fixedWin, MIN_RT_MS, minTrialsPerUnit, requireValid, dt_ms, g, ...
    nT, idx0T, idxPT, idxLT);

ens_E_incorrect = computeGroupEnsemble_(beh, S, cell_type, earlySessIdxList, false, ...
    Tpress, Tlick, winLeft, postLickMs, fixedWin, MIN_RT_MS, minTrialsPerUnit, requireValid, dt_ms, g, ...
    nT, idx0T, idxPT, idxLT);

ens_L_correct   = computeGroupEnsemble_(beh, S, cell_type, lateSessIdxList, true, ...
    Tpress, Tlick, winLeft, postLickMs, fixedWin, MIN_RT_MS, minTrialsPerUnit, requireValid, dt_ms, g, ...
    nT, idx0T, idxPT, idxLT);

ens_L_incorrect = computeGroupEnsemble_(beh, S, cell_type, lateSessIdxList, false, ...
    Tpress, Tlick, winLeft, postLickMs, fixedWin, MIN_RT_MS, minTrialsPerUnit, requireValid, dt_ms, g, ...
    nT, idx0T, idxPT, idxLT);

%% ---- X axis + event lines ----
xPlot = tgrid_ms / 1000;
cueLine   = 0;
pressLine = Tpress/1000;
lickLine  = Tlick/1000;

%% ---- PLOT 4 PANELS (keep figures open; do not close) ----
% 1) EARLY correct
fig1 = figure('Color','w','Position',figPos);
ax1 = axes(fig1); hold(ax1,'on');
plot(ax1, xPlot, ens_E_correct, 'LineWidth', lw);
xline(ax1, cueLine,   '--', 'Color', colCue,   'LineWidth', 2.5);
xline(ax1, pressLine, '--', 'Color', colPress, 'LineWidth', 2.5);
xline(ax1, lickLine,  '--', 'Color', colLick,  'LineWidth', 2.5);
xlabel(ax1, 'Warped time from Cue (s)'); ylabel(ax1, 'Z score');
title(ax1, sprintf('EARLY (sessions %d-%d): CORRECT trials', earlySessIdxList(1), earlySessIdxList(end)), ...
    'FontSize', fsTitle, 'FontWeight','bold');
set(ax1,'FontSize',fsAx,'TickDir','out','LineWidth',2); box(ax1,'off');
out1 = fullfile(outDir, sprintf('EARLY_sess%02d_%02d_CORRECT_ensembleZ.png', earlySessIdxList(1), earlySessIdxList(end)));
saveas(fig1, out1);

% 2) EARLY incorrect
fig2 = figure('Color','w','Position',figPos);
ax2 = axes(fig2); hold(ax2,'on');
plot(ax2, xPlot, ens_E_incorrect, 'LineWidth', lw);
xline(ax2, cueLine,   '--', 'Color', colCue,   'LineWidth', 2.5);
xline(ax2, pressLine, '--', 'Color', colPress, 'LineWidth', 2.5);
xline(ax2, lickLine,  '--', 'Color', colLick,  'LineWidth', 2.5);
xlabel(ax2, 'Warped time from Cue (s)'); ylabel(ax2, 'Z score');
title(ax2, sprintf('EARLY (sessions %d-%d): INCORRECT trials', earlySessIdxList(1), earlySessIdxList(end)), ...
    'FontSize', fsTitle, 'FontWeight','bold');
set(ax2,'FontSize',fsAx,'TickDir','out','LineWidth',2); box(ax2,'off');
out2 = fullfile(outDir, sprintf('EARLY_sess%02d_%02d_INCORRECT_ensembleZ.png', earlySessIdxList(1), earlySessIdxList(end)));
saveas(fig2, out2);

% 3) LATE correct
fig3 = figure('Color','w','Position',figPos);
ax3 = axes(fig3); hold(ax3,'on');
plot(ax3, xPlot, ens_L_correct, 'LineWidth', lw);
xline(ax3, cueLine,   '--', 'Color', colCue,   'LineWidth', 2.5);
xline(ax3, pressLine, '--', 'Color', colPress, 'LineWidth', 2.5);
xline(ax3, lickLine,  '--', 'Color', colLick,  'LineWidth', 2.5);
xlabel(ax3, 'Warped time from Cue (s)'); ylabel(ax3, 'Z score');
title(ax3, sprintf('LATE (sessions %d-%d): CORRECT trials', lateSessIdxList(1), lateSessIdxList(end)), ...
    'FontSize', fsTitle, 'FontWeight','bold');
set(ax3,'FontSize',fsAx,'TickDir','out','LineWidth',2); box(ax3,'off');
out3 = fullfile(outDir, sprintf('LATE_sess%02d_%02d_CORRECT_ensembleZ.png', lateSessIdxList(1), lateSessIdxList(end)));
saveas(fig3, out3);

% 4) LATE incorrect
fig4 = figure('Color','w','Position',figPos);
ax4 = axes(fig4); hold(ax4,'on');
plot(ax4, xPlot, ens_L_incorrect, 'LineWidth', lw);
xline(ax4, cueLine,   '--', 'Color', colCue,   'LineWidth', 2.5);
xline(ax4, pressLine, '--', 'Color', colPress, 'LineWidth', 2.5);
xline(ax4, lickLine,  '--', 'Color', colLick,  'LineWidth', 2.5);
xlabel(ax4, 'Warped time from Cue (s)'); ylabel(ax4, 'Z score');
title(ax4, sprintf('LATE (sessions %d-%d): INCORRECT trials', lateSessIdxList(1), lateSessIdxList(end)), ...
    'FontSize', fsTitle, 'FontWeight','bold');
set(ax4,'FontSize',fsAx,'TickDir','out','LineWidth',2); box(ax4,'off');
out4 = fullfile(outDir, sprintf('LATE_sess%02d_%02d_INCORRECT_ensembleZ.png', lateSessIdxList(1), lateSessIdxList(end)));
saveas(fig4, out4);

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', out1, out2, out3, out4);

% Save outputs to MAT
outmat = fullfile(outDir, 'GROUP_ENSEMBLES_EARLY_LATE_CORRECT_INCORRECT.mat');
save(outmat, 'ens_E_correct', 'ens_E_incorrect', 'ens_L_correct', 'ens_L_incorrect', ...
    'tgrid_ms', 'Tpress', 'Tlick', 'earlySessIdxList', 'lateSessIdxList', ...
    '-v7.3');
fprintf('Saved: %s\n', outmat);

fprintf('\nDONE. Figures are left OPEN and saved in: %s\n', outDir);

%% ================= LOCAL HELPERS =================

function ensZ = computeGroupEnsemble_(beh, S, cell_type, sessIdxList, wantCorrect, ...
    Tpress, Tlick, winLeft, postLickMs, fixedWin, MIN_RT_MS, minTrialsPerUnit, requireValid, dt_ms, g, ...
    nT, idx0T, idxPT, idxLT)

    % Returns 1 x nT ensemble z-scored curve for the specified group (sessions + correctness filter)

    allUnitCurves = []; % concat across sessions: nUnitsTotal x nT

    for sIdx = sessIdxList
        session = beh(sIdx);
        if ~isGoodTrialsStructFlexible_(session), continue; end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

        % ---- Determine which trials are correct/incorrect (per user format) ----
        isCorr = getCorrectVectorFromTrialsCell_(trials);
        if isempty(isCorr)
            % fallback: try struct field "correct" if trials are struct array
            isCorr = getCorrectVectorFromTrialsStruct_(trials);
        end
        if isempty(isCorr)
            continue;
        end

        % ---- Behavior-only keep trials with correctness constraint ----
        nTr = numel(trials);
        keepTrial = false(nTr,1);
        cueAbs   = nan(nTr,1);
        pressLat = nan(nTr,1);
        lickLat  = nan(nTr,1);

        for k = 1:nTr
            tr = getTrialStruct_(trials, k);
            if isempty(tr), continue; end

            % correctness filter
            if wantCorrect
                if ~(isCorr(k) == 1), continue; end
            else
                if ~(isCorr(k) == 0), continue; end
            end

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

        % ---- Per-unit curves for this session ----
        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end

            for u = 1:numel(uc)

                % classified-only
                ct = safe_get_celltype_(cell_type, sIdx, ch, u);
                if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                    continue;
                end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                % per-unit active trial filter: spike in [cue+fixedWin(1), cue+fixedWin(2)]
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

                    % trial-specific right edge
                    winRight_trial = tLick + postLickMs;
                    tgrid_trial = winLeft:dt_ms:winRight_trial;
                    nT_trial = numel(tgrid_trial);

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
                    y  = sm / (dt_ms/1000);  % Hz

                    % warp cue->press and press->lick; keep pre-cue and post-lick unwarped
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

                % held-out even trials mean
                idxPlot = 2:2:nTrU;
                if isempty(idxPlot), idxPlot = 1:nTrU; end
                muPlot = mean(warpedZ(idxPlot,:), 1);

                allUnitCurves(end+1,:) = muPlot; %#ok<AGROW>
            end
        end
    end

    if isempty(allUnitCurves)
        ensZ = nan(1, nT);
        return;
    end

    ens  = mean(allUnitCurves, 1);
    ensZ = zscoreTrial_(ens);
end

function isCorr = getCorrectVectorFromTrialsCell_(trials)
    % User-stated structure:
    %   trials is a CELL (one cell per session), and each cell contains a TABLE
    %   with column "correct" (rows = trials)
    %
    % BUT in your prior scripts, session.trials was iterated as numel(trials) and
    % accessed as trials(k) (struct array). To support BOTH without changing upstream code,
    % we implement flexible parsing:
    %
    % If trials is a cell array AND trials{1} is a table with var 'correct', we:
    %   - concatenate correct across cells if each cell is a table, OR
    %   - if trials itself is a single cell holding a table, use it.
    isCorr = [];

    if ~iscell(trials), return; end
    if isempty(trials), return; end

    % Case A: trials is {table} (single cell)
    if numel(trials) == 1 && istable(trials{1}) && any(strcmp(trials{1}.Properties.VariableNames,'correct'))
        v = trials{1}.correct;
        isCorr = double(v(:));
        isCorr(~isfinite(isCorr)) = nan;
        return;
    end

    % Case B: trials is cell array of per-trial structs (not tables) -> not handled here
    if ~istable(trials{1})
        return;
    end

    % Case C: cell array of tables; stack correct
    cc = [];
    for i = 1:numel(trials)
        if ~istable(trials{i}), continue; end
        if ~any(strcmp(trials{i}.Properties.VariableNames,'correct')), continue; end
        v = double(trials{i}.correct(:));
        cc = [cc; v]; %#ok<AGROW>
    end
    if ~isempty(cc)
        isCorr = cc;
    end
end

function isCorr = getCorrectVectorFromTrialsStruct_(trials)
    % Fallback if trials is a struct array with field 'correct'
    isCorr = [];
    if isstruct(trials) && isfield(trials,'correct')
        v = [trials.correct];
        isCorr = double(v(:));
    end
end

function tr = getTrialStruct_(trials, k)
    % Flexible trial accessor:
    % - If trials is struct array: return trials(k)
    % - If trials is cell array of structs: return trials{k}
    % - If trials is a cell array containing a table (user-described): cannot return struct;
    %   then return [] (trial timing fields not available in that table-only representation).
    tr = [];
    if isstruct(trials)
        if k <= numel(trials), tr = trials(k); end
        return;
    end
    if iscell(trials)
        if k <= numel(trials) && isstruct(trials{k})
            tr = trials{k};
            return;
        end
        % if it's a cell-of-table, we cannot retrieve cue/press/lick per trial from here
        tr = [];
        return;
    end
end

function z = zscoreTrial_(x)
    x = x(:)';
    mu = mean(x);
    sd = std(x);
    if ~isfinite(sd) || sd <= 0, sd = 1; end
    z = (x - mu) / sd;
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

function tf = isGoodTrialsStructFlexible_(session)
    tf = false;
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    % allow struct OR cell; timing fields checked per-trial
    tf = isstruct(session.trials) || iscell(session.trials);
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