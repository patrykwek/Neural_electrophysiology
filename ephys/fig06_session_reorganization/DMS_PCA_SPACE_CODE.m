%% ===== PCA SCATTER (OPTION A, CLOUD FIX): EARLY vs LATE; TRIAL×STAGE; SHARED SUBSPACE =====
% KEY FIX vs your line-like plot:
%   - NO pseudo-pop union with NaN imputation.
%   - Instead: per-session population vectors -> per-session PCA -> concatenate in common K0 space -> global PCA.
%
% WHAT THIS DOES:
%   1) For each session, build trial×stage population vectors using ONLY units recorded that session:
%        x(trial,stage) = [dFR_unit1, dFR_unit2, ...] (no NaNs)
%        dFR = FR(stage window) - FR(baseline window)
%   2) Z-score features within-session; compress to K0 dims via PCA (session-local latent)
%   3) Stack all sessions' K0-latents and fit ONE global PCA (shared subspace)
%   4) Plot Early vs Late panels, colored by stage, with stage centroids
%
% RELAXED trial inclusion (keeps 30 & 34):
%   requireValid=false, requireCorrect=false, MIN_RT_MS=0
% Windows clipped to [t_start,t_end] when present.

clear; clc;

%% ---- SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);

outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/DECODING_STAGE_SCALARONLY_RELAXED_WITHIN_GROUP';
if ~exist(outDir,'dir'), mkdir(outDir); end

% classification file
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;

% groups
earlySess = 1:8;
lateSess  = 29:36;

% RELAXED inclusion
requireValid   = false;
requireCorrect = false;
MIN_RT_MS      = 0;

% baseline relative to cue (ms)
BASELINE_WIN_MS = [-500 0];

% GLOBAL event windows file (robust search)
baseOut_forWins = '/Volumes/WD_BLACK/A_THESIS_FINAL';
candidates = { ...
    fullfile(baseOut_forWins,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'POPULATION_CHANGE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'DMS_SIMILARITY_HEATMAPS','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
    };

winFile = '';
for i = 1:numel(candidates)
    if exist(candidates{i},'file') == 2
        winFile = candidates{i};
        break;
    end
end

% PCA settings
N_PCS_GLOBAL = 3;

% K0 = per-session latent dimensionality BEFORE global PCA (cloud fix)
% We'll automatically set K0 = min(15, minUnitsAcrossSess-1), after we compute session unit counts.
K0_cap = 15;

% plotting / style
figPos = [220 140 1450 650];

lwAx    = 4.0;
fsAx    = 24;
fsTitle = 26;

% stage colors
colCue   = [0.85 0.10 0.10];   % red
colPress = [0.10 0.25 0.95];   % blue
colLick  = [0.10 0.75 0.20];   % green

msPt    = 8;
alphaPt = 0.20;
msCent  = 18;
lwCent  = 3.0;

% downsample points per panel per stage
maxPtsPerStagePerPanel = 6000;

rng(0);

% outputs
tagBase = sprintf('PCA_SCATTER_OPTIONA_CLOUDFIX_TRIALxSTAGE_SHARED_GLOBALSPACE_RELAXED_CLASSIFIEDONLY_Early%02d_%02d_Late%02d_%02d', ...
    earlySess(1), earlySess(end), lateSess(1), lateSess(end));

outPng = fullfile(outDir, [tagBase '.png']);
outFig = fullfile(outDir, [tagBase '.fig']);
outSvg = fullfile(outDir, [tagBase '.svg']);

%% ---- LOAD DATA ----
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);
fprintf('Loaded beh with %d sessions.\n', nSessions);

earlySess = earlySess(earlySess>=1 & earlySess<=nSessions);
lateSess  = lateSess(lateSess>=1  & lateSess<=nSessions);

sessUse = unique([earlySess(:); lateSess(:)]);
assert(~isempty(sessUse), 'No sessions selected.');

%% ---- LOAD WINDOWS ----
if ~isempty(winFile)
    W = load(winFile);
    assert(isfield(W,'GLOBAL_WINDOWS') && isstruct(W.GLOBAL_WINDOWS), 'GLOBAL_WINDOWS missing in: %s', winFile);
    GW = W.GLOBAL_WINDOWS;
    WIN_CUE_MS   = GW.cue_ms;
    WIN_PRESS_MS = GW.press_ms;
    WIN_LICK_MS  = GW.lick_ms;
    fprintf('\nLoaded GLOBAL windows from:\n  %s\n', winFile);
else
    WIN_CUE_MS   = [0 300];
    WIN_PRESS_MS = [-100 200];
    WIN_LICK_MS  = [-100 200];
    fprintf('\nGLOBAL windows file not found. Using defaults.\n');
end

fprintf('\n=== WINDOWS (RELATIVE TO EVENT) ===\n');
fprintf('Baseline (rel cue): [%d, %d] ms\n', round(BASELINE_WIN_MS(1)), round(BASELINE_WIN_MS(2)));
fprintf('Cue   (rel cue):    [%d, %d] ms\n', round(WIN_CUE_MS(1)),   round(WIN_CUE_MS(2)));
fprintf('Press (rel press):  [%d, %d] ms\n', round(WIN_PRESS_MS(1)), round(WIN_PRESS_MS(2)));
fprintf('Lick  (rel lick):   [%d, %d] ms\n', round(WIN_LICK_MS(1)),  round(WIN_LICK_MS(2)));

%% ============================================================
%  PASS 1: build per-session trial×stage matrices (no NaNs), track unit counts
%% ============================================================
Sess = struct();
unitCounts = nan(numel(sessUse),1);

for ii = 1:numel(sessUse)
    sIdx = sessUse(ii);

    [Xraw, yStage] = buildTrialStageVectors_oneSession_LOCAL_( ...
        beh, S, cell_type, sIdx, ...
        BASELINE_WIN_MS, WIN_CUE_MS, WIN_PRESS_MS, WIN_LICK_MS, ...
        requireValid, requireCorrect, MIN_RT_MS);

    if isempty(Xraw)
        fprintf('Session %02d: no usable samples.\n', sIdx);
        continue;
    end

    Sess(ii).sIdx   = sIdx;
    Sess(ii).Xraw   = Xraw;    % [nSamples x nUnitsThisSession]
    Sess(ii).yStage = yStage;  % [nSamples x 1], 1/2/3
    unitCounts(ii)  = size(Xraw,2);

    fprintf('Session %02d: samples=%d, units=%d\n', sIdx, size(Xraw,1), size(Xraw,2));
end

% keep only sessions that produced data
okSess = arrayfun(@(q) isfield(q,'Xraw') && ~isempty(q.Xraw), Sess);
Sess = Sess(okSess);
sessUse2 = arrayfun(@(q) q.sIdx, Sess);

assert(~isempty(Sess), 'No sessions produced any samples.');

unitCounts2 = arrayfun(@(q) size(q.Xraw,2), Sess);
minUnits = min(unitCounts2);
fprintf('\nMin units across usable sessions: %d\n', minUnits);

% choose K0 safely
K0 = min(K0_cap, max(2, minUnits-1));
fprintf('Using K0=%d per-session latent dims (then global PCA to %d PCs)\n', K0, N_PCS_GLOBAL);

%% ============================================================
%  PASS 2: within-session PCA to K0, then concatenate into common K0 space
%% ============================================================
Zall = [];      % [Ntotal x K0] concatenated session-latents
yStageAll = []; % [Ntotal x 1]
yGroupAll = []; % [Ntotal x 1] 1=early,2=late
sIdAll    = []; % [Ntotal x 1] session index (for debugging)

for ii = 1:numel(Sess)
    sIdx = Sess(ii).sIdx;
    X = Sess(ii).Xraw;
    ySt = Sess(ii).yStage;

    % z-score features within session (unit-wise)
    Xz = zscore(X, 0, 1);
    Xz(~isfinite(Xz)) = 0;

    % session PCA -> scores (latent)
    % note: we do not care about coeff identity across sessions; we only need a common K0 coordinate count.
    [~, score, ~] = pca(Xz, 'Algorithm','svd', 'NumComponents', K0);

    Zall = [Zall; score(:,1:K0)]; %#ok<AGROW>
    yStageAll = [yStageAll; ySt(:)]; %#ok<AGROW>

    if ismember(sIdx, earlySess)
        g = 1;
    elseif ismember(sIdx, lateSess)
        g = 2;
    else
        continue;
    end
    yGroupAll = [yGroupAll; g*ones(size(ySt(:)))]; %#ok<AGROW>
    sIdAll    = [sIdAll; sIdx*ones(size(ySt(:)))]; %#ok<AGROW>
end

assert(~isempty(Zall), 'No latent samples collected.');

fprintf('\nTOTAL latent samples: %d (each = one trial×stage)\n', size(Zall,1));
fprintf('Stage counts: cue=%d press=%d lick=%d\n', nnz(yStageAll==1), nnz(yStageAll==2), nnz(yStageAll==3));
fprintf('Group counts: early=%d late=%d\n', nnz(yGroupAll==1), nnz(yGroupAll==2));

%% ============================================================
%  GLOBAL PCA (shared subspace) on common K0-latent space
%% ============================================================
[coeffG, scoreG, ~, ~, explainedG, muG] = pca(Zall, 'Algorithm','svd', 'NumComponents', max(3,N_PCS_GLOBAL));
PC = scoreG(:,1:N_PCS_GLOBAL);

fprintf('\nGLOBAL PCA explained: PC1=%.2f%% PC2=%.2f%% PC3=%.2f%%\n', ...
    explainedG(1), explainedG(2), explainedG(3));

%% ---- SUBSAMPLE FOR PLOTTING (per panel per stage) ----
keepPlot = false(size(PC,1),1);
for g = 1:2
    for st = 1:3
        idx = find(yGroupAll==g & yStageAll==st);
        if isempty(idx), continue; end
        if numel(idx) > maxPtsPerStagePerPanel
            idx = idx(randperm(numel(idx), maxPtsPerStagePerPanel));
        end
        keepPlot(idx) = true;
    end
end

PCp = PC(keepPlot,:);
ysP = yStageAll(keepPlot);
ygP = yGroupAll(keepPlot);

%% ---- PLOT: TWO PANELS ----
fig = figure('Color','w','Position',figPos);

% shared limits
xl = [min(PCp(:,1)) max(PCp(:,1))];
yl = [min(PCp(:,2)) max(PCp(:,2))];
zl = [min(PCp(:,3)) max(PCp(:,3))];

padFrac = 0.06;
xl = xl + [-1 1]*padFrac*rangeSafe_(xl);
yl = yl + [-1 1]*padFrac*rangeSafe_(yl);
zl = zl + [-1 1]*padFrac*rangeSafe_(zl);

for panel = 1:2
    if panel==1
        panelTitle = sprintf('Early sessions (%d-%d)', earlySess(1), earlySess(end));
        grpVal = 1;
    else
        panelTitle = sprintf('Late sessions (%d-%d)', lateSess(1), lateSess(end));
        grpVal = 2;
    end

    ax = subplot(1,2,panel);
    hold(ax,'on');

    plotStageScatter_(ax, PCp, ysP, ygP, grpVal, 1, colCue,   msPt, alphaPt);
    plotStageScatter_(ax, PCp, ysP, ygP, grpVal, 2, colPress, msPt, alphaPt);
    plotStageScatter_(ax, PCp, ysP, ygP, grpVal, 3, colLick,  msPt, alphaPt);

    % centroids (computed on ALL points, not subsampled)
    plotStageCentroids_(ax, PC, yStageAll, yGroupAll, grpVal, colCue, colPress, colLick, msCent, lwCent);

    xlabel(ax, sprintf('PC1 (%.1f%%)', explainedG(1)));
    ylabel(ax, sprintf('PC2 (%.1f%%)', explainedG(2)));
    zlabel(ax, sprintf('PC3 (%.1f%%)', explainedG(3)));

    title(ax, panelTitle, 'FontWeight','bold', 'FontSize', fsTitle);

    set(ax, 'LineWidth', lwAx, 'TickDir','out', 'FontSize', fsAx);
    grid(ax,'on'); box(ax,'off');
    view(ax, 35, 20);

    xlim(ax, xl); ylim(ax, yl); zlim(ax, zl);
end

% legend proxies
axL = subplot(1,2,2);
h = gobjects(0,1); lab = strings(0,1);
h(end+1) = plot3(axL, nan,nan,nan, 'o', 'MarkerFaceColor', colCue,   'MarkerEdgeColor', colCue);   lab(end+1) = "Cue";
h(end+1) = plot3(axL, nan,nan,nan, 'o', 'MarkerFaceColor', colPress, 'MarkerEdgeColor', colPress); lab(end+1) = "Press";
h(end+1) = plot3(axL, nan,nan,nan, 'o', 'MarkerFaceColor', colLick,  'MarkerEdgeColor', colLick);  lab(end+1) = "Lick";
h(end+1) = plot3(axL, nan,nan,nan, 'p', 'MarkerFaceColor', [0 0 0],  'MarkerEdgeColor', [0 0 0]);  lab(end+1) = "Stage centroid";
legend(axL, h, cellstr(lab), 'Location','bestoutside');

sgtitle('PCA of trial×stage population vectors (cloud-fix: session PCA -> shared global PCA)', ...
    'FontSize', fsTitle, 'FontWeight','bold');

%% ---- SAVE ----
saveas(fig, outPng);
savefig(fig, outFig);
try, print(fig, outSvg, '-dsvg'); catch, end
fprintf('\nSaved:\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg);

fprintf('\nDONE.\n');

%% ========================= HELPERS =========================

function r = rangeSafe_(ab)
    r = ab(2)-ab(1);
    if ~isfinite(r) || r<=0, r = 1; end
end

function plotStageScatter_(ax, PC, yStage, yGroup, grpVal, stageVal, col, msPt, alphaPt)
    idx = (yGroup==grpVal) & (yStage==stageVal);
    if ~any(idx), return; end
    try
        s = scatter3(ax, PC(idx,1), PC(idx,2), PC(idx,3), msPt, 'filled', ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
        s.MarkerFaceAlpha = alphaPt;
        s.MarkerEdgeAlpha = alphaPt;
    catch
        plot3(ax, PC(idx,1), PC(idx,2), PC(idx,3), 'o', ...
            'MarkerSize', max(2, round(msPt/3)), ...
            'MarkerFaceColor', col, 'MarkerEdgeColor', col);
    end
end

function plotStageCentroids_(ax, PC, yStage, yGroup, grpVal, colCue, colPress, colLick, msCent, lwCent)
    cols = {colCue, colPress, colLick};
    for st = 1:3
        idx = (yGroup==grpVal) & (yStage==st);
        if ~any(idx), continue; end
        c = mean(PC(idx,1:3), 1, 'omitnan');
        plot3(ax, c(1), c(2), c(3), 'p', ...
            'MarkerSize', msCent, ...
            'MarkerFaceColor', [0 0 0], ...
            'MarkerEdgeColor', cols{st}, ...
            'LineWidth', lwCent);
    end
end

function [Xraw, yStage] = buildTrialStageVectors_oneSession_LOCAL_( ...
    beh, S, cell_type, sIdx, ...
    BASELINE_WIN_MS, WIN_CUE_MS, WIN_PRESS_MS, WIN_LICK_MS, ...
    requireValid, requireCorrect, MIN_RT_MS)

    Xraw = [];
    yStage = [];

    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials) || ~isstruct(session.trials)
        return;
    end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
        return;
    end

    % --- build local unit list for this session (classified-only) ---
    % store spikes in a flat list for speed
    spkList = cell(0,1);
    chList  = [];
    uList   = [];

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end
        for u = 1:numel(uc)
            ct = safe_get_celltype_(cell_type, sIdx, ch, u);
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                continue;
            end
            spk = uc{u};
            if isempty(spk) || ~isnumeric(spk), continue; end
            spk = double(spk(:));
            if isempty(spk) || ~any(isfinite(spk)), continue; end
            spkList{end+1,1} = spk; %#ok<AGROW>
            chList(end+1,1) = ch; %#ok<AGROW>
            uList(end+1,1)  = u; %#ok<AGROW>
        end
    end

    nU = numel(spkList);
    if nU < 2
        return;
    end

    % --- build 3 samples per trial: cue/press/lick ---
    for k = 1:numel(trials)
        tr = trials(k);

        if requireValid
            if ~isfield(tr,'valid') || ~tr.valid, continue; end
        end
        if requireCorrect
            if ~isfield(tr,'correct') || ~isfinite(double(tr.correct)) || double(tr.correct) ~= 1
                continue;
            end
        end

        if ~all(isfield(tr, {'cue','press','lick'})), continue; end
        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end

        % optional bounds for clipping
        tStart = -inf; tEnd = inf;
        if isfield(tr,'t_start') && isfinite(double(tr.t_start)), tStart = double(tr.t_start); end
        if isfield(tr,'t_end')   && isfinite(double(tr.t_end)),   tEnd   = double(tr.t_end);   end

        % baseline window relative cue (clipped)
        b0 = cue + BASELINE_WIN_MS(1);
        b1 = cue + BASELINE_WIN_MS(2);
        [b0c, b1c] = clipWin_(b0, b1, tStart, tEnd);
        if ~isfinite(b0c) || ~isfinite(b1c) || b1c <= b0c
            continue;
        end

        baseHz = nan(nU,1);
        for iU = 1:nU
            baseHz(iU) = rateHz_(spkList{iU}, b0c, b1c);
        end

        % CUE sample
        [w0, w1] = clipWin_(cue + WIN_CUE_MS(1), cue + WIN_CUE_MS(2), tStart, tEnd);
        if isfinite(w0) && isfinite(w1) && w1 > w0
            x = nan(1,nU);
            for iU = 1:nU
                stHz = rateHz_(spkList{iU}, w0, w1);
                x(iU) = stHz - baseHz(iU);
            end
            if nnz(isfinite(x)) >= 2
                Xraw(end+1,:) = x; %#ok<AGROW>
                yStage(end+1,1) = 1; %#ok<AGROW>
            end
        end

        % PRESS sample
        [w0, w1] = clipWin_(press + WIN_PRESS_MS(1), press + WIN_PRESS_MS(2), tStart, tEnd);
        if isfinite(w0) && isfinite(w1) && w1 > w0
            x = nan(1,nU);
            for iU = 1:nU
                stHz = rateHz_(spkList{iU}, w0, w1);
                x(iU) = stHz - baseHz(iU);
            end
            if nnz(isfinite(x)) >= 2
                Xraw(end+1,:) = x; %#ok<AGROW>
                yStage(end+1,1) = 2; %#ok<AGROW>
            end
        end

        % LICK sample
        [w0, w1] = clipWin_(lick + WIN_LICK_MS(1), lick + WIN_LICK_MS(2), tStart, tEnd);
        if isfinite(w0) && isfinite(w1) && w1 > w0
            x = nan(1,nU);
            for iU = 1:nU
                stHz = rateHz_(spkList{iU}, w0, w1);
                x(iU) = stHz - baseHz(iU);
            end
            if nnz(isfinite(x)) >= 2
                Xraw(end+1,:) = x; %#ok<AGROW>
                yStage(end+1,1) = 3; %#ok<AGROW>
            end
        end
    end

    % drop rows with too many NaNs (just in case)
    if ~isempty(Xraw)
        ok = sum(isfinite(Xraw),2) >= 2;
        Xraw = Xraw(ok,:);
        yStage = yStage(ok,:);
    end
end

function hz = rateHz_(spk_abs, t0, t1)
    if ~(isfinite(t0) && isfinite(t1)) || t1 <= t0
        hz = nan; return;
    end
    n = nnz(spk_abs >= t0 & spk_abs <= t1);
    hz = n / ((t1 - t0)/1000);
end

function [t0c, t1c] = clipWin_(t0, t1, tStart, tEnd)
    t0c = t0; t1c = t1;
    if isfinite(tStart), t0c = max(t0c, tStart); end
    if isfinite(tEnd),   t1c = min(t1c, tEnd);   end
end

function ct = safe_get_celltype_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
       ch   <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
       u    <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k = 1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k})), beh = S.(cands{k}); return; end
    end
    f = fieldnames(S);
    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v) && numel(v)>1, beh = v; return; end
    end
    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v), beh = v; return; end
    end
end

function spikes_by_ch = getSpikesForSession_(session, S, sessionIdx)
    spikes_by_ch = [];
    if isfield(session,'spikes') && ~isempty(session.spikes)
        spikes_by_ch = session.spikes; return
    end
    if isfield(S,'spikes_session') && ~isempty(S.spikes_session) && ...
            numel(S.spikes_session) >= sessionIdx && ~isempty(S.spikes_session{sessionIdx})
        spikes_by_ch = S.spikes_session{sessionIdx}; return
    end
    if isfield(S,'spikes_persession') && ~isempty(S.spikes_persession) && ...
            numel(S.spikes_persession) >= sessionIdx && ~isempty(S.spikes_persession{sessionIdx})
        spikes_by_ch = S.spikes_persession{sessionIdx}; return
    end
    if isfield(S,'spikes') && iscell(S.spikes) && ~isempty(S.spikes)
        spikes_by_ch = S.spikes;
    end
end