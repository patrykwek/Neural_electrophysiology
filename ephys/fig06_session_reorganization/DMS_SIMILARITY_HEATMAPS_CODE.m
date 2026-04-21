%% ===== DMS Similarity Heatmaps (Cue / Press / Lick) =====
% Session-wise population vectors (paper-style within-session population geometry),
% averaged across Early vs Late sessions, plotted as 6x6 block matrices.
%
% CHANGE (requested, only):
% - Equalize trial counts across L/C/R within each session BEFORE bootstrapping:
%   randomly subsample each cueType to nMin = min(nL,nC,nR), using correct-only trials.
%   Then bootstrap within these balanced trial sets.
%
% NEW CHANGE (requested, only):
% - Use ONLY CLASSIFIED units (cell_type in {1,2,3}) (empties and 0=Unclassified skipped)
%
% No other unnecessary changes.

clear; clc;

%% ---- SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/DMS_SIMILARITY_HEATMAPS';
if ~exist(outDir,'dir'), mkdir(outDir); end

% >>> ADDED: classification file (ONLY classified units) <<<
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skipped)

% Session groups (your convention)
earlySess = 1:8;
lateSess  = 29:36;

% Trial filters
MIN_RT_MS = 100;

% Baseline window (relative to cue) for z-scoring, ms
BASELINE_WIN_MS = [-100 -1];

% --- BOOTSTRAP SETTINGS ---
nBoot   = 10;   % number of bootstrap resamples per session
rngSeed = 0;     % reproducible

rng(rngSeed);    % seed ONCE

% -------------------------------------------------------------------------
% LOAD GLOBAL WINDOWS (ROBUST PATH SEARCH)
% -------------------------------------------------------------------------
baseOut_forWins = '/Volumes/WD_BLACK/A_THESIS_FINAL'; % keep consistent with your other scripts
candidates = { ...
    fullfile(baseOut_forWins,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut_forWins,'GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(outDir,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(fileparts(outDir),'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
    };

winFile = '';
for i = 1:numel(candidates)
    if exist(candidates{i},'file') == 2
        winFile = candidates{i};
        break;
    end
end

assert(~isempty(winFile), ...
    'Missing GLOBAL windows file. Looked in:\n  %s', strjoin(candidates, sprintf('\n  ')));

W = load(winFile);
assert(isfield(W,'GLOBAL_WINDOWS') && isstruct(W.GLOBAL_WINDOWS), ...
    'GLOBAL_WINDOWS missing in: %s', winFile);
GW = W.GLOBAL_WINDOWS;

fprintf('\nLoaded GLOBAL windows from:\n  %s\n', winFile);

% Windows are stored as [start end] in ms RELATIVE TO EACH EVENT (lag axis).
WIN_CUE_MS   = GW.cue_ms;    % relative to cue
WIN_PRESS_MS = GW.press_ms;  % relative to press
WIN_LICK_MS  = GW.lick_ms;   % relative to lick

fprintf('\n=== WINDOWS (RELATIVE TO EACH EVENT) ===\n');
fprintf('Cue window (rel cue):     [%d, %d] ms\n', round(WIN_CUE_MS(1)),   round(WIN_CUE_MS(2)));
fprintf('Press window (rel press): [%d, %d] ms\n', round(WIN_PRESS_MS(1)), round(WIN_PRESS_MS(2)));
fprintf('Lick window (rel lick):   [%d, %d] ms\n', round(WIN_LICK_MS(1)),  round(WIN_LICK_MS(2)));
% -------------------------------------------------------------------------

% Unit inclusion (per session)
minTrialsPerCond = 5;  % per cueType within session; unit must pass for ALL L/C/R to enter vectors

% Plot style (match your scripts)
figPos  = [120, 120, 1350, 520];
fsTitle = 22;
fsAx    = 16;
lwAx    = 2.5;
cmap    = parula(256);

cats = ["L","C","R"];
condLabels = ["Early-L","Early-C","Early-R","Late-L","Late-C","Late-R"];

%% ---- LOAD ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);
fprintf('Loaded beh with %d sessions.\n', nSessions);

% Clamp session lists
earlySess = earlySess(earlySess>=1 & earlySess<=nSessions);
lateSess  = lateSess(lateSess>=1  & lateSess<=nSessions);

%% ---- BUILD 6x6 MATRICES (SESSION-WISE POPULATION VECTORS; SPLIT-BOOTSTRAP; BALANCED TRIALS) ----
Rcue   = buildSessionWise6x6_(beh, S, cell_type, earlySess, lateSess, cats, 'cue',   WIN_CUE_MS,   BASELINE_WIN_MS, MIN_RT_MS, minTrialsPerCond, nBoot);
Rpress = buildSessionWise6x6_(beh, S, cell_type, earlySess, lateSess, cats, 'press', WIN_PRESS_MS, BASELINE_WIN_MS, MIN_RT_MS, minTrialsPerCond, nBoot);
Rlick  = buildSessionWise6x6_(beh, S, cell_type, earlySess, lateSess, cats, 'lick',  WIN_LICK_MS,  BASELINE_WIN_MS, MIN_RT_MS, minTrialsPerCond, nBoot);

%% ---- PRINT MATRICES ----
fprintf('\n============================================================\n');
fprintf('DMS similarity matrices (6x6): session-wise PV; zFR; correct-only; balanced L/C/R; split-bootstrap n=%d\n', nBoot);
fprintf('Labels: %s\n', strjoin(cellstr(condLabels), ', '));
fprintf('============================================================\n');

printMatrix_(Rcue,   condLabels, sprintf('CUE window   [%d %d] ms rel cue',   round(WIN_CUE_MS(1)),   round(WIN_CUE_MS(2))));
printMatrix_(Rpress, condLabels, sprintf('PRESS window [%d %d] ms rel press', round(WIN_PRESS_MS(1)), round(WIN_PRESS_MS(2))));
printMatrix_(Rlick,  condLabels, sprintf('LICK window  [%d %d] ms rel lick',  round(WIN_LICK_MS(1)),  round(WIN_LICK_MS(2))));

%% ---- PLOT: 3 HEATMAPS ----
fig = figure('Color','w','Position',figPos);

subplot(1,3,1);
plotHeat_(Rcue, condLabels, 'Cue window', fsTitle, fsAx, lwAx, cmap);

subplot(1,3,2);
plotHeat_(Rpress, condLabels, 'Press window', fsTitle, fsAx, lwAx, cmap);

subplot(1,3,3);
plotHeat_(Rlick, condLabels, 'Lick window', fsTitle, fsAx, lwAx, cmap);

sgtitle(sprintf('DMS similarity (session-wise PV; zFR; correct-only; balanced L/C/R; split-bootstrap n=%d) | Early %d-%d vs Late %d-%d', ...
    nBoot, earlySess(1), earlySess(end), lateSess(1), lateSess(end)), ...
    'FontSize', fsTitle, 'FontWeight','bold');

%% ---- SAVE ----
tag = sprintf('DMS_SimHeatmaps_SessionWiseZ_CORRECTONLY_BALANCED_SPLITBOOT%03d_Early%02d_%02d_Late%02d_%02d', ...
    nBoot, earlySess(1), earlySess(end), lateSess(1), lateSess(end));

outPng = fullfile(outDir, [tag '.png']);
outFig = fullfile(outDir, [tag '.fig']);
outSvg = fullfile(outDir, [tag '.svg']);

saveas(fig, outPng);
savefig(fig, outFig);
try
    print(fig, outSvg, '-dsvg');
catch ME
    warning('SVG save failed: %s', ME.message);
end

fprintf('\nSaved:\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg);

%% ========================= HELPERS =========================

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

function ct = safe_get_celltype_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
       ch   <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
       u    <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end

function R6 = buildSessionWise6x6_(beh, S, cell_type, earlySess, lateSess, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot)

    Re = [];
    Rl = [];

    ePool = struct('L',[],'C',[],'R',[]);
    lPool = struct('L',[],'C',[],'R',[]);

    % --- EARLY ---
    for sIdx = earlySess
        [Rboot, vL, vC, vR] = sessionPopulationVectors_splitBootCorr_balanced_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot);
        if isempty(Rboot), continue; end
        Re(:,:,end+1) = Rboot; %#ok<AGROW>
        ePool.L = [ePool.L; vL]; %#ok<AGROW>
        ePool.C = [ePool.C; vC]; %#ok<AGROW>
        ePool.R = [ePool.R; vR]; %#ok<AGROW>
    end

    % --- LATE ---
    for sIdx = lateSess
        [Rboot, vL, vC, vR] = sessionPopulationVectors_splitBootCorr_balanced_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot);
        if isempty(Rboot), continue; end
        Rl(:,:,end+1) = Rboot; %#ok<AGROW>
        lPool.L = [lPool.L; vL]; %#ok<AGROW>
        lPool.C = [lPool.C; vC]; %#ok<AGROW>
        lPool.R = [lPool.R; vR]; %#ok<AGROW>
    end

    assert(~isempty(Re), '%s: no usable EARLY sessions.', eventName);
    assert(~isempty(Rl), '%s: no usable LATE sessions.', eventName);

    R_early = mean(Re, 3, 'omitnan');
    R_late  = mean(Rl, 3, 'omitnan');

    % --- CROSS 3x3 (pooled + length-matched subsampling) ---
    Xearly = [ePool.L, ePool.C, ePool.R];
    Xlate  = [lPool.L, lPool.C, lPool.R];

    R_cross = nan(3,3);
    for i = 1:3
        for j = 1:3
            a = Xearly(:,i); a = a(isfinite(a));
            b = Xlate(:,j);  b = b(isfinite(b));
            n = min(numel(a), numel(b));
            assert(n > 10, '%s: too few pooled units for cross block.', eventName);
            ia = randperm(numel(a), n);
            ib = randperm(numel(b), n);
            R_cross(i,j) = corr(a(ia), b(ib), 'Rows','pairwise');
        end
    end

    R6 = nan(6,6);
    R6(1:3,1:3) = R_early;
    R6(4:6,4:6) = R_late;
    R6(1:3,4:6) = R_cross;
    R6(4:6,1:3) = R_cross';
end

function [Rboot, vL, vC, vR] = sessionPopulationVectors_splitBootCorr_balanced_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot)
% Requested change:
% - Balance L/C/R trial counts within this session:
%     nMin = min(nL,nC,nR)
%   Randomly subsample each cueType to nMin (WITHOUT replacement) to define
%   the working set of trials.
% - Then do split-bootstrap WITH replacement within each balanced set.
%
% NEW change:
% - ONLY CLASSIFIED units (cell_type in {1,2,3})

    Rboot = [];
    vL = []; vC = []; vR = [];

    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), return; end

    % ---- trial keep mask (correct-only) ----
    nT = numel(trials);
    keep = false(nT,1);
    cueType = strings(nT,1);
    cueAbs  = nan(nT,1);
    pressAbs= nan(nT,1);
    lickAbs = nan(nT,1);

    for k = 1:nT
        tr = trials(k);

        if ~isfield(tr,'valid') || ~tr.valid, continue; end
        if ~isfield(tr,'correct') || ~isfinite(double(tr.correct)) || double(tr.correct) ~= 1
            continue;
        end
        if ~all(isfield(tr, {'cue','press','lick','cueType'})), continue; end

        cue   = double(tr.cue);
        press = double(tr.press);
        lick  = double(tr.lick);
        if any(~isfinite([cue press lick])), continue; end

        rt = press - cue;
        if rt < MIN_RT_MS, continue; end

        ct = upper(strtrim(string(tr.cueType)));
        if ~ismember(ct, cats), continue; end

        keep(k) = true;
        cueType(k) = ct;
        cueAbs(k)  = cue;
        pressAbs(k)= press;
        lickAbs(k) = lick;
    end

    idxKeep = find(keep);
    if isempty(idxKeep), return; end

    switch lower(eventName)
        case 'cue'
            evtAbs = cueAbs;
        case 'press'
            evtAbs = pressAbs;
        case 'lick'
            evtAbs = lickAbs;
        otherwise
            error('Unknown eventName: %s', eventName);
    end

    ctKeep = cueType(idxKeep);

    idxL_all = find(ctKeep=="L");
    idxC_all = find(ctKeep=="C");
    idxR_all = find(ctKeep=="R");

    if numel(idxL_all) < minTrialsPerCond || numel(idxC_all) < minTrialsPerCond || numel(idxR_all) < minTrialsPerCond
        return;
    end

    % ---- BALANCE TRIAL COUNTS ACROSS L/C/R (requested) ----
    nMin = min([numel(idxL_all), numel(idxC_all), numel(idxR_all)]);
    if nMin < minTrialsPerCond
        return;
    end

    % Random subsample WITHOUT replacement to define balanced trial pools
    idxL = idxL_all(randperm(numel(idxL_all), nMin));
    idxC = idxC_all(randperm(numel(idxC_all), nMin));
    idxR = idxR_all(randperm(numel(idxR_all), nMin));

    % ---- collect per-unit zEvt across kept trials ----
    unitZ = {};
    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)

            % >>> ONLY CLASSIFIED UNITS (1/2/3) <<<
            ct = safe_get_celltype_(cell_type, sIdx, ch, u);
            if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                continue;
            end

            spk_abs = uc{u};
            if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
            spk_abs = double(spk_abs(:));

            baseFR = nan(numel(idxKeep),1);
            evtFR  = nan(numel(idxKeep),1);

            for ii = 1:numel(idxKeep)
                t = idxKeep(ii);
                cue0 = cueAbs(t);
                evt0 = evtAbs(t);

                b0 = cue0 + winBase(1);
                b1 = cue0 + winBase(2);
                if ~(isfinite(b0) && isfinite(b1) && b1 > b0), continue; end
                nB = nnz(spk_abs >= b0 & spk_abs <= b1);
                baseFR(ii) = nB / ((b1-b0)/1000);

                w0 = evt0 + winEvt(1);
                w1 = evt0 + winEvt(2);
                if ~(isfinite(w0) && isfinite(w1) && w1 > w0), continue; end
                nE = nnz(spk_abs >= w0 & spk_abs <= w1);
                evtFR(ii) = nE / ((w1-w0)/1000);
            end

            goodB = isfinite(baseFR);
            if nnz(goodB) < 5, continue; end

            muB = mean(baseFR(goodB), 'omitnan');
            sdB = std(baseFR(goodB), 0, 'omitnan');
            if ~isfinite(sdB) || sdB <= 0, continue; end

            zEvt = (evtFR - muB) ./ sdB;

            if nnz(isfinite(zEvt(idxL))) < minTrialsPerCond, continue; end
            if nnz(isfinite(zEvt(idxC))) < minTrialsPerCond, continue; end
            if nnz(isfinite(zEvt(idxR))) < minTrialsPerCond, continue; end

            unitZ{end+1,1} = zEvt; %#ok<AGROW>
        end
    end

    nUnits = numel(unitZ);
    if nUnits < 5, return; end

    % Non-boot pooled vectors (for cross-block pooling) using BALANCED trials
    vL = nan(nUnits,1); vC = nan(nUnits,1); vR = nan(nUnits,1);
    for k = 1:nUnits
        z = unitZ{k};
        vL(k) = mean(z(idxL), 'omitnan');
        vC(k) = mean(z(idxC), 'omitnan');
        vR(k) = mean(z(idxR), 'omitnan');
    end
    if any(~isfinite([vL;vC;vR])), return; end

    % ---- Split-bootstrap correlation matrices (WITH replacement within balanced sets) ----
    Rstack = nan(3,3,nBoot);

    for b = 1:nBoot
        % independent resamples for set-1
        rL1 = idxL(randi(numel(idxL), numel(idxL), 1));
        rC1 = idxC(randi(numel(idxC), numel(idxC), 1));
        rR1 = idxR(randi(numel(idxR), numel(idxR), 1));

        % independent resamples for set-2
        rL2 = idxL(randi(numel(idxL), numel(idxL), 1));
        rC2 = idxC(randi(numel(idxC), numel(idxC), 1));
        rR2 = idxR(randi(numel(idxR), numel(idxR), 1));

        A = nan(nUnits,3); % set-1 vectors for [L C R]
        B = nan(nUnits,3); % set-2 vectors for [L C R]

        for k = 1:nUnits
            z = unitZ{k};
            A(k,1) = mean(z(rL1), 'omitnan');
            A(k,2) = mean(z(rC1), 'omitnan');
            A(k,3) = mean(z(rR1), 'omitnan');

            B(k,1) = mean(z(rL2), 'omitnan');
            B(k,2) = mean(z(rC2), 'omitnan');
            B(k,3) = mean(z(rR2), 'omitnan');
        end

        if any(~isfinite(A(:))) || any(~isfinite(B(:)))
            continue;
        end

        % Corr between independent estimates: 3x3
        Rb = corr(A, B, 'Rows','pairwise'); % columns of A vs columns of B

        if all(isfinite(Rb(:)))
            % keep symmetric similarity-style matrix
            Rb = 0.5*(Rb + Rb');
            Rstack(:,:,b) = Rb;
        end
    end

    if all(all(all(~isfinite(Rstack))))
        return;
    end

    Rboot = mean(Rstack, 3, 'omitnan');
end

function printMatrix_(R, labels, titleStr)
    fprintf('\n--- %s ---\n', titleStr);

    fprintf('%12s', '');
    for j = 1:numel(labels)
        fprintf('%12s', labels(j));
    end
    fprintf('\n');

    for i = 1:numel(labels)
        fprintf('%12s', labels(i));
        for j = 1:numel(labels)
            fprintf('%12.3f', R(i,j));
        end
        fprintf('\n');
    end
end

function plotHeat_(R, labels, ttl, fsTitle, fsAx, lwAx, cmap)
    imagesc(R);
    axis square;
    colormap(cmap);
    caxis([-1 1]);
    cb = colorbar;
    cb.TickDirection = 'out';
    cb.LineWidth = lwAx;

    xticks(1:numel(labels));
    yticks(1:numel(labels));
    xticklabels(labels);
    yticklabels(labels);
    xtickangle(45);

    set(gca, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    box off;

    title(ttl, 'FontSize', fsTitle, 'FontWeight','bold');
end