%% ===== DMS Similarity: SAME-vs-DIFF barplots (Cue / Press / Lick) =====
% Produces 3 separate figures (cue, press, lick).
%
% For each event:
%   - Re-compute session-wise 3x3 similarity matrices (Early block and Late block logic),
%     using:
%       * correct-only trials
%       * balanced L/C/R trial counts within session (subsample to nMin)
%       * split-bootstrap (nBoot) to estimate the 3x3 within-session similarity
%       * ONLY CLASSIFIED units (cell_type in {1,2,3})
%   - For each SESSION, compute:
%       SAME  = mean(diag(R3))
%       DIFF  = mean(upper off-diagonals: L-C, L-R, C-R)
%   - Plot 4 bars:
%       Early SAME, Early DIFF, Late SAME, Late DIFF
%     Y-axis: correlation coefficient r
%   - Student t-tests:
%       (1) Early SAME vs Early DIFF  (paired t-test across early sessions)
%       (2) Late  SAME vs Late  DIFF  (paired t-test across late sessions)
%       (3) Early SAME vs Late SAME   (two-sample t-test)
%       (4) Early DIFF vs Late DIFF   (two-sample t-test)
%
% No unnecessary changes beyond adding these barplots + needed recomputation.

clear; clc;

%% ---- SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
outDir  = '/Volumes/WD_BLACK/A_THESIS_FINAL/DMS_SIMILARITY_HEATMAPS';
if ~exist(outDir,'dir'), mkdir(outDir); end

% >>> REQUIRED: classification output (ONLY classified units) <<<
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
nBoot   = 10;   % bootstrap resamples per session
rngSeed = 0;
rng(rngSeed);

% -------------------------------------------------------------------------
% LOAD GLOBAL WINDOWS (ROBUST PATH SEARCH)
% -------------------------------------------------------------------------
baseOut_forWins = '/Volumes/WD_BLACK/A_THESIS_FINAL';
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

WIN_CUE_MS   = GW.cue_ms;
WIN_PRESS_MS = GW.press_ms;
WIN_LICK_MS  = GW.lick_ms;

fprintf('\n=== WINDOWS (RELATIVE TO EACH EVENT) ===\n');
fprintf('Cue window (rel cue):     [%d, %d] ms\n', round(WIN_CUE_MS(1)),   round(WIN_CUE_MS(2)));
fprintf('Press window (rel press): [%d, %d] ms\n', round(WIN_PRESS_MS(1)), round(WIN_PRESS_MS(2)));
fprintf('Lick window (rel lick):   [%d, %d] ms\n', round(WIN_LICK_MS(1)),  round(WIN_LICK_MS(2)));
% -------------------------------------------------------------------------

% Unit inclusion (per session)
minTrialsPerCond = 5;  % per cueType within session; unit must pass for ALL L/C/R to enter vectors

% Plot style (kept simple; tweak if you want)
figPos  = [120, 120, 900, 520];
fsTitle = 20;
fsAx    = 16;
lwAx    = 2.5;

cats = ["L","C","R"];

% Bar colors (4 bars): EarlySame, EarlyDiff, LateSame, LateDiff
colEarlySame = [0.10 0.55 0.95];
colEarlyDiff = [0.10 0.75 0.55];
colLateSame  = [0.90 0.20 0.10];
colLateDiff  = [0.95 0.70 0.10];

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

assert(~isempty(earlySess) && ~isempty(lateSess), 'Early/Late lists empty after clamping.');

%% ---- COMPUTE + PLOT FOR EACH EVENT ----
eventList = { ...
    struct('name','cue',   'win',WIN_CUE_MS,   'ttl','Cue window'), ...
    struct('name','press', 'win',WIN_PRESS_MS, 'ttl','Press window'), ...
    struct('name','lick',  'win',WIN_LICK_MS,  'ttl','Lick window') ...
    };

for ev = 1:numel(eventList)
    evName = eventList{ev}.name;
    evWin  = eventList{ev}.win;
    evTtl  = eventList{ev}.ttl;

    % Per-session SAME/DIFF metrics
    [E_same, E_diff, L_same, L_diff, usedEarly, usedLate] = ...
        computeSameDiff_bySession_(beh, S, cell_type, earlySess, lateSess, cats, ...
            evName, evWin, BASELINE_WIN_MS, MIN_RT_MS, minTrialsPerCond, nBoot);

    fprintf('\n=== %s: sessions used (early=%d, late=%d) ===\n', upper(evName), numel(usedEarly), numel(usedLate));

    % ---------- STATS (Student t-tests) ----------
    % 1) Early: SAME vs DIFF (paired)
    p_E_same_diff = nan; t_E = nan; df_E = nan;
    if numel(E_same) >= 2 && numel(E_diff) == numel(E_same)
        [~, p_E_same_diff, ~, st] = ttest(E_same, E_diff);
        t_E = st.tstat; df_E = st.df;
    end

    % 2) Late: SAME vs DIFF (paired)
    p_L_same_diff = nan; t_L = nan; df_L = nan;
    if numel(L_same) >= 2 && numel(L_diff) == numel(L_same)
        [~, p_L_same_diff, ~, st] = ttest(L_same, L_diff);
        t_L = st.tstat; df_L = st.df;
    end

    % 3) Early SAME vs Late SAME (two-sample)
    p_same_E_vs_L = nan; t_s = nan; df_s = nan;
    if numel(E_same) >= 2 && numel(L_same) >= 2
        [~, p_same_E_vs_L, ~, st] = ttest2(E_same, L_same);
        t_s = st.tstat; df_s = st.df;
    end

    % 4) Early DIFF vs Late DIFF (two-sample)
    p_diff_E_vs_L = nan; t_d = nan; df_d = nan;
    if numel(E_diff) >= 2 && numel(L_diff) >= 2
        [~, p_diff_E_vs_L, ~, st] = ttest2(E_diff, L_diff);
        t_d = st.tstat; df_d = st.df;
    end

    fprintf('\n--- %s t-tests ---\n', upper(evName));
    fprintf('Early SAME vs DIFF (paired): t(%s)=%.4f, p=%.6g\n', num2str(df_E), t_E, p_E_same_diff);
    fprintf('Late  SAME vs DIFF (paired): t(%s)=%.4f, p=%.6g\n', num2str(df_L), t_L, p_L_same_diff);
    fprintf('SAME: Early vs Late (2-sample): t(%s)=%.4f, p=%.6g\n', num2str(df_s), t_s, p_same_E_vs_L);
    fprintf('DIFF: Early vs Late (2-sample): t(%s)=%.4f, p=%.6g\n', num2str(df_d), t_d, p_diff_E_vs_L);

    % ---------- BARPLOT (4 bars) ----------
    fig = figure('Color','w','Position',figPos);
    ax = axes(fig); hold(ax,'on');

    vals = [mean(E_same,'omitnan'), mean(E_diff,'omitnan'), mean(L_same,'omitnan'), mean(L_diff,'omitnan')];
    ns   = [numel(E_same),         numel(E_diff),         numel(L_same),         numel(L_diff)];
    sems = [std(E_same,0,'omitnan')/sqrt(max(1,ns(1))), ...
            std(E_diff,0,'omitnan')/sqrt(max(1,ns(2))), ...
            std(L_same,0,'omitnan')/sqrt(max(1,ns(3))), ...
            std(L_diff,0,'omitnan')/sqrt(max(1,ns(4)))];

    xlab = categorical({'Early SAME','Early DIFF','Late SAME','Late DIFF'}, ...
                       {'Early SAME','Early DIFF','Late SAME','Late DIFF'}, ...
                       'Ordinal', true);

    b = bar(ax, xlab, vals);
    b.FaceColor = 'flat';
    b.CData(1,:) = colEarlySame;
    b.CData(2,:) = colEarlyDiff;
    b.CData(3,:) = colLateSame;
    b.CData(4,:) = colLateDiff;

    errorbar(ax, b.XEndPoints, vals, sems, 'k.', 'LineWidth', 2);

    % Jitter points (session-level)
    rng(0);
    jitter = 0.10;
    ptSize = 30;

    % plot session points under each bar (gray)
    plotJitter_(ax, b.XEndPoints(1), E_same, jitter, ptSize);
    plotJitter_(ax, b.XEndPoints(2), E_diff, jitter, ptSize);
    plotJitter_(ax, b.XEndPoints(3), L_same, jitter, ptSize);
    plotJitter_(ax, b.XEndPoints(4), L_diff, jitter, ptSize);

    ylabel(ax, 'Correlation r', 'FontSize', fsAx);
    title(ax, sprintf('%s: SAME vs DIFF cue geometry (session as unit)', evTtl), ...
        'FontSize', fsTitle, 'FontWeight','bold');

    set(ax, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
    box(ax,'off');

    % Reasonable y-lims with headroom for stars
    yAll = [E_same(:); E_diff(:); L_same(:); L_diff(:)];
    yAll = yAll(isfinite(yAll));
    if isempty(yAll), yAll = 0; end
    yMin = min(yAll); yMax = max(yAll);
    yRange = yMax - yMin;
    if ~isfinite(yRange) || yRange==0, yRange = 1; end
    ylim(ax, [max(-1, yMin - 0.10*yRange), min(1, yMax + 0.25*yRange)]);

    % ---------- Significance bars (stars) ----------
    % We draw up to 4 comparisons:
    %   A: Early SAME vs Early DIFF  (1 vs 2)
    %   B: Late  SAME vs Late  DIFF  (3 vs 4)
    %   C: Early SAME vs Late SAME   (1 vs 3)
    %   D: Early DIFF vs Late DIFF   (2 vs 4)
    comps = [1 2; 3 4; 1 3; 2 4];
    pvals = [p_E_same_diff; p_L_same_diff; p_same_E_vs_L; p_diff_E_vs_L];
    compNames = {'E SAME vs E DIFF','L SAME vs L DIFF','E SAME vs L SAME','E DIFF vs L DIFF'};

    baseY = max(vals + sems) + 0.08*yRange;
    hBar  = 0.03*yRange;
    stepY = 0.07*yRange;

    drawCount = 0;
    for k = 1:4
        p = pvals(k);
        if ~(isfinite(p) && p < 0.05), continue; end
        drawCount = drawCount + 1;

        iA = comps(k,1); iB = comps(k,2);
        xA = b.XEndPoints(iA);
        xB = b.XEndPoints(iB);
        yLevel = baseY + (drawCount-1)*stepY;

        plot(ax, [xA xA xB xB], [yLevel yLevel+hBar yLevel+hBar yLevel], 'k-', 'LineWidth', 2);

        if p < 0.001
            s = '***';
        elseif p < 0.01
            s = '**';
        else
            s = '*';
        end

        text(ax, mean([xA xB]), yLevel+hBar+0.01*yRange, s, ...
            'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
            'FontSize', fsTitle, 'FontWeight','bold');

        fprintf('%s: p=%.6g\n', compNames{k}, p);
    end

    % ---------- SAVE ----------
    tagBase = sprintf('DMS_SameDiff_Bars_%s_CORRECTONLY_BALANCED_CLASSIFIEDONLY_BOOT%03d_Early%02d_%02d_Late%02d_%02d', ...
        upper(evName), nBoot, earlySess(1), earlySess(end), lateSess(1), lateSess(end));

    outPng = fullfile(outDir, [tagBase '.png']);
    outFig = fullfile(outDir, [tagBase '.fig']);
    outSvg = fullfile(outDir, [tagBase '.svg']);

    saveas(fig, outPng);
    savefig(fig, outFig);
    try
        print(fig, outSvg, '-dsvg');
    catch
        saveas(fig, outSvg);
    end

    % Save underlying numbers too
    outMat = fullfile(outDir, [tagBase '.mat']);
    save(outMat, ...
        'E_same','E_diff','L_same','L_diff','usedEarly','usedLate', ...
        'evName','evWin','BASELINE_WIN_MS','MIN_RT_MS','minTrialsPerCond','nBoot', ...
        'p_E_same_diff','p_L_same_diff','p_same_E_vs_L','p_diff_E_vs_L');

    fprintf('Saved:\n  %s\n  %s\n  %s\n  %s\n', outPng, outFig, outSvg, outMat);
end

%% ========================= HELPERS =========================

function plotJitter_(ax, x0, y, jitter, ptSize)
    y = y(:);
    y = y(isfinite(y));
    if isempty(y), return; end
    xj = x0 + (rand(numel(y),1)-0.5)*2*jitter;
    scatter(ax, xj, y, ptSize, 'filled', ...
        'MarkerFaceColor', [0.35 0.35 0.35], 'MarkerFaceAlpha', 0.35, ...
        'MarkerEdgeAlpha', 0);
end

function [E_same, E_diff, L_same, L_diff, usedEarly, usedLate] = computeSameDiff_bySession_( ...
    beh, S, cell_type, earlySess, lateSess, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot)

    E_same = []; E_diff = [];
    L_same = []; L_diff = [];
    usedEarly = [];
    usedLate  = [];

    % Early sessions
    for sIdx = earlySess(:)'
        [R3] = sessionR3_splitBoot_balanced_classified_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot);
        if isempty(R3), continue; end
        [sameVal, diffVal] = R3_to_same_diff_(R3);
        if ~isfinite(sameVal) || ~isfinite(diffVal), continue; end
        E_same(end+1,1) = sameVal; %#ok<AGROW>
        E_diff(end+1,1) = diffVal; %#ok<AGROW>
        usedEarly(end+1,1) = sIdx; %#ok<AGROW>
    end

    % Late sessions
    for sIdx = lateSess(:)'
        [R3] = sessionR3_splitBoot_balanced_classified_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot);
        if isempty(R3), continue; end
        [sameVal, diffVal] = R3_to_same_diff_(R3);
        if ~isfinite(sameVal) || ~isfinite(diffVal), continue; end
        L_same(end+1,1) = sameVal; %#ok<AGROW>
        L_diff(end+1,1) = diffVal; %#ok<AGROW>
        usedLate(end+1,1) = sIdx; %#ok<AGROW>
    end
end

function [sameVal, diffVal] = R3_to_same_diff_(R3)
    % SAME = mean(diagonal)
    sameVal = mean(diag(R3), 'omitnan');
    % DIFF = mean(unique off-diagonals (upper triangle))
    off = [R3(1,2), R3(1,3), R3(2,3)];
    diffVal = mean(off, 'omitnan');
end

function R3 = sessionR3_splitBoot_balanced_classified_(beh, S, cell_type, sIdx, cats, eventName, winEvt, winBase, MIN_RT_MS, minTrialsPerCond, nBoot)
% Returns a 3x3 similarity matrix for a single session using:
%   - correct-only + valid
%   - balanced L/C/R trials (subsample without replacement to nMin)
%   - ONLY CLASSIFIED units (1/2/3)
%   - split-bootstrap correlation (two independent resamples per cue type)
%
% If unusable, returns [].

    R3 = [];

    session = beh(sIdx);
    if ~isfield(session,'trials') || isempty(session.trials), return; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), return; end

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

    % Balance L/C/R
    nMin = min([numel(idxL_all), numel(idxC_all), numel(idxR_all)]);
    if nMin < minTrialsPerCond, return; end

    idxL = idxL_all(randperm(numel(idxL_all), nMin));
    idxC = idxC_all(randperm(numel(idxC_all), nMin));
    idxR = idxR_all(randperm(numel(idxR_all), nMin));

    % Collect per-unit zEvt across kept trials (ONLY classified units)
    unitZ = {};
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

    % Split-bootstrap R3
    Rstack = nan(3,3,nBoot);

    for b = 1:nBoot
        rL1 = idxL(randi(numel(idxL), numel(idxL), 1));
        rC1 = idxC(randi(numel(idxC), numel(idxC), 1));
        rR1 = idxR(randi(numel(idxR), numel(idxR), 1));

        rL2 = idxL(randi(numel(idxL), numel(idxL), 1));
        rC2 = idxC(randi(numel(idxC), numel(idxC), 1));
        rR2 = idxR(randi(numel(idxR), numel(idxR), 1));

        A = nan(nUnits,3);
        B = nan(nUnits,3);

        for k = 1:nUnits
            z = unitZ{k};
            A(k,1) = mean(z(rL1), 'omitnan');
            A(k,2) = mean(z(rC1), 'omitnan');
            A(k,3) = mean(z(rR1), 'omitnan');

            B(k,1) = mean(z(rL2), 'omitnan');
            B(k,2) = mean(z(rC2), 'omitnan');
            B(k,3) = mean(z(rR2), 'omitnan');
        end

        if any(~isfinite(A(:))) || any(~isfinite(B(:))), continue; end

        Rb = corr(A, B, 'Rows','pairwise'); % 3x3
        if all(isfinite(Rb(:)))
            Rb = 0.5*(Rb + Rb'); % symmetric
            Rstack(:,:,b) = Rb;
        end
    end

    if all(all(all(~isfinite(Rstack))))
        return;
    end

    R3 = mean(Rstack, 3, 'omitnan');
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

function ct = safe_get_celltype_(cell_type, sIdx, ch, u)
    ct = [];
    if sIdx <= numel(cell_type) && ~isempty(cell_type{sIdx}) && iscell(cell_type{sIdx}) && ...
       ch   <= numel(cell_type{sIdx}) && ~isempty(cell_type{sIdx}{ch}) && iscell(cell_type{sIdx}{ch}) && ...
       u    <= numel(cell_type{sIdx}{ch})
        ct = cell_type{sIdx}{ch}{u};
    end
end