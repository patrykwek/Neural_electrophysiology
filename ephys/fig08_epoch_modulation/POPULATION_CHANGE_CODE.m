%% ========= FIG_SESSION_ACTIVITY_CHANGES__MI_AND_QUALITY__GLOBALWINDOWS.m =========
% Per-session analysis of activity changes using GLOBAL window intervals
% (applied exactly like your previous code: per-trial windows relative to cue/press/lick
% on a common cue-aligned SDF axis, no warping).
%
% Produces 3 publication-style panels:
% A) Per-session median |MI| across units (3 epoch lines; optional %correct overlay if detectable)
% B) Fraction of significantly modulated units per session (3 epoch lines + binomial Wilson 95% CI)
% C) Unit count per session + mean baseline firing rate (sanity checks)
%
% Notes:
% - MI is computed per trial as (FR_epoch - FR_baseline) / (FR_epoch + FR_baseline + eps)
% - Baseline window is fixed relative to cue: [-500 -100] ms
% - Per-unit MI summary = median(|MI_trial|) across trials
% - Per-session MI summary = median(per-unit median(|MI|)) across units
% - Significance per unit: signrank(MI_trial, 0), alpha=0.05 (two-sided)

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE';
outPng  = fullfile(baseOut, 'SESSION_MI_ACTIVITY_CHANGES.png');

% >>> REQUIRED: load classification output <<<
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';
assert(exist(cellTypeFile,'file')==2, 'Missing cell type file: %s', cellTypeFile);
CT = load(cellTypeFile);
assert(isfield(CT,'cell_type') && iscell(CT.cell_type), 'cell_type missing/wrong type in: %s', cellTypeFile);
cell_type = CT.cell_type;  % cell_type{sess}{ch}{u}: 0=Uncl, 1=MSN, 2=FSI, 3=TAN; []=no spikes (skipped)

% Analysis window used to BUILD SDFs (kept, but we will use an expanded fixed window below)
preCueMs   = 500;
postLickMs = 3000;

% FIG2-like trial selection window (relative to cue)
fixedWin = [-1000, 6000];

% Trial filters
MIN_RT_MS    = 100;
requireValid = true;

% PER-UNIT inclusion
minTrialsPerUnit = 10;

% SDF settings
dt_ms        = 10;
gaussSigmaMs = 25;

% Learning transition session (from Weibull)
sStar = 8;

% Baseline window for MI (relative to cue)
baselineWin_relCue = [-500, -100];  % ms

% Plot style
figPos  = [120, 120, 1200, 850];
fsTitle = 18;
fsAx    = 14;
lwAx    = 2.0;
lwLine  = 2.5;

% Line colors
colCue   = [0.15 0.50 0.80];
colPress = [0.85 0.35 0.15];
colLick  = [0.20 0.65 0.25];

rng(0);

%% ---- LOAD BEHAVIOR/SPIKES ----
assert(exist(matFile,'file')==2, 'File not found: %s', matFile);
S   = load(matFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'No behavior struct found.');
nSessions = numel(beh);

%% ---- PRECOMPUTE GAUSS KERNEL (unit area; Hz after dividing by dt) ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- LOAD GLOBAL WINDOWS (ROBUST PATH SEARCH) ----
candidates = { ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile('/Volumes/WD_BLACK/A_THESIS_FINAL','GENERAL_CHANGE_EARLY_LATE','GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(baseOut,'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat'), ...
    fullfile(fileparts(outPng),'GLM','GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat') ...
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

% Windows are stored as [start end] in ms RELATIVE TO EACH EVENT (lag axis).
wCue_rel   = GW.cue_ms;    % [ms, ms] relative to cue
wPress_rel = GW.press_ms;  % [ms, ms] relative to press
wLick_rel  = GW.lick_ms;   % [ms, ms] relative to lick

fprintf('\nLoaded GLOBAL windows from:\n  %s\n', winFile);
fprintf('Cue window (rel cue):     [%d, %d] ms\n', round(wCue_rel(1)),   round(wCue_rel(2)));
fprintf('Press window (rel press): [%d, %d] ms\n', round(wPress_rel(1)), round(wPress_rel(2)));
fprintf('Lick window (rel lick):   [%d, %d] ms\n', round(wLick_rel(1)),  round(wLick_rel(2)));

%% ---- BUILD A SINGLE CUE-ALIGNED TIME AXIS THAT GUARANTEES WINDOW COVERAGE ----
minLag = min([wCue_rel(1), wPress_rel(1), wLick_rel(1), baselineWin_relCue(1)]);
maxLag = max([wCue_rel(2), wPress_rel(2), wLick_rel(2), baselineWin_relCue(2)]);

winLeft  = fixedWin(1) + minLag;
winRight = fixedWin(2) + maxLag;

winLeft  = min(winLeft, -preCueMs);
winRight = max(winRight, fixedWin(2));

nT = numel(winLeft:dt_ms:winRight);
fprintf('\nSDF axis: [%d, %d] ms (dt=%d ms, nT=%d)\n', round(winLeft), round(winRight), dt_ms, nT);

%% ---- COMPUTE PER-SESSION METRICS ----
% Outputs (per session):
% medAbsMI(sess,3)   : median across units of per-unit median(|MI_trial|)
% fracSig(sess,3)    : fraction significant units (signrank MI_trial vs 0)
% nUnits(sess)       : number of included units
% baseFR(sess)       : mean across units of per-unit mean baseline FR (Hz)
% perfPct(sess)      : optional, percent-correct estimate if detectable, else NaN

[medAbsMI, fracSig, fracSig_CIlo, fracSig_CIhi, nUnits, baseFR, perfPct] = ...
    computeSessionMIandQuality_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue);

%% ---- SMOOTHING (ADDITION ONLY) ----
SMOOTH_WIN = 5;  % sessions (odd recommended)
fracSig_s = fracSig;
fracSig_s(:,1) = smoothdata(fracSig(:,1), 'movmean', SMOOTH_WIN);
fracSig_s(:,2) = smoothdata(fracSig(:,2), 'movmean', SMOOTH_WIN);
fracSig_s(:,3) = smoothdata(fracSig(:,3), 'movmean', SMOOTH_WIN);

medAbsMI_s = medAbsMI;
medAbsMI_s(:,1) = smoothdata(medAbsMI(:,1), 'movmean', SMOOTH_WIN);
medAbsMI_s(:,2) = smoothdata(medAbsMI(:,2), 'movmean', SMOOTH_WIN);
medAbsMI_s(:,3) = smoothdata(medAbsMI(:,3), 'movmean', SMOOTH_WIN);

baseFR_s = baseFR;
baseFR_s(:) = smoothdata(baseFR(:), 'movmean', SMOOTH_WIN);

%% ---- PLOT ----
x = (1:nSessions)';

epochLabels = {'Cue','Lever press','Reward lick'};

% ---------- Plot 1: % units modulated by the event ----------
figA = figure('Color','w','Position',figPos); axA = axes(figA); hold(axA,'on');

plot(axA, x, 100*fracSig_s(:,1), '-', 'LineWidth', lwLine, 'Color', colCue);
plot(axA, x, 100*fracSig_s(:,2), '-', 'LineWidth', lwLine, 'Color', colPress);
plot(axA, x, 100*fracSig_s(:,3), '-', 'LineWidth', lwLine, 'Color', colLick);

xline(axA, sStar+0.5, '--', 'LineWidth', 1.8);

ylabel(axA, '% units modulated');
xlabel(axA, 'Session');
title(axA, 'Event modulation prevalence', 'FontSize', fsTitle, 'FontWeight','bold');

set(axA, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
box(axA,'off');
xlim(axA, [1 nSessions]);
ylim(axA, [0 100]);

legA = legend(axA, epochLabels, 'Location','northwest');
set(legA, 'Box','off');

saveas(figA, outPng);  % keep same output name exactly as provided
fprintf('\nSaved: %s\n', outPng);

% ---------- Plot 2: Per-session median |MI| across units ----------
figB = figure('Color','w','Position',figPos); axB = axes(figB); hold(axB,'on');

plot(axB, x, medAbsMI_s(:,1), '-', 'LineWidth', lwLine, 'Color', colCue);
plot(axB, x, medAbsMI_s(:,2), '-', 'LineWidth', lwLine, 'Color', colPress);
plot(axB, x, medAbsMI_s(:,3), '-', 'LineWidth', lwLine, 'Color', colLick);

xline(axB, sStar+0.5, '--', 'LineWidth', 1.8);

ylabel(axB, 'Median |MI|');
xlabel(axB, 'Session');
title(axB, 'Per-session modulation magnitude', 'FontSize', fsTitle, 'FontWeight','bold');

set(axB, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
box(axB,'off');
xlim(axB, [1 nSessions]);

legB = legend(axB, epochLabels, 'Location','northwest');
set(legB, 'Box','off');

% Optional overlay of performance (if detectable) -- unchanged logic, now per-figure
if any(isfinite(perfPct))
    yyaxis(axB, 'right');
    plot(axB, x, perfPct, '-', 'LineWidth', 1.8);
    ylabel(axB, '% correct');
    set(axB, 'YColor', 'k');
    yyaxis(axB, 'left');
end

% ---------- Plot 3: Fraction significant with Wilson CI ----------
figC = figure('Color','w','Position',figPos); axC = axes(figC); hold(axC,'on');

shadeCI_(axC, x, fracSig(:,1), fracSig_CIlo(:,1), fracSig_CIhi(:,1));
shadeCI_(axC, x, fracSig(:,2), fracSig_CIlo(:,2), fracSig_CIhi(:,2));
shadeCI_(axC, x, fracSig(:,3), fracSig_CIlo(:,3), fracSig_CIhi(:,3));

plot(axC, x, fracSig_s(:,1), '-', 'LineWidth', lwLine, 'Color', colCue);
plot(axC, x, fracSig_s(:,2), '-', 'LineWidth', lwLine, 'Color', colPress);
plot(axC, x, fracSig_s(:,3), '-', 'LineWidth', lwLine, 'Color', colLick);

xline(axC, sStar+0.5, '--', 'LineWidth', 1.8);

ylabel(axC, 'Fraction significant');
xlabel(axC, 'Session');
title(axC, 'Significantly modulated units', 'FontSize', fsTitle, 'FontWeight','bold');

set(axC, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
box(axC,'off');
xlim(axC, [1 nSessions]);
ylim(axC, [0 1]);

legC = legend(axC, epochLabels, 'Location','northwest');
set(legC, 'Box','off');

% ---------- Plot 4: Unit count + baseline FR ----------
figD = figure('Color','w','Position',figPos); axD = axes(figD); hold(axD,'on');

yyaxis(axD, 'left');
plot(axD, x, nUnits, '-', 'LineWidth', lwLine);
ylabel(axD, 'Units/session');

yyaxis(axD, 'right');
plot(axD, x, baseFR_s, '-', 'LineWidth', lwLine);
ylabel(axD, 'Baseline firing rate (Hz)');

xline(axD, sStar+0.5, '--', 'LineWidth', 1.8);

xlabel(axD, 'Session');
title(axD, 'Recording yield and baseline rate', 'FontSize', fsTitle, 'FontWeight','bold');

set(axD, 'FontSize', fsAx, 'LineWidth', lwAx, 'TickDir','out');
box(axD,'off');
xlim(axD, [1 nSessions]);

%% ================= LOCAL HELPERS =================

function [medAbsMI, fracSig, CIlo, CIhi, nUnits, baseFR, perfPct] = ...
    computeSessionMIandQuality_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue)

    nSessions = numel(beh);

    medAbsMI = nan(nSessions,3);
    fracSig  = nan(nSessions,3);
    CIlo     = nan(nSessions,3);
    CIhi     = nan(nSessions,3);
    nUnits   = nan(nSessions,1);
    baseFR   = nan(nSessions,1);
    perfPct  = nan(nSessions,1);

    alpha = 0.05;

    for sIdx = 1:nSessions
        session = beh(sIdx);
        if ~isGoodTrialsStruct_(session)
            continue;
        end
        trials = session.trials;

        spikes_by_ch = getSpikesForSession_(session, S, sIdx);
        if isempty(spikes_by_ch) || ~iscell(spikes_by_ch)
            continue;
        end

        % ---- Behavior-only keepTrial (same logic as your code) ----
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
        if isempty(idxKeep)
            continue;
        end

        cueAbsK   = cueAbs(idxKeep);
        pressLatK = pressLat(idxKeep);
        lickLatK  = lickLat(idxKeep);

        % Optional: attempt to compute %correct if there is an obvious field
        perfPct(sIdx) = tryComputePerfPct_(trials(idxKeep));

        % Per-unit summaries
        unitMedAbsMI = nan(0,3);
        unitSig      = false(0,3);
        unitBaseFR   = nan(0,1);

        for ch = 1:numel(spikes_by_ch)
            uc = spikes_by_ch{ch};
            if ~iscell(uc) || isempty(uc), continue; end

            for u = 1:numel(uc)
                % only classified units
                ct = safe_get_celltype_(cell_type, sIdx, ch, u);
                if isempty(ct) || ~isnumeric(ct) || ~isscalar(ct) || ~ismember(ct, [1 2 3])
                    continue;
                end

                spk_abs = uc{u};
                if isempty(spk_abs) || ~isnumeric(spk_abs), continue; end
                spk_abs = double(spk_abs(:));

                % Active-trial filter: >=1 spike in [cue+fixedWin]
                activeTrials = false(numel(idxKeep),1);
                spk = spk_abs;
                if ~issorted(spk), spk = sort(spk); end
                j = 1;
                nSpk = numel(spk);

                for iTr = 1:numel(idxKeep)
                    cue0 = cueAbsK(iTr);
                    w0 = cue0 + fixedWin(1);
                    w1 = cue0 + fixedWin(2);

                    while j <= nSpk && spk(j) < w0
                        j = j + 1;
                    end
                    if j <= nSpk && spk(j) <= w1
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

                MI_cue_tr   = nan(nTrU,1);
                MI_press_tr = nan(nTrU,1);
                MI_lick_tr  = nan(nTrU,1);
                baseFR_tr   = nan(nTrU,1);

                for iTr = 1:nTrU
                    cue0   = cueU(iTr);
                    tPress = rtU(iTr);
                    tLick  = lkU(iTr);

                    % SDF (Hz)
                    spk_rel = spk_abs - cue0;
                    spk_rel = spk_rel(spk_rel >= winLeft & spk_rel <= winRight);

                    counts = zeros(1, nT);
                    if ~isempty(spk_rel)
                        jj = round((spk_rel - winLeft)/dt_ms) + 1;
                        jj = jj(jj >= 1 & jj <= nT);
                        counts = accumarray(jj(:), 1, [nT 1], @sum, 0).';
                    end

                    sm = conv(counts, g, 'same');
                    y  = sm / (dt_ms/1000);

                    % Indices for baseline (relative cue) and epochs (relative events)
                    idxBase  = msWindowToIdx_(baselineWin_relCue, winLeft, dt_ms, nT);

                    idxCue   = msWindowToIdx_([0      + wCue_rel(1),   0      + wCue_rel(2)],   winLeft, dt_ms, nT);
                    idxPress = msWindowToIdx_([tPress + wPress_rel(1), tPress + wPress_rel(2)], winLeft, dt_ms, nT);
                    idxLick  = msWindowToIdx_([tLick  + wLick_rel(1),  tLick  + wLick_rel(2)],  winLeft, dt_ms, nT);

                    frBase  = mean(y(idxBase(1):idxBase(2)),   'omitnan');
                    frCue   = mean(y(idxCue(1):idxCue(2)),     'omitnan');
                    frPress = mean(y(idxPress(1):idxPress(2)), 'omitnan');
                    frLick  = mean(y(idxLick(1):idxLick(2)),   'omitnan');

                    baseFR_tr(iTr) = frBase;

                    % MI per trial (signed), baseline-referenced
                    MI_cue_tr(iTr)   = (frCue   - frBase) ./ (frCue   + frBase + eps);
                    MI_press_tr(iTr) = (frPress - frBase) ./ (frPress + frBase + eps);
                    MI_lick_tr(iTr)  = (frLick  - frBase) ./ (frLick  + frBase + eps);
                end

                % Per-unit summaries
                mCue   = median(abs(MI_cue_tr),   'omitnan');
                mPress = median(abs(MI_press_tr), 'omitnan');
                mLick  = median(abs(MI_lick_tr),  'omitnan');

                % Significance per unit
                sigCue   = false;
                sigPress = false;
                sigLick  = false;

                if nnz(isfinite(MI_cue_tr)) >= 3
                    sigCue = (signrank(MI_cue_tr(isfinite(MI_cue_tr)), 0, 'alpha', alpha) < alpha);
                end
                if nnz(isfinite(MI_press_tr)) >= 3
                    sigPress = (signrank(MI_press_tr(isfinite(MI_press_tr)), 0, 'alpha', alpha) < alpha);
                end
                if nnz(isfinite(MI_lick_tr)) >= 3
                    sigLick = (signrank(MI_lick_tr(isfinite(MI_lick_tr)), 0, 'alpha', alpha) < alpha);
                end

                unitMedAbsMI(end+1,:) = [mCue mPress mLick]; %#ok<AGROW>
                unitSig(end+1,:)      = [sigCue sigPress sigLick]; %#ok<AGROW>
                unitBaseFR(end+1,1)   = mean(baseFR_tr, 'omitnan'); %#ok<AGROW>
            end
        end

        if isempty(unitMedAbsMI)
            continue;
        end

        nU = size(unitMedAbsMI,1);
        nUnits(sIdx) = nU;

        medAbsMI(sIdx,:) = median(unitMedAbsMI, 1, 'omitnan');
        baseFR(sIdx)     = mean(unitBaseFR, 'omitnan');

        for e = 1:3
            k = nnz(unitSig(:,e));
            fracSig(sIdx,e) = k / nU;

            [lo, hi] = wilsonCI_(k, nU, 0.05);
            CIlo(sIdx,e) = lo;
            CIhi(sIdx,e) = hi;
        end
    end
end

function [lo, hi] = wilsonCI_(k, n, alpha)
    if n <= 0
        lo = nan; hi = nan; return;
    end
    z = norminv(1 - alpha/2);
    phat = k / n;
    denom = 1 + (z^2)/n;
    center = (phat + (z^2)/(2*n)) / denom;
    half = (z/denom) * sqrt( (phat*(1-phat)/n) + (z^2)/(4*n^2) );
    lo = max(0, center - half);
    hi = min(1, center + half);
end

function shadeCI_(ax, x, y, lo, hi)
    % Minimal CI shading (no explicit color control beyond default face alpha)
    ok = isfinite(x) & isfinite(lo) & isfinite(hi);
    if nnz(ok) < 2, return; end
    xx = x(ok);
    ll = lo(ok);
    hh = hi(ok);
    patch(ax, [xx; flipud(xx)], [ll; flipud(hh)], 1, ...
        'FaceAlpha', 0.08, 'EdgeColor', 'none');
end

function pct = tryComputePerfPct_(trials)
    % Optional best-effort: returns %correct if there is a clearly named boolean field.
    pct = nan;
    if isempty(trials) || ~isstruct(trials), return; end
    candidates = {'correct','hit','success','rewarded','choice_correct','isCorrect'};
    for i = 1:numel(candidates)
        f = candidates{i};
        if isfield(trials, f)
            v = [trials.(f)];
            if isnumeric(v) || islogical(v)
                v = double(v(:));
                v = v(isfinite(v));
                if ~isempty(v)
                    pct = 100 * mean(v > 0);
                    return;
                end
            end
        end
    end
end

function idx12 = msWindowToIdx_(w_ms, winLeft_ms, dt_ms, nT)
    if any(~isfinite(w_ms))
        idx12 = [1 1];
        return;
    end
    i1 = round((w_ms(1) - winLeft_ms)/dt_ms) + 1;
    i2 = round((w_ms(2) - winLeft_ms)/dt_ms) + 1;
    i1 = max(1, min(nT, i1));
    i2 = max(1, min(nT, i2));
    if i2 < i1, tmp = i1; i1 = i2; i2 = tmp; end
    idx12 = [i1 i2];
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