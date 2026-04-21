%% ========= FIG_SESSION_EPOCH_FIRINGRATE_HZ__GLOBALWINDOWS_STYLED__A_THESIS_FINAL.m =========
% Same computation as your Hz-per-session script (cue/press/lick epoch FR in Hz using GLOBAL windows),
% but renders the plot in the style of:
%   PLOT_PERCENT_CORRECT_WEIBULL_FIT_ASYMPTOTE_STYLED__A_THESIS_FINAL.m
%
% Plot:
%   - y-axis: Epoch firing rate (Hz)
%   - x-axis: Session
%   - 3 lines: Cue / Lever press / Reward lick
%   - Manual ticks: label every 5th session, custom tick marks and lowered labels
%   - (Optional) learning transition line at sStar+0.5
%
% Saves:
%   PNG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_EPOCH_FIRINGRATE_HZ_STYLED.png
%   SVG: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_EPOCH_FIRINGRATE_HZ_STYLED.svg
%   MAT: /Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE/SESSION_EPOCH_FIRINGRATE_HZ_STYLED.mat
%
% MAT contains:
%   medFR, medFR_s, nUnits, baseFR, perfPct, fixedWin, dt_ms, gaussSigmaMs, minTrialsPerUnit, sStar,
%   wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue, winLeft, winRight, nT, SMOOTH_WIN

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL/POPULATION_CHANGE';

outPng  = fullfile(baseOut, 'SESSION_EPOCH_FIRINGRATE_HZ_STYLED.png');
outSvg  = fullfile(baseOut, 'SESSION_EPOCH_FIRINGRATE_HZ_STYLED.svg');
outMat  = fullfile(baseOut, 'SESSION_EPOCH_FIRINGRATE_HZ_STYLED.mat');

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

% Baseline window (relative to cue) -- computed but not the plotted quantity
baselineWin_relCue = [-500, -100];  % ms

rng(0);

%% ---- STYLE (matched to Weibull-style script) ----
figPos = [120, 120, 900, 520];

dataLW      = 5;    % line width
dataMS      = 10;   % marker size
markerLW    = 2.5;  % marker edge width
eventLineLW = 5;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

axesTickLW  = 4.0;

%% ---- TICK / LABEL GEOMETRY ----
majorLenFrac      = 0.060;   % tick length as fraction of y-span
minorLenFrac      = 0.032;
tickLabelDownFrac = 0.0005;  % move tick numbers down (fraction of y-span)
xLabelDownFrac    = 0.08;    % move "Session" further down (fraction of y-span)

%% ---- COLORS (keep your original epoch colors) ----
colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

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

%% ---- COMPUTE PER-SESSION METRICS (Hz) ----
[medFR, nUnits, baseFR, perfPct] = ...
    computeSessionEpochFR_Hz_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue);

%% ---- SMOOTHING (ADDITION ONLY; matches your previous Hz script) ----
SMOOTH_WIN = 5;  % sessions (odd recommended)
medFR_s = medFR;
medFR_s(:,1) = smoothdata(medFR(:,1), 'movmean', SMOOTH_WIN);
medFR_s(:,2) = smoothdata(medFR(:,2), 'movmean', SMOOTH_WIN);
medFR_s(:,3) = smoothdata(medFR(:,3), 'movmean', SMOOTH_WIN);

%% ---- PLOT (STYLED) ----
xAll = (1:nSessions)';

fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

% y-limits based on data (clean, stable)
yAll = medFR_s(:);
yAll = yAll(isfinite(yAll));
if isempty(yAll)
    yl = [0 1];
else
    yMin = min(yAll);
    yMax = max(yAll);
    pad  = 0.08 * max(eps, (yMax - yMin));
    yl   = [max(0, yMin - pad), yMax + pad];
    if yl(2) <= yl(1), yl = [max(0, yl(1)-1), yl(1)+1]; end
end

set(ax, 'XLim',[1 nSessions], 'YLim',yl);

% -------------------- BACKGROUND SHADING --------------------
% early (yellowish): [1, sStar+0.5]
% late  (greenish):  [sStar+0.5, nSessions]
yl_patch = yl;

patch(ax, [1 (sStar+0.5) (sStar+0.5) 1], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
    [1 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % yellowish

patch(ax, [(sStar+0.5) nSessions nSessions (sStar+0.5)], [yl_patch(1) yl_patch(1) yl_patch(2) yl_patch(2)], ...
    [0 1 0], 'FaceAlpha', 0.10, 'EdgeColor','none');   % greenish
% ------------------------------------------------------------

% --- Lines (thick, styled) ---
plot(ax, xAll, medFR_s(:,1), '-', 'LineWidth', dataLW, 'Color', colCue);
plot(ax, xAll, medFR_s(:,2), '-', 'LineWidth', dataLW, 'Color', colPress);
plot(ax, xAll, medFR_s(:,3), '-', 'LineWidth', dataLW, 'Color', colLick);

% --- Open circles (markers) ---
plot(ax, xAll, medFR_s(:,1), 'o', 'LineStyle','none', ...
    'MarkerSize', dataMS, 'MarkerFaceColor', colCue, 'MarkerEdgeColor', colCue, 'LineWidth', markerLW);
plot(ax, xAll, medFR_s(:,2), 'o', 'LineStyle','none', ...
    'MarkerSize', dataMS, 'MarkerFaceColor', colPress, 'MarkerEdgeColor', colPress, 'LineWidth', markerLW);
plot(ax, xAll, medFR_s(:,3), 'o', 'LineStyle','none', ...
    'MarkerSize', dataMS, 'MarkerFaceColor', colLick, 'MarkerEdgeColor', colLick, 'LineWidth', markerLW);

% --- Event line (learning transition) ---
xline(ax, sStar+0.5, '--', 'LineWidth', eventLineLW, 'Color', [0.2 0.2 0.2]);

hxlab = xlabel(ax, 'Session', 'FontSize', labelFontSize);
ylabel(ax, 'Epoch firing rate (Hz)', 'FontSize', labelFontSize);
title(ax, 'Per-session epoch firing rate (GLOBAL windows)', ...
    'FontWeight','bold', 'FontSize', titleFontSize);

set(ax, 'FontSize', tickFontSize, 'Box','off', 'XLim',[1 nSessions], 'YLim',yl);
set(ax, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

%% ---- TICKS: label every 5th; draw ticks + labels manually (Weibull style) ----
ax.XTick = 1:nSessions;
labIdx = 5:5:nSessions;

% blank built-in labels; we draw them manually
ax.XTickLabel = repmat({''}, 1, nSessions);
xtickangle(ax, 0);

% Hide axis' own tick marks
set(ax, 'TickLength', [0 0]);

% Manual ticks (DATA units)
y0 = ax.YLim(1);
yr = diff(ax.YLim);

majorLen = majorLenFrac * yr;
minorLen = minorLenFrac * yr;

tickCol = [0 0 0];

for s = 1:nSessions
    if ismember(s, labIdx)
        L  = majorLen;
        lw = axesTickLW;
    else
        L  = minorLen;
        lw = max(1.5, axesTickLW * 0.55);
    end
    line(ax, [s s], [y0, y0 - L], 'Color', tickCol, 'LineWidth', lw, 'Clipping','off');
end

% Manual tick labels (moved down)
yTickText = y0 - (majorLen + tickLabelDownFrac*yr);
for s = labIdx
    text(ax, s, yTickText, num2str(s), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize', tickFontSize, ...
        'Color', [0 0 0], ...
        'Clipping','off');
end

% ---- MOVE X-LABEL DOWN BELOW TICK NUMBERS ----
xCenter = mean(ax.XLim);
yXlab   = y0 - (majorLen + xLabelDownFrac*yr);
set(hxlab, 'Units','data', 'Position',[xCenter, yXlab, 0]);

% Give more bottom room so nothing clips
ax.Position = [0.12 0.22 0.84 0.70];

% Ensure limits stay fixed
set(ax, 'XLim',[1 nSessions], 'YLim',yl);

% (CHANGED) legend removed

%% ---- SAVE OUTPUTS ----
saveas(fig, outPng);
saveas(fig, outSvg);

save(outMat, ...
    'medFR','medFR_s','nUnits','baseFR','perfPct', ...
    'fixedWin','dt_ms','gaussSigmaMs','minTrialsPerUnit','sStar', ...
    'MIN_RT_MS','requireValid','baselineWin_relCue', ...
    'wCue_rel','wPress_rel','wLick_rel', ...
    'winLeft','winRight','nT','SMOOTH_WIN');

fprintf('\nSaved: %s\n', outPng);
fprintf('Saved: %s\n', outSvg);
fprintf('Saved: %s\n', outMat);

%% ================= LOCAL HELPERS =================

function [medFR, nUnits, baseFR, perfPct] = ...
    computeSessionEpochFR_Hz_(beh, S, cell_type, ...
        fixedWin, requireValid, MIN_RT_MS, minTrialsPerUnit, ...
        dt_ms, g, winLeft, winRight, nT, ...
        wCue_rel, wPress_rel, wLick_rel, baselineWin_relCue)

    nSessions = numel(beh);

    medFR   = nan(nSessions,3);
    nUnits  = nan(nSessions,1);
    baseFR  = nan(nSessions,1);
    perfPct = nan(nSessions,1);

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

        % ---- Behavior-only keepTrial (same logic as your MI code) ----
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
        unitMedFR  = nan(0,3);
        unitBaseFR = nan(0,1);

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

                FR_cue_tr   = nan(nTrU,1);
                FR_press_tr = nan(nTrU,1);
                FR_lick_tr  = nan(nTrU,1);
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

                    baseFR_tr(iTr)   = frBase;
                    FR_cue_tr(iTr)   = frCue;
                    FR_press_tr(iTr) = frPress;
                    FR_lick_tr(iTr)  = frLick;
                end

                % Per-unit summaries (Hz): median across trials
                mCue   = median(FR_cue_tr,   'omitnan');
                mPress = median(FR_press_tr, 'omitnan');
                mLick  = median(FR_lick_tr,  'omitnan');

                unitMedFR(end+1,:) = [mCue mPress mLick]; %#ok<AGROW>
                unitBaseFR(end+1,1) = mean(baseFR_tr, 'omitnan'); %#ok<AGROW>
            end
        end

        if isempty(unitMedFR)
            continue;
        end

        nU = size(unitMedFR,1);
        nUnits(sIdx) = nU;

        % Per-session summary: median across units
        medFR(sIdx,:) = median(unitMedFR, 1, 'omitnan');
        baseFR(sIdx)  = mean(unitBaseFR, 'omitnan');
    end
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