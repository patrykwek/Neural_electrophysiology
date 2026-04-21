%% ========= FIG4C_LIKE_AVG_TRIAL_FR__EARLY_LATE__CLASSIFIEDONLY.m =========
% Replicates the "Average trial FR" style of Fig 4c:
%   LEFT: distribution of per-unit Average Trial FR (log x-axis), Early vs Late overlay
%   RIGHT: bar plot of group mean FR with error bars
%
% Early: sessions 1–8
% Late : sessions 29–36
%
% Uses your pipeline conventions:
%   - keep trials: valid, RT>=MIN_RT (no press/lick fixedWin restriction)
%   - per-unit inclusion: >= minTrialsPerUnit trials that contain >=1 spike in [cue+fixedWin(1), cue+lickLat]
%   - ONLY classified units (cell_type in {1,2,3})
%
% Per-unit "Average trial FR" definition:
%   For each kept trial:
%       w0 = cue + fixedWin(1)
%       w1 = cue + lickLat
%       FR_trial = spikes_in_window / ((w1-w0)/1000)
%   UnitAvgTrialFR = mean(FR_trial across kept trials)
%
% Output:
%   /Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_NEURAL_SIMILARITY
% Saves: MAT + PNG + FIG + SVG

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/2_FIRING_RATE_SIMILARITY';
if ~exist(outDir,'dir'), mkdir(outDir); end

% Early/Late sessions
stageRanges = [1 8; 29 36];

% Trial filters
fixedWin = [-1000, 6000];   % ms relative to cue (left edge); right edge uses lick latency
MIN_RT_MS = 100;
requireValid = true;

% Unit inclusion
minTrialsPerUnit = 10;

% Plot bins (log space)
nBins = 35;                 % similar feel to paper
frMin = 1e-2;               % Hz
frMax = 1e2;                % Hz (adjust if needed)
useProbability = false;     % paper shows "Count"; set true for probability

% Style (keep close to your figure style)
figPos = [150 150 1150 520];
titleFontSize = 26;
labelFontSize = 22;
tickFontSize  = 18;
axesTickLW    = 2.5;
axesTickLen   = [0.02 0.02];

% Colors (keep your plot style: early blue, late orange)
colEarly = [0.0 0.45 0.74];
colLate  = [0.85 0.33 0.10];

lineW = 3;

tag = sprintf('FIG4C_like_AvgTrialFR_Early%02d_%02d_Late%02d_%02d', ...
    stageRanges(1,1),stageRanges(1,2),stageRanges(2,1),stageRanges(2,2));

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
cell_type = CT.cell_type;

%% ---- SESSION GROUPS ----
e0 = max(1, min(nSessions, stageRanges(1,1)));
e1 = max(1, min(nSessions, stageRanges(1,2)));
l0 = max(1, min(nSessions, stageRanges(2,1)));
l1 = max(1, min(nSessions, stageRanges(2,2)));

earlySess = e0:e1;
lateSess  = l0:l1;

fprintf('Early sessions: %d-%d | Late sessions: %d-%d\n', e0,e1,l0,l1);

%% ---- COLLECT PER-UNIT Avg Trial FR for Early and Late ----
FR_early = [];
FR_late  = [];

nUnitsEarly = 0;
nUnitsLate  = 0;

for sIdx = [earlySess, lateSess]

    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch) || ~iscell(spikes_by_ch), continue; end

    % ---- Keep trials: valid + RT >= MIN ----
    keepTrial = false(numel(trials),1);
    cueAbs    = nan(numel(trials),1);
    pressLat  = nan(numel(trials),1);
    lickLat   = nan(numel(trials),1);

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

        lr = lick - cue;

        keepTrial(k) = true;
        cueAbs(k)    = cue;
        pressLat(k)  = rt; %#ok<NASGU>
        lickLat(k)   = lr;
    end

    idxKeep = find(keepTrial);
    if isempty(idxKeep), continue; end

    cueAbsK  = cueAbs(idxKeep);
    lickLatK = lickLat(idxKeep);

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc) || isempty(uc), continue; end

        for u = 1:numel(uc)

            % ONLY CLASSIFIED UNITS (1/2/3)
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

            % Trial-wise FR values for this unit
            fr_trials = nan(numel(idxKeep),1);
            activeMask = false(numel(idxKeep),1);

            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                lr0  = lickLatK(iTr);

                w0 = cue0 + fixedWin(1);
                w1 = cue0 + lr0;

                if ~(isfinite(w0) && isfinite(w1) && w1 > w0)
                    continue;
                end

                nsp = nnz(spk_abs >= w0 & spk_abs <= w1);
                dur_s = (w1 - w0)/1000;

                fr_trials(iTr) = nsp / dur_s;
                if nsp > 0
                    activeMask(iTr) = true;
                end
            end

            if nnz(activeMask) < minTrialsPerUnit
                continue;
            end

            fr_unit = mean(fr_trials(isfinite(fr_trials)), 'omitnan');
            if ~isfinite(fr_unit) || fr_unit <= 0
                continue;
            end

            if ismember(sIdx, earlySess)
                FR_early(end+1,1) = fr_unit; %#ok<AGROW>
                nUnitsEarly = nUnitsEarly + 1;
            else
                FR_late(end+1,1) = fr_unit; %#ok<AGROW>
                nUnitsLate = nUnitsLate + 1;
            end
        end
    end
end

assert(~isempty(FR_early) && ~isempty(FR_late), 'No units collected for early or late.');

fprintf('Collected units: Early=%d | Late=%d\n', nUnitsEarly, nUnitsLate);
fprintf('Median FR: Early=%.3f Hz | Late=%.3f Hz\n', median(FR_early), median(FR_late));

%% ---- SUMMARY STATS FOR BAR PLOT ----
% Use mean +/- SEM like typical bar plots; paper often uses mean.
mE = mean(FR_early);
mL = mean(FR_late);
semE = std(FR_early) / sqrt(numel(FR_early));
semL = std(FR_late) / sqrt(numel(FR_late));

% Optional: nonparametric p-value (Mann–Whitney) for reporting
p_ranksum = ranksum(FR_early, FR_late);

%% ---- FIGURE (2 panels) ----
fig = figure('Color','w','Position',figPos);

% --- LEFT: distribution (log x) ---
ax1 = subplot(1,2,1); hold(ax1,'on');

edges = logspace(log10(frMin), log10(frMax), nBins+1);

if useProbability
    normMode = 'probability';
    ylab = 'Probability';
else
    normMode = 'count';
    ylab = 'Count';
end

h1 = histogram(ax1, FR_early, edges, 'Normalization',normMode, 'DisplayStyle','stairs', 'LineWidth',lineW);
h2 = histogram(ax1, FR_late,  edges, 'Normalization',normMode, 'DisplayStyle','stairs', 'LineWidth',lineW);

% Set colors
h1.EdgeColor = colEarly;
h2.EdgeColor = colLate;

set(ax1, 'XScale','log');

xlabel(ax1, 'Average trial FR (Hz)', 'FontSize', labelFontSize);
ylabel(ax1, ylab, 'FontSize', labelFontSize);

title(ax1, 'Average trial FR', 'FontSize', titleFontSize, 'FontWeight','bold');

legend(ax1, {sprintf('Early (sessions %d-%d)', e0,e1), sprintf('Late (sessions %d-%d)', l0,l1)}, ...
    'Location','northwest');

set(ax1, 'FontSize', tickFontSize, ...
    'TickDir','out', 'LineWidth', axesTickLW, 'TickLength', axesTickLen, 'Box','off');

% --- RIGHT: bar plot mean +/- SEM ---
ax2 = subplot(1,2,2); hold(ax2,'on');

x = [1 2];
bar(ax2, 1, mE, 0.75, 'FaceColor', colEarly, 'EdgeColor','none');
bar(ax2, 2, mL, 0.75, 'FaceColor', colLate,  'EdgeColor','none');

errorbar(ax2, x, [mE mL], [semE semL], 'k', 'LineStyle','none', 'LineWidth',2.5, 'CapSize',14);

% Keep axis aesthetics similar to your style
xlim(ax2, [0.4 2.6]);
xticks(ax2, [1 2]);
xticklabels(ax2, {'Early','Late'});
ylabel(ax2, 'Average FR (Hz)', 'FontSize', labelFontSize);

title(ax2, sprintf('Mean \\pm SEM | ranksum p=%.3g', p_ranksum), ...
    'FontSize', titleFontSize, 'FontWeight','bold');

set(ax2, 'FontSize', tickFontSize, ...
    'TickDir','out', 'LineWidth', axesTickLW, 'TickLength', axesTickLen, 'Box','off');

%% ---- SAVE OUTPUTS ----
outMat = fullfile(outDir, [tag '_data.mat']);
save(outMat, 'FR_early','FR_late','earlySess','lateSess','minTrialsPerUnit','fixedWin','MIN_RT_MS','requireValid', ...
    'mE','mL','semE','semL','p_ranksum', '-v7.3');
fprintf('Saved MAT: %s\n', outMat);

outPng = fullfile(outDir, [tag '.png']);
saveas(fig, outPng);
fprintf('Saved PNG: %s\n', outPng);

outFig = fullfile(outDir, [tag '.fig']);
savefig(fig, outFig);
fprintf('Saved FIG: %s\n', outFig);

outSvg = fullfile(outDir, [tag '.svg']);
try
    print(fig, outSvg, '-dsvg');
    fprintf('Saved SVG: %s\n', outSvg);
catch ME
    warning('SVG save failed: %s', ME.message);
end

fprintf('DONE. Outputs in: %s\n', outDir);

%% ================= HELPERS =================

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