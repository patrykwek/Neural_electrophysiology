%% ========= ADDON_GLOBAL_WINDOW_VALIDATION_PLOTS__PLOT1_PLOT2_PLOT4A.m =========
% PURPOSE:
%   Produce the 3 “validity of chosen windows” plots AFTER you already ran the GLM script
%   that SAVED:
%       outDir/GLM_SESSIONLEVEL_OUTPUTS_for_window_validity.mat
%
% INPUT:
%   GLM_SESSIONLEVEL_OUTPUTS_for_window_validity.mat (contains GLM_OUTPUTS struct)
%
% OUTPUTS (saved into the same outDir):
%   PLOT1_pooled_kernels_with_global_windows.png/.svg/.fig
%   PLOT2_evidence_tscore_with_windows.png/.svg/.fig
%   PLOT4A_threshold_sweep_window_stability.png/.svg/.fig
%   GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat
%
% ADDED (requested):
%   PLOT_E_EarlyLate_GLM_R2_bar.png/.svg/.fig
%   - Early vs Late session-level GLM R^2 bar (mean±SEM) + grey dots
%   - Nonparametric permutation test on difference of means (Early-Late)
%
% ADDED (requested):
%   PLOT_E3_R2_vs_nUnits_scatter.png/.svg/.fig
%   - Scatter sessR2_full vs sess_nUnitsUsed
%   - Spearman rho + p printed
%
% --- MINIMAL FIX FOR R^2 (ADDED) ---
%   1) Prefer loading the *SIMPLIFIED* file if it exists (prevents old non-R^2 file being used).
%   2) Enforce that sessR2_full looks like real R^2 (not >1). If not, throw an error.
%
% --- OPTION 1 IMPLEMENTATION (ADDED) ---
%   Use session-level PSTH R^2 if present: GLM_OUTPUTS.sessR2_psth
%   Otherwise fall back to sessR2_full.
%   Output file names are NOT changed.
%
% --- MINIMAL ADDITION (ADDED) ---
%   Compute overall mean ± SEM across all finite session-level R^2 values
%   so you can report:
%   “Session-level fits were broadly positive overall (mean ± SEM: XXX; range: ... )”

clear; clc;

%% ---------------- USER: POINT TO SAVED OUTPUTS ----------------
outDir = '/Volumes/WD_BLACK/A_THESIS_FINAL/GLM';

% ---- MINIMAL FIX: Prefer the simplified file FIRST ----
glmSaveFile_simplified = fullfile(outDir, 'GLM_SESSIONLEVEL_OUTPUTS_for_window_validity_SIMPLIFIED.mat');
glmSaveFile_canonical  = fullfile(outDir, 'GLM_SESSIONLEVEL_OUTPUTS_for_window_validity.mat');

if exist(glmSaveFile_simplified,'file') == 2
    glmSaveFile = glmSaveFile_simplified;
elseif exist(glmSaveFile_canonical,'file') == 2
    glmSaveFile = glmSaveFile_canonical;
else
    error('Missing saved GLM outputs: %s (and %s)', glmSaveFile_canonical, glmSaveFile_simplified);
end

tmp = load(glmSaveFile);
assert(isfield(tmp,'GLM_OUTPUTS') && isstruct(tmp.GLM_OUTPUTS), 'GLM_OUTPUTS struct missing in: %s', glmSaveFile);
G = tmp.GLM_OUTPUTS;

% ---- Pull required variables from the saved struct ----
sessKernels_cue   = G.sessKernels_cue;
sessKernels_press = G.sessKernels_press;
sessKernels_lick  = G.sessKernels_lick;

lags_cue_ms   = G.lags_cue_ms;
lags_press_ms = G.lags_press_ms;
lags_lick_ms  = G.lags_lick_ms;

% ---- Choose metric for plots E and E3: prefer PSTH R^2 ----
if isfield(G,'sessR2_psth')
    sessR2_metric = G.sessR2_psth;
    metricNameForPrint = 'sessR2_psth';
    metricLabel = 'PSTH R^2';
    metricTitle = 'PSTH fit quality (session-level)';
else
    if isfield(G,'sessR2_full')
        sessR2_metric = G.sessR2_full;
    else
        sessR2_metric = nan(size(sessKernels_cue,1),1);
    end
    metricNameForPrint = 'sessR2_full';
    metricLabel = 'GLM R^2';
    metricTitle = 'GLM fit quality (session-level)';
end

if isfield(G,'sess_nUnitsUsed')
    sess_nUnitsUsed = G.sess_nUnitsUsed;
else
    sess_nUnitsUsed = nan(size(sessKernels_cue,1),1);
end

assert(exist(outDir,'dir')==7, 'outDir does not exist: %s', outDir);

% ---- Optional provenance used only for printing/record keeping ----
if isfield(G,'dt_ms'), dt_ms = G.dt_ms; else, dt_ms = nan; end %#ok<NASGU>

% ---- MINIMAL R^2 SANITY CHECK (keep original behavior for sessR2_full) ----
% If we're using sessR2_full, ensure it looks like R^2 (not >1).
fprintf('\nLoaded GLM outputs: %s\n', glmSaveFile);
fprintf('Using metric for Plot E/E3: %s\n', metricNameForPrint);

if strcmp(metricNameForPrint,'sessR2_full')
    finiteR2 = sessR2_metric(isfinite(sessR2_metric));
    if ~isempty(finiteR2)
        fprintf('sessR2_full sanity: min=%.6g, max=%.6g\n', min(finiteR2), max(finiteR2));
        if max(finiteR2) > 1.05
            error(['sessR2_full exceeds 1 (max=%.6g). This is NOT valid R^2. ', ...
                   'You are likely loading an older file or a different metric. ', ...
                   'Rerun the GLM script and ensure the addon loads *_SIMPLIFIED.mat.'], max(finiteR2));
        end
    end
else
    finiteR2 = sessR2_metric(isfinite(sessR2_metric));
    if ~isempty(finiteR2)
        fprintf('%s sanity: min=%.6g, max=%.6g\n', metricNameForPrint, min(finiteR2), max(finiteR2));
    end
end
% -----------------------------------------

%% ---------------- MINIMAL ADDITION: OVERALL SESSION-LEVEL SUMMARY ----------------
allR2 = sessR2_metric(isfinite(sessR2_metric));
if ~isempty(allR2)
    overallMeanR2 = mean(allR2,'omitnan');
    overallSemR2  = std(allR2,'omitnan') / sqrt(numel(allR2));
    overallMinR2  = min(allR2);
    overallMaxR2  = max(allR2);
    overallNR2    = numel(allR2);

    fprintf('\n=== Overall %s across all sessions ===\n', metricNameForPrint);
    fprintf('nSessions = %d\n', overallNR2);
    fprintf('mean = %.6f\n', overallMeanR2);
    fprintf('SEM  = %.6f\n', overallSemR2);
    fprintf('range = [%.6f, %.6f]\n', overallMinR2, overallMaxR2);
else
    overallMeanR2 = nan;
    overallSemR2  = nan;
    overallMinR2  = nan;
    overallMaxR2  = nan;
    overallNR2    = 0;

    fprintf('\n=== Overall %s across all sessions ===\n', metricNameForPrint);
    fprintf('No finite session-level values found.\n');
end
% -------------------------------------------------------------------------------

%% ---------------- SETTINGS FOR WINDOW DEFINITION ----------------
tau_main    = 1.5;    % LOWER threshold => wider windows (was 2.0)
minContigMs = 50;     % LOWER min duration => allows broader candidate segments (was 80)
postOnly    = false;  % ALLOW negative lags (was true)

cueRestrictAroundZero  = false;
lickRestrictAroundZero = false;

tauSweep    = [1.0 1.5 2.0 2.5 3.0];

% Plot styling (match your earlier figures)
titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 28;
axesTickLW    = 4.0;
traceLW       = 6;

% Colors (use your existing if present)
if ~exist('colCue','var'),   colCue   = [1 0 0]; end
if ~exist('colPress','var'), colPress = [0.10 0.55 0.95]; end
if ~exist('colLick','var'),  colLick  = [0.15 0.70 0.20]; end

%% ---------------- GATHER SESSIONS USED (POOLING) ----------------
Kcue_all   = sessKernels_cue;    Kcue_all   = Kcue_all(any(isfinite(Kcue_all),2),:);
Kpress_all = sessKernels_press;  Kpress_all = Kpress_all(any(isfinite(Kpress_all),2),:);
Klick_all  = sessKernels_lick;   Klick_all  = Klick_all(any(isfinite(Klick_all),2),:);

assert(~isempty(Kcue_all) && ~isempty(Kpress_all) && ~isempty(Klick_all), ...
    'No session kernels available for pooled window definition.');

nSessCue   = size(Kcue_all,1);
nSessPress = size(Kpress_all,1);
nSessLick  = size(Klick_all,1);

fprintf('Pooling sessions: Cue=%d, Press=%d, Lick=%d\n', nSessCue, nSessPress, nSessLick);

%% ---------------- COMPUTE POOLED MEAN + SEM (SESSION-LEVEL) ----------------
[muCue,   semCue]   = mean_sem_rows_(Kcue_all);
[muPress, semPress] = mean_sem_rows_(Kpress_all);
[muLick,  semLick]  = mean_sem_rows_(Klick_all);

tCue   = abs(muCue)   ./ (semCue   + eps);
tPress = abs(muPress) ./ (semPress + eps);
tLick  = abs(muLick)  ./ (semLick  + eps);

%% ---------------- DEFINE GLOBAL WINDOWS FROM POOLED KERNELS ----------------
[wCue,   idxCueWin]   = defineWindow_from_tscore_(tCue,   lags_cue_ms,   tau_main, minContigMs, postOnly, ...
    cueRestrictAroundZero, 0, inf, false);
[wLick,  idxLickWin]  = defineWindow_from_tscore_(tLick,  lags_lick_ms,  tau_main, minContigMs, postOnly, ...
    lickRestrictAroundZero, 0, inf, false);

[wPress, idxPressWin] = defineWindow_from_tscore_(tPress, lags_press_ms, tau_main, minContigMs, postOnly, ...
    false, 0, inf, false);

fprintf('\n=== GLOBAL windows @ tau=%.2f, minContig=%d ms (postOnly=%d) ===\n', tau_main, minContigMs, postOnly);
fprintf('Cue:   [%d, %d] ms\n', wCue(1),   wCue(2));
fprintf('Press: [%d, %d] ms\n', wPress(1), wPress(2));
fprintf('Lick:  [%d, %d] ms\n', wLick(1),  wLick(2));

%% ---------------- SAVE GLOBAL WINDOWS FOR NEXT STEPS ----------------
winFile = fullfile(outDir, 'GLOBAL_KERNEL_WINDOWS_from_pooled_tscore.mat');
GLOBAL_WINDOWS = struct();
GLOBAL_WINDOWS.method = 'pooled_session_kernels__t=abs(mu)/SEM__best_contiguous_suprathreshold';
GLOBAL_WINDOWS.tau_main = tau_main;
GLOBAL_WINDOWS.minContigMs = minContigMs;
GLOBAL_WINDOWS.postOnly = postOnly;

GLOBAL_WINDOWS.cue_ms   = wCue;
GLOBAL_WINDOWS.press_ms = wPress;
GLOBAL_WINDOWS.lick_ms  = wLick;

GLOBAL_WINDOWS.idxCueWin   = idxCueWin;
GLOBAL_WINDOWS.idxPressWin = idxPressWin;
GLOBAL_WINDOWS.idxLickWin  = idxLickWin;

GLOBAL_WINDOWS.lags_cue_ms   = lags_cue_ms;
GLOBAL_WINDOWS.lags_press_ms = lags_press_ms;
GLOBAL_WINDOWS.lags_lick_ms  = lags_lick_ms;

GLOBAL_WINDOWS.source_file = glmSaveFile;
save(winFile, 'GLOBAL_WINDOWS');
fprintf('Saved global windows: %s\n', winFile);

%% ============================================================
%  PLOT 1: Pooled kernel mean ± SEM with chosen window overlaid
%% ============================================================
figP1 = figure('Color','w','Position',[120 120 1500 520]);

subplot(1,3,1);
plot_kernel_pooled_with_window_(gca, lags_cue_ms, muCue, semCue, idxCueWin, colCue, ...
    sprintf('Cue (pooled session kernel)  n=%d', nSessCue), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

subplot(1,3,2);
plot_kernel_pooled_with_window_(gca, lags_press_ms, muPress, semPress, idxPressWin, colPress, ...
    sprintf('Press (pooled session kernel)  n=%d', nSessPress), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

subplot(1,3,3);
plot_kernel_pooled_with_window_(gca, lags_lick_ms, muLick, semLick, idxLickWin, colLick, ...
    sprintf('Lick (pooled session kernel)  n=%d', nSessLick), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

sgtitle(sprintf('PLOT 1: Global windows from pooled session-kernels | t(l)=|\\mu|/SEM, \\tau=%.2f', tau_main), ...
    'FontWeight','bold','FontSize',titleFontSize);

outP1_png = fullfile(outDir, 'PLOT1_pooled_kernels_with_global_windows.png');
outP1_svg = fullfile(outDir, 'PLOT1_pooled_kernels_with_global_windows.svg');
outP1_fig = fullfile(outDir, 'PLOT1_pooled_kernels_with_global_windows.fig');
saveas(figP1, outP1_png);
saveas(figP1, outP1_svg);
savefig(figP1, outP1_fig);
close(figP1);
fprintf('Saved: %s\n', outP1_png);
fprintf('Saved: %s\n', outP1_svg);
fprintf('Saved: %s\n', outP1_fig);

%% ============================================================
%  PLOT 2: Evidence trace t(lag)=|mu|/SEM with threshold + window
%% ============================================================
figP2 = figure('Color','w','Position',[120 120 1500 520]);

subplot(1,3,1);
plot_tscore_with_window_(gca, lags_cue_ms, tCue, tau_main, idxCueWin, colCue, ...
    sprintf('Cue: t(l)=|\\mu|/SEM  n=%d', nSessCue), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

subplot(1,3,2);
plot_tscore_with_window_(gca, lags_press_ms, tPress, tau_main, idxPressWin, colPress, ...
    sprintf('Press: t(l)=|\\mu|/SEM  n=%d', nSessPress), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

subplot(1,3,3);
plot_tscore_with_window_(gca, lags_lick_ms, tLick, tau_main, idxLickWin, colLick, ...
    sprintf('Lick: t(l)=|\\mu|/SEM  n=%d', nSessLick), ...
    titleFontSize, labelFontSize, tickFontSize, axesTickLW, traceLW);

sgtitle(sprintf('PLOT 2: Window-definition evidence | threshold \\tau=%.2f, minContig=%d ms', tau_main, minContigMs), ...
    'FontWeight','bold','FontSize',titleFontSize);

outP2_png = fullfile(outDir, 'PLOT2_evidence_tscore_with_windows.png');
outP2_svg = fullfile(outDir, 'PLOT2_evidence_tscore_with_windows.svg');
outP2_fig = fullfile(outDir, 'PLOT2_evidence_tscore_with_windows.fig');
saveas(figP2, outP2_png);
saveas(figP2, outP2_svg);
savefig(figP2, outP2_fig);
close(figP2);
fprintf('Saved: %s\n', outP2_png);
fprintf('Saved: %s\n', outP2_svg);
fprintf('Saved: %s\n', outP2_fig);

%% ============================================================
%  PLOT 4A: Threshold sweep stability of window start/end
%% ============================================================
winCueSweep   = nan(numel(tauSweep),2);
winPressSweep = nan(numel(tauSweep),2);
winLickSweep  = nan(numel(tauSweep),2);

for i=1:numel(tauSweep)
    tau = tauSweep(i);

    [w,~] = defineWindow_from_tscore_(tCue,   lags_cue_ms,   tau, minContigMs, postOnly, ...
        cueRestrictAroundZero, 0, inf, false);
    winCueSweep(i,:) = w;

    [w,~] = defineWindow_from_tscore_(tPress, lags_press_ms, tau, minContigMs, postOnly, ...
        false, 0, inf, false);
    winPressSweep(i,:) = w;

    [w,~] = defineWindow_from_tscore_(tLick,  lags_lick_ms,  tau, minContigMs, postOnly, ...
        lickRestrictAroundZero, 0, inf, false);
    winLickSweep(i,:) = w;
end

figP4A = figure('Color','w','Position',[120 120 1500 520]);

subplot(1,3,1); hold on;
plot(tauSweep, winCueSweep(:,1),  '-o', 'LineWidth', traceLW);
plot(tauSweep, winCueSweep(:,2),  '-o', 'LineWidth', traceLW);
title('Cue window stability', 'FontWeight','bold','FontSize',titleFontSize);
xlabel('\tau', 'FontSize', labelFontSize);
ylabel('Window (ms)', 'FontSize', labelFontSize);
legend({'start','end'}, 'Location','best');
set(gca,'FontSize',tickFontSize,'LineWidth',axesTickLW,'TickDir','out','Box','off');

subplot(1,3,2); hold on;
plot(tauSweep, winPressSweep(:,1), '-o', 'LineWidth', traceLW);
plot(tauSweep, winPressSweep(:,2), '-o', 'LineWidth', traceLW);
title('Press window stability', 'FontWeight','bold','FontSize',titleFontSize);
xlabel('\tau', 'FontSize', labelFontSize);
ylabel('Window (ms)', 'FontSize', labelFontSize);
legend({'start','end'}, 'Location','best');
set(gca,'FontSize',tickFontSize,'LineWidth',axesTickLW,'TickDir','out','Box','off');

subplot(1,3,3); hold on;
plot(tauSweep, winLickSweep(:,1), '-o', 'LineWidth', traceLW);
plot(tauSweep, winLickSweep(:,2), '-o', 'LineWidth', traceLW);
title('Lick window stability', 'FontWeight','bold','FontSize',titleFontSize);
xlabel('\tau', 'FontSize', labelFontSize);
ylabel('Window (ms)', 'FontSize', labelFontSize);
legend({'start','end'}, 'Location','best');
set(gca,'FontSize',tickFontSize,'LineWidth',axesTickLW,'TickDir','out','Box','off');

sgtitle('PLOT 4A: Window stability under threshold sweep', ...
    'FontWeight','bold','FontSize',titleFontSize);

outP4A_png = fullfile(outDir, 'PLOT4A_threshold_sweep_window_stability.png');
outP4A_svg = fullfile(outDir, 'PLOT4A_threshold_sweep_window_stability.svg');
outP4A_fig = fullfile(outDir, 'PLOT4A_threshold_sweep_window_stability.fig');
saveas(figP4A, outP4A_png);
saveas(figP4A, outP4A_svg);
savefig(figP4A, outP4A_fig);
close(figP4A);
fprintf('Saved: %s\n', outP4A_png);
fprintf('Saved: %s\n', outP4A_svg);
fprintf('Saved: %s\n', outP4A_fig);

%% ============================================================
%  PLOT E (ADDED): Early vs Late session-level metric bar + dots
%   - Uses PSTH R^2 if available; otherwise GLM R^2
%   - Output file names unchanged
%% ============================================================
earlyGroup = 1:8;
lateGroup  = 29:36;

earlyUse = earlyGroup(earlyGroup>=1 & earlyGroup<=numel(sessR2_metric));
lateUse  = lateGroup(lateGroup>=1  & lateGroup<=numel(sessR2_metric));

earlyR2 = sessR2_metric(earlyUse); earlyR2 = earlyR2(isfinite(earlyR2));
lateR2  = sessR2_metric(lateUse);  lateR2  = lateR2(isfinite(lateR2));

colEarlyShade = [1 1 0]; % yellow
colLateShade  = [0 1 0]; % green
greyDot       = [0.6 0.6 0.6];
jit           = 0.10;

figE = figure('Color','w','Position',[120,120,520,700]);
axE = axes(figE); hold(axE,'on');

mE = mean(earlyR2,'omitnan');
mL = mean(lateR2,'omitnan');
semE = std(earlyR2,'omitnan') / sqrt(max(1,numel(earlyR2)));
semL = std(lateR2,'omitnan')  / sqrt(max(1,numel(lateR2)));

xPos = [1 2];
barW = 0.60;

b = bar(axE, xPos, [mE mL], barW, 'FaceColor','flat', 'EdgeColor',[0 0 0], 'LineWidth', axesTickLW);
b.CData(1,:) = colEarlyShade;
b.CData(2,:) = colLateShade;
b.FaceAlpha  = 0.25;

errorbar(axE, xPos, [mE mL], [semE semL], 'k', 'LineStyle','none', 'LineWidth', axesTickLW, 'CapSize', 18);

rng(0);
if ~isempty(earlyR2)
    plot(axE, 1 + (rand(size(earlyR2))-0.5)*2*jit, earlyR2, 'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end
if ~isempty(lateR2)
    plot(axE, 2 + (rand(size(lateR2))-0.5)*2*jit,  lateR2,  'o', 'MarkerSize', 8, ...
        'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);
end

set(axE, 'XLim',[0.5 2.5], 'XTick',[1 2], 'XTickLabel',{'Early','Late'});
ylabel(axE, metricLabel, 'FontSize', labelFontSize);
title(axE, metricTitle, 'FontWeight','bold', 'FontSize', titleFontSize);

set(axE, 'FontSize', tickFontSize, 'Box','off');
set(axE, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

% -------------------------------------------------------------------------
% CHANGE REQUESTED:
% Replace the two-sample label-shuffle test with the SAME permutation type
% used in your split-half script: a sign-flip (paired) permutation test.
%
% Here: treat Early/Late as paired by position within the group (1..N),
% and run a sign-flip test on the paired differences.
% -------------------------------------------------------------------------
if ~isempty(earlyR2) && ~isempty(lateR2)

    nPerm = 5000;

    nPair = min(numel(earlyR2), numel(lateR2));
    earlyP = earlyR2(1:nPair);
    lateP  = lateR2(1:nPair);

    diffs = lateP(:) - earlyP(:);          % Late - Early
    diffs = diffs(isfinite(diffs));

    Tobs = mean(diffs, 'omitnan');

    rng(0);
    Tperm = nan(nPerm,1);
    for k = 1:nPerm
        flip = (rand(size(diffs)) > 0.5) * 2 - 1; % +/- 1
        Tperm(k) = mean(diffs .* flip, 'omitnan');
    end

    p_perm = (1 + nnz(abs(Tperm) >= abs(Tobs))) / (nPerm + 1);

    addSigIfNeeded_(axE, 1, 2, p_perm, [mE mL], axesTickLW, tickFontSize);

    fprintf('\n=== %s Early vs Late (paired sign-flip permutation test) ===\n', metricNameForPrint);
    fprintf('Early n=%d | Late n=%d | paired n=%d\n', numel(earlyR2), numel(lateR2), nPair);
    fprintf('Observed mean diff (Late-Early) = %.6f\n', Tobs);
    fprintf('Permutation p = %.6g (nPerm=%d)\n', p_perm, nPerm);
end
% -------------------------------------------------------------------------

% Axis limits: allow negative, cap at 1 (PSTH R^2 can be negative; GLM R^2 can be negative too)
yAll = [earlyR2(:); lateR2(:); (mE+semE); (mL+semL)];
yAll = yAll(isfinite(yAll));
if isempty(yAll)
    ylim(axE, [-0.2 1]);
else
    yMin = min(yAll);
    yMax = max(yAll);
    pad  = 0.12 * max(eps, (yMax - yMin));
    yl   = [yMin - pad, min(1, yMax + pad)];
    if yl(2) <= yl(1), yl = [yl(1)-0.1, min(1, yl(1)+0.1)]; end
    ylim(axE, yl);
end

outE_png = fullfile(outDir, 'PLOT_E_EarlyLate_GLM_R2_bar.png');
outE_svg = fullfile(outDir, 'PLOT_E_EarlyLate_GLM_R2_bar.svg');
outE_fig = fullfile(outDir, 'PLOT_E_EarlyLate_GLM_R2_bar.fig');
saveas(figE, outE_png);
saveas(figE, outE_svg);
savefig(figE, outE_fig);
close(figE);
fprintf('Saved: %s\n', outE_png);
fprintf('Saved: %s\n', outE_svg);
fprintf('Saved: %s\n', outE_fig);

%% ============================================================
%  PLOT E3 (ADDED): metric vs sess_nUnitsUsed scatter (+ Spearman)
%   - Uses PSTH R^2 if available; otherwise GLM R^2
%   - Output file names unchanged
%% ============================================================
ok = isfinite(sessR2_metric) & isfinite(sess_nUnitsUsed);
xU = sess_nUnitsUsed(ok);
yR = sessR2_metric(ok);

figE3 = figure('Color','w','Position',[120,120,720,620]);
axE3 = axes(figE3); hold(axE3,'on');

plot(axE3, xU, yR, 'o', 'MarkerSize', 9, ...
    'MarkerFaceColor', greyDot, 'MarkerEdgeColor', greyDot, 'LineWidth', 1.0);

xlabel(axE3, 'Units used (n)', 'FontSize', labelFontSize);
ylabel(axE3, metricLabel, 'FontSize', labelFontSize);
title(axE3, sprintf('%s vs unit count (session-level)', metricTitle), 'FontWeight','bold', 'FontSize', titleFontSize);

set(axE3, 'FontSize', tickFontSize, 'Box','off');
set(axE3, 'TickDir','out', 'LineWidth', axesTickLW, 'Layer','top');

if ~isempty(xU) && ~isempty(yR)
    [rho, pSpearman] = corr(xU(:), yR(:), 'Type','Spearman', 'Rows','complete');
    fprintf('\n=== %s vs nUnitsUsed (Spearman) ===\n', metricNameForPrint);
    fprintf('nSessions = %d\n', numel(xU));
    fprintf('Spearman rho = %.4f, p = %.6g\n', rho, pSpearman);

    txt = sprintf('\\rho = %.2f, p = %.3g', rho, pSpearman);
    text(axE3, 0.02, 0.98, txt, 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', 'FontSize', 18);
end

if isempty(yR)
    ylim(axE3, [-0.2 1]);
else
    ylim(axE3, [min(yR)-0.05 min(1, max(yR)+0.05)]);
end

outE3_png = fullfile(outDir, 'PLOT_E3_R2_vs_nUnits_scatter.png');
outE3_svg = fullfile(outDir, 'PLOT_E3_R2_vs_nUnits_scatter.svg');
outE3_fig = fullfile(outDir, 'PLOT_E3_R2_vs_nUnits_scatter.fig');
saveas(figE3, outE3_png);
saveas(figE3, outE3_svg);
savefig(figE3, outE3_fig);
close(figE3);
fprintf('Saved: %s\n', outE3_png);
fprintf('Saved: %s\n', outE3_svg);
fprintf('Saved: %s\n', outE3_fig);

fprintf('\nDONE: validity plots produced.\n');

%% =============================== HELPERS ===============================

function [mu, sem] = mean_sem_rows_(K)
    mu  = mean(K, 1, 'omitnan');
    sem = std(K, 0, 1, 'omitnan') ./ sqrt(size(K,1) + eps);
end

function [w_ms, idxWin] = defineWindow_from_tscore_(tScore, lags_ms, tau, minContigMs, postOnly, ...
                                                   restrictAroundZero, centerMs, maxAbsLagMs, requireCrossZero)
    lags_ms = lags_ms(:)';
    tScore  = tScore(:)';

    if postOnly
        maskLag = (lags_ms >= 0);
    else
        maskLag = true(size(lags_ms));
    end

    if restrictAroundZero
        maskLag = maskLag & (lags_ms >= (centerMs - maxAbsLagMs)) & (lags_ms <= (centerMs + maxAbsLagMs));
    end

    good = (tScore >= tau) & maskLag;

    if ~any(good)
        w_ms = [nan nan];
        idxWin = [];
        return;
    end

    d = diff([0 good 0]);
    segStart = find(d==1);
    segEnd   = find(d==-1) - 1;

    bestScore = -inf;
    bestIdx   = [];

    for i=1:numel(segStart)
        idx = segStart(i):segEnd(i);

        dur = lags_ms(idx(end)) - lags_ms(idx(1));
        if dur < minContigMs
            continue;
        end

        if requireCrossZero
            if ~(lags_ms(idx(1)) <= 0 && lags_ms(idx(end)) >= 0)
                continue;
            end
        end

        score = sum(tScore(idx), 'omitnan');

        if score > bestScore
            bestScore = score;
            bestIdx = idx;
        end
    end

    if isempty(bestIdx)
        w_ms = [nan nan];
        idxWin = [];
        return;
    end

    idxWin = bestIdx;
    w_ms   = [lags_ms(idxWin(1)), lags_ms(idxWin(end))];
end

function plot_kernel_pooled_with_window_(ax, lags_ms, mu, sem, idxWin, col, ttl, titleFS, labelFS, tickFS, axLW, traceLW)
    hold(ax,'on');

    t = lags_ms(:)'/1000;

    patch(ax, [t fliplr(t)], [(mu-sem) fliplr(mu+sem)], col, 'FaceAlpha',0.15, 'EdgeColor','none');
    plot(ax, t, mu, '-', 'Color', col, 'LineWidth', traceLW);

    xline(ax, 0, 'k--', 'LineWidth', 2);
    yline(ax, 0, 'k-',  'LineWidth', 1.5);

    if ~isempty(idxWin)
        xs = [t(idxWin(1)) t(idxWin(end))];
        yl = ylim(ax);
        patch(ax, [xs(1) xs(2) xs(2) xs(1)], [yl(1) yl(1) yl(2) yl(2)], ...
            col, 'FaceAlpha',0.08, 'EdgeColor','none');
        ylim(ax, yl);
    end

    title(ax, ttl, 'FontWeight','bold','FontSize',titleFS);
    xlabel(ax, 'Lag (s)', 'FontSize', labelFS);
    ylabel(ax, 'Kernel (Hz/event)', 'FontSize', labelFS);
    set(ax,'FontSize',tickFS,'LineWidth',axLW,'TickDir','out','Box','off');
end

function plot_tscore_with_window_(ax, lags_ms, tScore, tau, idxWin, col, ttl, titleFS, labelFS, tickFS, axLW, traceLW)
    hold(ax,'on');

    t = lags_ms(:)'/1000;

    plot(ax, t, tScore, '-', 'Color', col, 'LineWidth', traceLW);
    yline(ax, tau, 'k--', 'LineWidth', 2);
    xline(ax, 0,   'k--', 'LineWidth', 2);

    if ~isempty(idxWin)
        xs = [t(idxWin(1)) t(idxWin(end))];
        yl = ylim(ax);
        patch(ax, [xs(1) xs(2) xs(2) xs(1)], [yl(1) yl(1) yl(2) yl(2)], ...
            col, 'FaceAlpha',0.08, 'EdgeColor','none');
        ylim(ax, yl);
    end

    title(ax, ttl, 'FontWeight','bold','FontSize',titleFS);
    xlabel(ax, 'Lag (s)', 'FontSize', labelFS);
    ylabel(ax, '|mean|/SEM', 'FontSize', labelFS);
    set(ax,'FontSize',tickFS,'LineWidth',axLW,'TickDir','out','Box','off');
end

function addSigIfNeeded_(ax, x1, x2, p, means, lw, fs)
    if ~isfinite(p) || p >= 0.05, return; end
    stars = pToStars_(p);
    yl = ylim(ax);
    yMaxData = max(means);
    yRange = yl(2) - yl(1);
    y = yMaxData + 0.08*yRange;
    h = 0.03*yRange;

    if y + h + 0.05*yRange > yl(2)
        ylim(ax, [yl(1), y + h + 0.08*yRange]);
        yl = ylim(ax); yRange = yl(2)-yl(1);
        y = yMaxData + 0.08*yRange;
        h = 0.03*yRange;
    end

    plot(ax, [x1 x1 x2 x2], [y y+h y+h y], 'k-', 'LineWidth', lw);
    text(ax, mean([x1 x2]), y+h + 0.01*yRange, stars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize', fs, 'FontWeight','bold', 'Color', [0 0 0]);
end

function s = pToStars_(p)
    if p < 0.001, s='***';
    elseif p < 0.01, s='**';
    else, s='*';
    end
end

% NOTE: permTestDiffMeans_ helper left intact for compatibility, but no longer
% used for Plot E (now uses sign-flip paired permutation to match your template script).
function [p, Tobs] = permTestDiffMeans_(a, b, nPerm)
    a = a(:); b = b(:);
    a = a(isfinite(a));
    b = b(isfinite(b));

    Tobs = mean(a,'omitnan') - mean(b,'omitnan');

    pool = [a; b];
    nA = numel(a);
    nP = numel(pool);

    rng(0);
    Tnull = nan(nPerm,1);
    for k = 1:nPerm
        perm = pool(randperm(nP));
        A = perm(1:nA);
        B = perm(nA+1:end);
        Tnull(k) = mean(A,'omitnan') - mean(B,'omitnan');
    end

    p = (1 + nnz(abs(Tnull) >= abs(Tobs))) / (nPerm + 1);
end