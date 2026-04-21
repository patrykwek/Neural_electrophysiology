function CELLTYPE_CODE()
%% Classify cells and plot (MSN/FSI): 3D feature space + population mean waveforms
% This is your classifier, with ONE ADDITION:
%   - Print, to the Command Window, counts of classified vs not-classified (with reasons)
%     for TWO session intervals:
%         (A) sessions 1–8
%         (B) sessions 9–50   (clamped to nSess)
%
% IMPORTANT INTERPRETATION:
%   - "Ignored / no spikes" = ct == [] (you set these to [] and skip them entirely)
%   - "Unclassified"        = ct == 0  (has spikes, but failed classification; reason_code set)
%   - "Classified"          = ct in {1,2}
%
% Printed categories per interval:
%   1) Total units encountered (channels x units) with waveform containers present
%   2) No spikes (ct==[]) counts
%   3) Classified counts (MSN/FSI + total)
%   4) Unclassified counts (ct==0) broken down by reason_code
%
% INPUT:
behMatFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';

clearvars -except behMatFile; clc;

%% ------------------ INPUTS ------------------
wfFieldCandidates = {'avgwaveforms','finalavgwaveforms','waveforms','spikewaveforms'};

% Outputs
baseDirWF_out = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION';
figDir = fullfile(baseDirWF_out,'1_FIGURE_CLASSIFICATION');
if ~exist(figDir,'dir'), mkdir(figDir); end

% ---- Canonical colors ----
col_msn  = [0.90 0.35 0.60];
col_fsi  = [0.35 0.55 0.75];
col_tan  = [0.95 0.55 0.20];
col_uncl = [0.65 0.65 0.65];
col_msn_3d = col_msn;
col_fsi_3d = col_fsi;
col_tan_3d = col_tan;

% feature extraction
fs_assume     = 30000;   % Hz
interp_len_ft = 1000;    % for feature timing

% waveform processing
contact_len         = 64;
interp_len_plot     = 200;
gauss_window_pts    = 15;

% ------------------ REFINED QC / REJECTION SETTINGS ------------------
min_peak_to_peak_uVlike = 1e-6;
min_trough_frac_of_ptp  = 0.10;

max_asymmetry = 2.0;
min_right_bins_for_ratio = 3;
min_left_bins_for_ratio  = 3;

min_edge_margin_bins_ft  = 0.10;
min_trough_depth_z = 0.5;

% ---- ADDED BACK: rounded-left rejection (left side too wide) ----
max_left_to_right_ratio = 1.2;

% ------------------ ADDED: pre-trough bump rejection gate ------------------
pretrough_win_bins = [10 200];      % bins before trough: 10–50
pretrough_peak_frac = 0.20;        % reject if prePeak > 0.35 * depth
% --------------------------------------------------------------------

% ------------------ NEW: "shaky waveform" + "fast downstroke" gates ------------------
max_rough_norm = 0.1;  % reject if rms(diff(wf_plot))/ptp_plot > this
min_left_hw_bins = 15;   % reject if left half-width (bins) is smaller than this
% ------------------------------------------------------------------------------------

% FR gates
FR_MSN_MIN = 0.1;   % (CHANGE 1) minimum MSN firing rate
FR_MSN_MAX = 10;
FR_FSI_MIN = 1.5;  FR_FSI_MAX = 1000;
FR_TAN_MIN = 0;    FR_TAN_MAX = 6;

% (CHANGE 2) ISI contamination exclusion rule:
% exclude units if >1% of ISIs are shorter than 2 ms
ISI_REFRACT_MS = 2;
MAX_REFRACT_FRAC = 0.01;

% diagnostics
plot_uncl_in_3d = false;             % do not plot unclassified units in 3D
max_uncl_waveforms_to_plot = 80;
save_diag_fig = true;

REAS = reason_map_();

%% ------------------ STYLE (single font sizes for ALL plots) ------------------
figPos = [120, 120, 900, 520];

FONT_SZ = 20;            % <-- single font size used everywhere
titleFontSize = FONT_SZ;
labelFontSize = FONT_SZ;
tickFontSize  = FONT_SZ;

axesTickLW    = 4.0;

%% ------------------ LOAD BEHAVIOR STRUCT ------------------
assert(exist(behMatFile,'file')==2, 'Missing MAT:\n%s', behMatFile);
S   = load(behMatFile);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in MAT.');
nSess = numel(beh);
fprintf('Behavior sessions: %d\n', nSess);

%% ------------------ OUTPUT CONTAINERS ------------------
cell_type       = cell(1, nSess); % []=no spikes (ignored), 0=Uncl, 1=MSN, 2=FSI
peak_to_valley  = cell(1, nSess); % us
hwhm            = cell(1, nSess); % us
amp             = cell(1, nSess); % trough amplitude
fr_hz           = cell(1, nSess); % Hz (rough)
isi_ms_cell     = cell(1, nSess); % ms
active_frac     = cell(1, nSess);

uncl_reason_code = cell(1, nSess);
uncl_reason_text = cell(1, nSess);

% For 3D plot (MSN/FSI)
PTV_ms   = [];
HWHM_ms  = [];
FR_hz_3d = [];
Type     = [];

% Unclassified overlay for 3D (kept, but not used when plot_uncl_in_3d=false)
PTV_ms_uncl  = [];
HWHM_ms_uncl = [];
FR_hz_uncl   = [];
UnclCode3D   = [];

% Mean waveform plots
Wf_MSN      = {};
Wf_FSI      = {};
Wf_TAN      = {};
Wf_UNCL     = {};
Wf_UNCL_code= [];

%% ------------------ MAIN LOOP ------------------
wfFieldUsedPrinted = false;

for sess = 1:nSess
    % Pull spikes
    spikes_by_ch = [];
    if isfield(beh(sess),'spikes') && iscell(beh(sess).spikes)
        spikes_by_ch = beh(sess).spikes;
    end

    % Pull waveforms
    wfs_sess = [];
    wfField = '';
    for k=1:numel(wfFieldCandidates)
        f = wfFieldCandidates{k};
        if isfield(beh(sess),f) && iscell(beh(sess).(f))
            wfs_sess = beh(sess).(f);
            wfField = f;
            break
        end
    end

    if isempty(wfs_sess)
        cell_type{sess}        = {};
        peak_to_valley{sess}   = {};
        hwhm{sess}             = {};
        amp{sess}              = {};
        fr_hz{sess}            = {};
        isi_ms_cell{sess}      = {};
        active_frac{sess}      = {};
        uncl_reason_code{sess} = {};
        uncl_reason_text{sess} = {};
        continue
    end

    if ~wfFieldUsedPrinted
        fprintf('Using waveform field "%s" from beh struct.\n', wfField);
        wfFieldUsedPrinted = true;
    end

    % allocate
    cell_type{sess}        = cell(size(wfs_sess));
    peak_to_valley{sess}   = cell(size(wfs_sess));
    hwhm{sess}             = cell(size(wfs_sess));
    amp{sess}              = cell(size(wfs_sess));
    fr_hz{sess}            = cell(size(wfs_sess));
    isi_ms_cell{sess}      = cell(size(wfs_sess));
    active_frac{sess}      = cell(size(wfs_sess));
    uncl_reason_code{sess} = cell(size(wfs_sess));
    uncl_reason_text{sess} = cell(size(wfs_sess));

    for ch = 1:numel(wfs_sess)
        units_ch = wfs_sess{ch};

        if isempty(units_ch)
            cell_type{sess}{ch}        = {};
            peak_to_valley{sess}{ch}   = {};
            hwhm{sess}{ch}             = {};
            amp{sess}{ch}              = {};
            fr_hz{sess}{ch}            = {};
            isi_ms_cell{sess}{ch}      = {};
            active_frac{sess}{ch}      = {};
            uncl_reason_code{sess}{ch} = {};
            uncl_reason_text{sess}{ch} = {};
            continue
        end

        nU = numel(units_ch);
        cell_type{sess}{ch}        = cell(1,nU);
        peak_to_valley{sess}{ch}   = cell(1,nU);
        hwhm{sess}{ch}             = cell(1,nU);
        amp{sess}{ch}              = cell(1,nU);
        fr_hz{sess}{ch}            = cell(1,nU);
        isi_ms_cell{sess}{ch}      = cell(1,nU);
        active_frac{sess}{ch}      = cell(1,nU);
        uncl_reason_code{sess}{ch} = cell(1,nU);
        uncl_reason_text{sess}{ch} = cell(1,nU);

        for u = 1:nU
            wf = units_ch{u};

            % ---- Require spike times ----
            hasSpikes = false;
            spk = [];
            if iscell(spikes_by_ch) && ch <= numel(spikes_by_ch) && ...
                    iscell(spikes_by_ch{ch}) && u <= numel(spikes_by_ch{ch})
                spk = spikes_by_ch{ch}{u};
                hasSpikes = ~isempty(spk);
            end
            if ~hasSpikes
                % ignore completely (no spikes)
                cell_type{sess}{ch}{u}        = [];
                peak_to_valley{sess}{ch}{u}   = [];
                hwhm{sess}{ch}{u}             = [];
                amp{sess}{ch}{u}              = [];
                fr_hz{sess}{ch}{u}            = [];
                isi_ms_cell{sess}{ch}{u}      = [];
                active_frac{sess}{ch}{u}      = [];
                uncl_reason_code{sess}{ch}{u} = [];
                uncl_reason_text{sess}{ch}{u} = [];
                continue
            end

            % ---- Empty waveform -> Unclassified ----
            if isempty(wf)
                cell_type{sess}{ch}{u}        = 0;
                active_frac{sess}{ch}{u}      = NaN;
                uncl_reason_code{sess}{ch}{u} = REAS.EMPTY_WAVEFORM;
                uncl_reason_text{sess}{ch}{u} = 'empty waveform';
                continue
            end

            % ---- ISI from spikes (timestamps assumed ms) ----
            isi_ms = NaN;
            spk = double(spk(:));
            if numel(spk) > 1
                isi_ms = mean(diff(spk));
            end
            isi_ms_cell{sess}{ch}{u} = isi_ms;

            % Exclude units with >1% ISIs < 2 ms
            if numel(spk) > 2
                isis = diff(spk); % ms
                frac_short = mean(isis < ISI_REFRACT_MS);
                if isfinite(frac_short) && frac_short > MAX_REFRACT_FRAC
                    cell_type{sess}{ch}{u}        = 0;
                    uncl_reason_code{sess}{ch}{u} = REAS.NO_CLASS_MATCH;
                    uncl_reason_text{sess}{ch}{u} = sprintf('ISI viol: %.2f%% < %d ms', 100*frac_short, ISI_REFRACT_MS);
                    continue
                end
            end

            % ---- FR from spikes (rough) ----
            fr_val = NaN;
            if numel(spk) > 1
                dur_ms = spk(end) - spk(1);
                if dur_ms > 0
                    fr_val = (numel(spk)-1) / (dur_ms/1000); % Hz
                end
            end
            fr_hz{sess}{ch}{u} = fr_val;
            active_frac{sess}{ch}{u} = NaN;

            % ---- Best contact ----
            v = wf(:).';
            selected = v;
            if numel(v) == 4*contact_len
                wf_mat = reshape(v, contact_len, 4)';          % 4 x 64
                [~, idxRow] = min(min(wf_mat, [], 2));
                selected = wf_mat(idxRow, :);
            elseif numel(v) == 2*contact_len
                wf_mat = reshape(v, contact_len, 2)';          % 2 x 64
                [~, idxRow] = min(min(wf_mat, [], 2));
                selected = wf_mat(idxRow, :);
            end

            aligned = selected; % no trough-centering

            % Must be trough-dominant
            if abs(min(aligned)) <= abs(max(aligned))
                cell_type{sess}{ch}{u}        = 0;
                uncl_reason_code{sess}{ch}{u} = REAS.NOT_TROUGH_DOMINANT;
                uncl_reason_text{sess}{ch}{u} = 'abs(min)<=abs(max)';
                continue
            end

            ptp = max(aligned) - min(aligned);
            if ~isfinite(ptp) || ptp <= min_peak_to_peak_uVlike
                cell_type{sess}{ch}{u}        = 0;
                uncl_reason_code{sess}{ch}{u} = REAS.NOT_TROUGH_DOMINANT;
                uncl_reason_text{sess}{ch}{u} = 'flat/invalid waveform after polarity normalize';
                maybe_store_uncl_(interp1(1:numel(aligned), aligned, linspace(1, numel(aligned), interp_len_plot), 'linear', 'extrap'), ...
                    uncl_reason_code{sess}{ch}{u});
                continue
            end
            tr = min(aligned);
            if abs(tr) < min_trough_frac_of_ptp * ptp
                cell_type{sess}{ch}{u}        = 0;
                uncl_reason_code{sess}{ch}{u} = REAS.NOT_TROUGH_DOMINANT;
                uncl_reason_text{sess}{ch}{u} = sprintf('weak trough: |trough| < %.2f*ptp', min_trough_frac_of_ptp);
                maybe_store_uncl_(interp1(1:numel(aligned), aligned, linspace(1, numel(aligned), interp_len_plot), 'linear', 'extrap'), ...
                    uncl_reason_code{sess}{ch}{u});
                continue
            end

            % ---- Interpolate ----
            t_src    = 1:numel(aligned);
            wf_plot  = interp1(t_src, aligned, linspace(1, numel(aligned), interp_len_plot), 'linear', 'extrap');
            x_dst_ft = linspace(1, numel(aligned), interp_len_ft);
            big_i    = interp1(t_src, aligned, x_dst_ft, 'linear', 'extrap');

            [~, trough_bin] = min(big_i);

            micros_per_bin = (1e6 / fs_assume) * (numel(aligned) / interp_len_ft);

            [~, peak_rel] = max(big_i(trough_bin:end));
            peak_bin      = trough_bin + peak_rel - 1;
            ptv_us        = (peak_bin - trough_bin) * micros_per_bin;

            base_n   = max(5, round(0.1 * interp_len_ft));
            baseline = mean(big_i(1:base_n));
            min_tr   = min(big_i);
            depth    = baseline - min_tr;

            % pre-trough bump rejection gate
            pre_hi = max(1, trough_bin - pretrough_win_bins(1)); % trough-10
            pre_lo = max(1, trough_bin - pretrough_win_bins(2)); % trough-50
            if pre_hi >= pre_lo
                prePeak = max(big_i(pre_lo:pre_hi));
                if isfinite(prePeak) && isfinite(depth) && (prePeak > pretrough_peak_frac * depth)
                    cell_type{sess}{ch}{u}        = 0;
                    uncl_reason_code{sess}{ch}{u} = REAS.NO_CLASS_MATCH;
                    uncl_reason_text{sess}{ch}{u} = sprintf('pretrough bump: prePeak=%.4g > %.2f*depth', prePeak, pretrough_peak_frac);
                    maybe_store_uncl_(wf_plot, uncl_reason_code{sess}{ch}{u});
                    continue
                end
            end

            half_min = min_tr + depth/2;

            [~, left_bin]  = min(abs(big_i(1:trough_bin) - half_min));
            [~, right_rel] = min(abs(big_i(trough_bin:end) - half_min));
            right_bin      = trough_bin + right_rel - 1;
            hwhm_us        = (right_bin - left_bin) * micros_per_bin;

            % rounded-left gate
            left_hw_bins  = trough_bin - left_bin;
            right_hw_bins = right_bin  - trough_bin;
            lr_ratio      = left_hw_bins / max(right_hw_bins, 1);

            if lr_ratio > max_left_to_right_ratio
                cell_type{sess}{ch}{u}        = 0;
                uncl_reason_code{sess}{ch}{u} = REAS.ROUNDED_LEFT;
                uncl_reason_text{sess}{ch}{u} = sprintf('lr_ratio=%.3f>%.3f', lr_ratio, max_left_to_right_ratio);
                maybe_store_uncl_(wf_plot, uncl_reason_code{sess}{ch}{u});
                continue
            end

            peak_to_valley{sess}{ch}{u} = ptv_us;
            hwhm{sess}{ch}{u}           = hwhm_us;
            amp{sess}{ch}{u}            = min_tr;

            % ---- FR gates ----
            FRokFSI = (isnan(fr_val) || (fr_val >= FR_FSI_MIN && fr_val <= FR_FSI_MAX));
            FRokMSN = (isnan(fr_val) || (fr_val >= FR_MSN_MIN && fr_val <= FR_MSN_MAX));

            % ---- Type rules ----
            ct = 0;
            if isfinite(hwhm_us) && isfinite(ptv_us) && isfinite(isi_ms)
                if (hwhm_us < 250) && (ptv_us < 350) && FRokFSI %FSI
                    ct = 2;
                elseif (hwhm_us > 75) && (hwhm_us < 300) && (ptv_us > 350) && FRokMSN %MSN
                    ct = 1;
                else
                    ct = 0;
                    uncl_reason_code{sess}{ch}{u} = REAS.NO_CLASS_MATCH;
                    uncl_reason_text{sess}{ch}{u} = 'features did not satisfy MSN/FSI rules';
                end
            else
                ct = 0;
                uncl_reason_code{sess}{ch}{u} = REAS.MISSING_FEATURES;
                uncl_reason_text{sess}{ch}{u} = 'non-finite (hwhm/ptv/isi)';
            end
            cell_type{sess}{ch}{u} = ct;

            if ct==0
                maybe_store_uncl_(wf_plot, uncl_reason_code{sess}{ch}{u});
            end

            % ---- Collect 3D features for MSN/FSI ----
            if (ct==1 || ct==2) && isfinite(ptv_us) && isfinite(hwhm_us) && isfinite(fr_val)
                PTV_ms   = [PTV_ms;  ptv_us/1000];  %#ok<AGROW>
                HWHM_ms  = [HWHM_ms; hwhm_us/1000]; %#ok<AGROW>
                FR_hz_3d = [FR_hz_3d; fr_val];      %#ok<AGROW>
                Type     = [Type;    ct];           %#ok<AGROW>
            end

            % ---- Collect normalized waveforms for population mean ----
            if ct==1 || ct==2
                wf_norm = zscore(wf_plot, 0, 2);
                if ct==1
                    Wf_MSN{end+1} = wf_norm; %#ok<AGROW>
                elseif ct==2
                    Wf_FSI{end+1} = wf_norm; %#ok<AGROW>
                else
                    Wf_TAN{end+1} = wf_norm; %#ok<AGROW>
                end
            end
        end
    end
end

fprintf('Waveforms collected -> MSN: %d, FSI: %d, Uncl-wf-saved: %d\n', numel(Wf_MSN), numel(Wf_FSI), numel(Wf_UNCL));

%% ------------------ PRINT UNCLASSIFIED DIAGNOSTICS (ALL) ------------------
print_uncl_diagnostics_(cell_type, uncl_reason_code, REAS);

%% ================== NEW: PRINT CELL TYPES PER SESSION + MEAN/SEM ==================
classified_per_sess = zeros(1, nSess);
msn_per_sess = zeros(1, nSess);
fsi_per_sess = zeros(1, nSess);

for s = 1:nSess
    nCls = 0;
    nMSN_s = 0;
    nFSI_s = 0;

    if isempty(cell_type{s})
        classified_per_sess(s) = 0;
        msn_per_sess(s) = 0;
        fsi_per_sess(s) = 0;
        continue;
    end

    for ch = 1:numel(cell_type{s})
        if isempty(cell_type{s}{ch}), continue; end
        for u = 1:numel(cell_type{s}{ch})
            ct = cell_type{s}{ch}{u};
            if isempty(ct), continue; end
            if ~isnumeric(ct) || ~isscalar(ct), continue; end

            if ct==1
                nMSN_s = nMSN_s + 1;
                nCls = nCls + 1;
            elseif ct==2
                nFSI_s = nFSI_s + 1;
                nCls = nCls + 1;
            end
        end
    end

    classified_per_sess(s) = nCls;
    msn_per_sess(s) = nMSN_s;
    fsi_per_sess(s) = nFSI_s;
end

fprintf('\n================= CELL TYPES PER SESSION =================\n');
for s = 1:nSess
    fprintf('Session %d: MSN = %d, FSI = %d, Total classified = %d\n', ...
        s, msn_per_sess(s), fsi_per_sess(s), classified_per_sess(s));
end
m_cls = mean(classified_per_sess);
sem_cls = std(classified_per_sess) / sqrt(numel(classified_per_sess));
fprintf('Mean classified per session: %.4f\n', m_cls);
fprintf('SEM classified per session : %.4f\n', sem_cls);
fprintf('==========================================================\n\n');

%% ================== NEW: INTERVAL COUNTS (REQUESTED) ==================
intervals = [1 8; 9 38];

fprintf('\n================= CLASSIFICATION SUMMARY BY INTERVAL =================\n');
for ii = 1:size(intervals,1)
    a = intervals(ii,1);
    b = intervals(ii,2);
    a = max(1, min(nSess, a));
    b = max(1, min(nSess, b));
    if b < a
        fprintf('Interval %d invalid after clamp. Skipping.\n', ii);
        continue;
    end
    sess_idx = a:b;

    stats = interval_stats_(cell_type, uncl_reason_code, sess_idx, REAS);

    fprintf('\n-- Sessions %d-%d --\n', a, b);
    fprintf('Total unit slots (ch x unit, waveform-present) : %d\n', stats.totalSlots);
    fprintf('No spikes (ct==[])                            : %d\n', stats.noSpikes);
    fprintf('Classified (ct in {1,2})                      : %d\n', stats.classifiedTotal);
    fprintf('  MSN (1): %d | FSI (2): %d\n', stats.nMSN, stats.nFSI);
    fprintf('Unclassified (ct==0; has spikes)              : %d\n', stats.unclassifiedTotal);

    if stats.unclassifiedTotal > 0
        print_reason_table_(stats.unclCodes, REAS, sprintf('Unclassified reasons (Sessions %d-%d)', a, b));
    else
        fprintf('No unclassified units with spikes in this interval.\n');
    end
end
fprintf('======================================================================\n\n');

%% ------------------ FIG 1: 3D feature space (MSN/FSI only) ------------------
f1 = figure('Color','w','Position',figPos);
ax1 = axes('Parent',f1); hold(ax1,'on');

scatter3(ax1, HWHM_ms(Type==1), PTV_ms(Type==1), FR_hz_3d(Type==1), 30, col_msn_3d, 'filled', 'DisplayName','MSN');
scatter3(ax1, HWHM_ms(Type==2), PTV_ms(Type==2), FR_hz_3d(Type==2), 30, col_fsi_3d, 'filled', 'DisplayName','FSI');

hx = xlabel(ax1,'Valley width (ms)','FontSize',labelFontSize);
hy = ylabel(ax1,'Peak to valley (ms)','FontSize',labelFontSize);
hz = zlabel(ax1,'Firing rate (Hz)','FontSize',labelFontSize);
title(ax1,'Feature Space (MSN/FSI)','FontWeight','bold','FontSize',titleFontSize);

set(ax1,'ZScale','log','FontSize',tickFontSize,'Box','off','TickDir','out', ...
    'LineWidth',axesTickLW,'Layer','top');
grid(ax1,'on'); view(ax1,135,25); axis(ax1,'square');

% Make Peak-to-valley (Y) grid/ticks at every 0.1
yl = ylim(ax1);
yMaxTick = ceil(yl(2)*5)/5;
yMinTick = floor(yl(1)*5)/5;
ax1.YTick = yMinTick:0.2:yMaxTick;

% Keep label rotation behavior (do not change)
set([hx hy hz], 'RotationMode','auto');
xtickangle(ax1, 0);
ytickangle(ax1, 0);

% (CHANGE) Remove legend entirely
% lg = legend(ax1,'Location','northeastoutside');
% set(lg,'FontSize',tickFontSize);

export_or_save_(f1, fullfile(figDir,'feature_space_3d_MSN_FSI.png'), 300);

%% ------------------ FIG 2: population mean waveforms (REVERTED) ------------------
centroid_fun = @(M) smoothdata(mean(M,1,'omitnan'), 'gaussian', gauss_window_pts);
MsnMat = cell2mat_or_empty_(Wf_MSN);
FsiMat = cell2mat_or_empty_(Wf_FSI);
TanMat = cell2mat_or_empty_(Wf_TAN); %#ok<NASGU>

orig_samples = contact_len;
total_ms     = (orig_samples / fs_assume) * 1000;
t_ms         = linspace(-total_ms/2, total_ms/2, interp_len_plot);

f2 = figure('Color','w','Position',[120 120 1000 380]);

subplot(1,2,1); hold on; axis off
if ~isempty(MsnMat)
    msn = centroid_fun(MsnMat);
    [~, cidx] = min(msn);
    msn = circshift(msn, round(numel(msn)/2) - cidx);
    plot([t_ms(1) t_ms(end)], [0 0], 'k:', 'LineWidth', 1);
    plot(t_ms, msn, 'Color', col_msn, 'LineWidth', 2.5);
    % (CHANGE) no "n=" printed
    title('MSN', 'FontSize', titleFontSize, 'FontWeight','bold');
    draw_ms_scalebar_(t_ms, msn, 1.0, tickFontSize);
end

subplot(1,2,2); hold on; axis off
if ~isempty(FsiMat)
    fsi = centroid_fun(FsiMat);
    [~, cidx] = min(fsi);
    fsi = circshift(fsi, round(numel(fsi)/2) - cidx);
    plot([t_ms(1) t_ms(end)], [0 0], 'k:', 'LineWidth', 1);
    plot(t_ms, fsi, 'Color', col_fsi, 'LineWidth', 2.5);
    % (CHANGE) no "n=" printed
    title('FSI', 'FontSize', titleFontSize, 'FontWeight','bold');
    draw_ms_scalebar_(t_ms, fsi, 1.0, tickFontSize);
end

sgtitle('Population Mean Waveforms — MSN & FSI','FontWeight','bold','FontSize',titleFontSize);
export_or_save_(f2, fullfile(figDir,'population_mean_waveforms_MSN_FSI.png'), 300);

%% ------------------ FIG 2b: example unclassified waveforms (REVERTED) ------------------
if ~isempty(Wf_UNCL)
    UnclMat = cell2mat_or_empty_(Wf_UNCL);
    f2b = figure('Color','w','Position',[120 120 980 420]); hold on;
    plot(UnclMat', 'Color', [col_uncl 0.12]);
    yline(0,'k:');
    title(sprintf('Unclassified example waveforms (n=%d shown)', size(UnclMat,1)), ...
        'FontWeight','bold','FontSize',titleFontSize);
    xlabel('Interpolated timepoints','FontSize',labelFontSize);
    ylabel('Z-scored amplitude','FontSize',labelFontSize);
    set(gca,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,'Box','off','Layer','top');
    grid on; axis tight;
    if save_diag_fig
        export_or_save_(f2b, fullfile(figDir,'unclassified_waveforms_examples.png'), 300);
    end
end

%% ------------------ FIG 3: examples + centroid (MSN/FSI) ------------------
f3 = figure('Color','w','Position',figPos);

subplot(1,2,1); hold on;
if ~isempty(MsnMat)
    msn_centroid = centroid_fun(MsnMat);
    [~, cidx]    = min(msn_centroid);
    msn_centroid = circshift(msn_centroid, round(numel(msn_centroid)/2) - cidx);

    cvals = nan(size(MsnMat,1),1);
    for ii = 1:size(MsnMat,1)
        a = MsnMat(ii,:);
        b = msn_centroid(:).';
        if all(isfinite(a)) && all(isfinite(b))
            r = corrcoef(a, b);
            cvals(ii) = r(1,2);
        end
    end

    corr_min = 0.85;
    keep = isfinite(cvals) & (cvals >= corr_min);

    plot(MsnMat(keep,:)', 'Color', [col_msn 0.12]);
    plot(msn_centroid, 'Color', col_msn, 'LineWidth', 5);
end
title(sprintf('MSNs (n=%d)', size(MsnMat,1)), 'FontWeight','bold','FontSize',titleFontSize);
yline(0,'k:','LineWidth',2.5); grid on; box on;
xlabel('Interpolated timepoints (center-aligned)','FontSize',labelFontSize);
ylabel('Z-scored amplitude','FontSize',labelFontSize);
set(gca,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,'Box','off','Layer','top');
axis tight;

subplot(1,2,2); hold on;
if ~isempty(FsiMat)
    plot(FsiMat', 'Color', [col_fsi 0.12]);
    fsi_centroid = centroid_fun(FsiMat);
    [~, cidx]    = min(fsi_centroid);
    fsi_centroid = circshift(fsi_centroid, round(numel(fsi_centroid)/2) - cidx);
    plot(fsi_centroid, 'Color', col_fsi, 'LineWidth', 5);
end
title(sprintf('FSIs (n=%d)', size(FsiMat,1)), 'FontWeight','bold','FontSize',titleFontSize);
yline(0,'k:','LineWidth',2.5); grid on; box on;
xlabel('Interpolated timepoints (center-aligned)','FontSize',labelFontSize);
set(gca,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,'Box','off','Layer','top');
axis tight;

sgtitle('Waveforms (Examples + Smoothed Centroid) — MSN & FSI','FontWeight','bold','FontSize',titleFontSize);
export_or_save_(f3, fullfile(figDir,'waveforms_examples_plus_centroid_MSN_FSI.png'), 300);

%% ------------------ FIG 4: counts (MSN/FSI) ------------------
[~, countMSN, countFSI, countTAN] = count_types_(cell_type); %#ok<NASGU>
vals   = [countMSN, countFSI];
labels = {'MSN','FSI'};

f4 = figure('Color','w','Position',figPos);
ax4 = axes(f4); hold(ax4,'on');

b = bar(ax4, vals, 'FaceColor','flat', 'BarWidth', 0.6);
b.CData(1,:) = col_msn;
b.CData(2,:) = col_fsi;

yMax = max(vals) + 50;
ylim(ax4, [0 yMax]);

set(ax4,'XTick',1:2,'XTickLabel',labels,'FontSize',tickFontSize, ...
    'TickDir','out','LineWidth',axesTickLW,'Box','off','Layer','top');
ylabel(ax4,'Number of cells','FontSize',labelFontSize);
title(ax4,'Cell-type counts (MSN / FSI)','FontWeight','bold','FontSize',titleFontSize);

grid(ax4,'on');
ax4.XGrid = 'off';
ax4.YGrid = 'on';

offset = max([1, vals])*0.03;
text(ax4, 1:2, vals+offset, string(vals), ...
    'HorizontalAlignment','center','FontWeight','bold','FontSize',tickFontSize);

export_or_save_(f4, fullfile(figDir,'celltype_histogram_counts_MSN_FSI.png'), 300);

%% ------------------ FIG 5: stacked histograms (same fields, but ISI/FR are from spikes) ------------------
ptv_msn=[]; ptv_fsi=[]; ptv_tan=[];
hwhm_msn=[]; hwhm_fsi=[]; hwhm_tan=[];
fr_msn=[]; fr_fsi=[]; fr_tan=[];
isi_msn=[]; isi_fsi=[]; isi_tan=[];

for s = 1:numel(cell_type)
    if isempty(cell_type{s}), continue; end
    for ch = 1:numel(cell_type{s})
        if isempty(cell_type{s}{ch}), continue; end
        for u = 1:numel(cell_type{s}{ch})
            vct = cell_type{s}{ch}{u};
            if isempty(vct) || ~isscalar(vct), continue; end
            if ~ismember(vct,[1 2]), continue; end

            ptv = try_get_num_(peak_to_valley,s,ch,u)/1000;
            hwh = try_get_num_(hwhm,         s,ch,u)/1000;
            frv = try_get_num_(fr_hz,        s,ch,u);
            isi = try_get_num_(isi_ms_cell,  s,ch,u);

            if isfinite(ptv)
                switch vct
                    case 1, ptv_msn(end+1)=ptv; %#ok<AGROW>
                    case 2, ptv_fsi(end+1)=ptv; %#ok<AGROW>
                    case 3, ptv_tan(end+1)=ptv; %#ok<AGROW>
                end
            end
            if isfinite(hwh)
                switch vct
                    case 1, hwhm_msn(end+1)=hwh; %#ok<AGROW>
                    case 2, hwhm_fsi(end+1)=hwh; %#ok<AGROW>
                    case 3, hwhm_tan(end+1)=hwh; %#ok<AGROW>
                end
            end
            if isfinite(frv)
                switch vct
                    case 1, fr_msn(end+1)=frv; %#ok<AGROW>
                    case 2, fr_fsi(end+1)=frv; %#ok<AGROW>
                    case 3, fr_tan(end+1)=frv; %#ok<AGROW>
                end
            end
            if isfinite(isi)
                switch vct
                    case 1, isi_msn(end+1)=isi; %#ok<AGROW>
                    case 2, isi_fsi(end+1)=isi; %#ok<AGROW>
                    case 3, isi_tan(end+1)=isi; %#ok<AGROW>
                end
            end
        end
    end
end

make_stacked_histograms_(...
    {ptv_msn, hwhm_msn, fr_msn,  isi_msn}, ...
    {ptv_fsi, hwhm_fsi, fr_fsi,  isi_fsi}, ...
    {'Peak-to-Valley Time (ms)','HWHM (ms)','Avg Firing Rate (Hz)','Avg ISI (ms)'}, ...
    {'Peak-to-Valley','HWHM','Avg Firing Rate','Avg ISI'}, ...
    [col_msn; col_fsi], figDir, titleFontSize, labelFontSize, tickFontSize, axesTickLW, figPos);

%% ------------------ SAVE outputs ------------------
compatMat = fullfile(baseDirWF_out, 'celltype_all_sessions.mat');
cell_type = sanitize_celltype_(cell_type);
save(compatMat, 'cell_type', '-v7.3');

fullMat = fullfile(baseDirWF_out, 'celltype_all_sessions_full.mat');
save(fullMat, 'cell_type','peak_to_valley','hwhm','amp','fr_hz','isi_ms_cell','active_frac', ...
    'uncl_reason_code','uncl_reason_text', ...
    'fs_assume','interp_len_ft', ...
    'FR_MSN_MAX','FR_FSI_MIN','FR_FSI_MAX','FR_TAN_MIN','FR_TAN_MAX', ...
    'behMatFile', ...
    'max_asymmetry','min_edge_margin_bins_ft','min_trough_depth_z', ...
    'min_trough_frac_of_ptp', ...
    'max_left_to_right_ratio', ...
    '-v7.3');

fprintf('Saved raster-compatible cell types to %s\n', compatMat);
fprintf('Saved full outputs to %s\n', fullMat);

%% nested helper
    function maybe_store_uncl_(wf_plot_local, rc)
        if numel(Wf_UNCL) >= max_uncl_waveforms_to_plot, return; end
        wf_norm_u = zscore(wf_plot_local, 0, 2);
        Wf_UNCL{end+1} = wf_norm_u; %#ok<AGROW>
        if isempty(rc) || ~isscalar(rc), rc = REAS.UNKNOWN; end
        Wf_UNCL_code(end+1) = rc; %#ok<AGROW>
    end
end

%% ======================== LOCAL FUNCTIONS ========================

function stats = interval_stats_(cell_type, uncl_reason_code, sess_idx, REAS)
    stats = struct();
    stats.totalSlots = 0;
    stats.noSpikes = 0;
    stats.classifiedTotal = 0;
    stats.unclassifiedTotal = 0;
    stats.nMSN = 0;
    stats.nFSI = 0;
    stats.nTAN = 0;
    stats.unclCodes = [];

    for s = sess_idx
        if s < 1 || s > numel(cell_type), continue; end
        if isempty(cell_type{s}), continue; end

        for ch = 1:numel(cell_type{s})
            if isempty(cell_type{s}{ch}), continue; end

            for u = 1:numel(cell_type{s}{ch})
                stats.totalSlots = stats.totalSlots + 1;

                ct = cell_type{s}{ch}{u};

                if isempty(ct)
                    stats.noSpikes = stats.noSpikes + 1;
                    continue
                end
                if ~isnumeric(ct) || ~isscalar(ct)
                    ct = 0;
                end

                if ismember(ct, [1 2])
                    stats.classifiedTotal = stats.classifiedTotal + 1;
                    if ct==1, stats.nMSN = stats.nMSN + 1; end
                    if ct==2, stats.nFSI = stats.nFSI + 1; end
                elseif ct==0
                    stats.unclassifiedTotal = stats.unclassifiedTotal + 1;
                    rc = [];
                    try
                        rc = uncl_reason_code{s}{ch}{u};
                    catch
                        rc = [];
                    end
                    if isempty(rc) || ~isnumeric(rc) || ~isscalar(rc)
                        rc = REAS.UNKNOWN;
                    end
                    stats.unclCodes(end+1) = rc; %#ok<AGROW>
                else
                    stats.unclassifiedTotal = stats.unclassifiedTotal + 1;
                    stats.unclCodes(end+1) = REAS.UNKNOWN; %#ok<AGROW>
                end
            end
        end
    end
end

function REAS = reason_map_()
    REAS = struct();
    REAS.UNKNOWN             = 0;
    REAS.EMPTY_WAVEFORM      = 1;
    REAS.NOT_TROUGH_DOMINANT = 2;
    REAS.ROUNDED_LEFT        = 3;
    REAS.MISSING_FEATURES    = 4;
    REAS.NO_CLASS_MATCH      = 5;
end

function print_uncl_diagnostics_(cell_type, uncl_reason_code, REAS)
    fprintf('\n================= UNCLASSIFIED DIAGNOSTICS =================\n');
    codes_all = collect_uncl_codes_(cell_type, uncl_reason_code, [], REAS);
    print_reason_table_(codes_all, REAS, 'All sessions (combined)');
    fprintf('============================================================\n\n');
end

function codes = collect_uncl_codes_(cell_type, uncl_reason_code, sess_idx, REAS)
    codes = [];
    if isempty(sess_idx), sess_idx = 1:numel(cell_type); end
    for s = sess_idx
        if isempty(cell_type{s}), continue; end
        for ch = 1:numel(cell_type{s})
            if isempty(cell_type{s}{ch}), continue; end
            for u = 1:numel(cell_type{s}{ch})
                ct = cell_type{s}{ch}{u};
                if isempty(ct), continue; end
                if ~isnumeric(ct) || ~isscalar(ct), ct = 0; end
                if ct ~= 0, continue; end
                rc = uncl_reason_code{s}{ch}{u};
                if isempty(rc) || ~isnumeric(rc) || ~isscalar(rc), rc = REAS.UNKNOWN; end
                codes(end+1) = rc; %#ok<AGROW>
            end
        end
    end
end

function print_reason_table_(codes, REAS, header)
    fprintf('\n-- %s --\n', header);
    total = numel(codes);
    if total==0
        fprintf('No unclassified units with spike-times in this subset.\n');
        return
    end

    ordered = {'EMPTY_WAVEFORM','NOT_TROUGH_DOMINANT','ROUNDED_LEFT','MISSING_FEATURES','NO_CLASS_MATCH','UNKNOWN'};
    fprintf('%-22s  %8s  %8s\n', 'Reason', 'Count', 'Percent');
    fprintf('%s\n', repmat('-',1,44));

    for k = 1:numel(ordered)
        nm = ordered{k};
        code = REAS.(nm);
        c = sum(codes==code);
        fprintf('%-22s  %8d  %7.2f%%\n', nm, c, 100*c/max(1,total));
    end
    fprintf('%s\n', repmat('-',1,44));
    fprintf('%-22s  %8d\n', 'TOTAL Unclassified', total);
end

function export_or_save_(fig, outPath, dpi)
    % Save PNG (existing behavior)
    try
        exportgraphics(fig, outPath, 'Resolution', dpi);
    catch
        saveas(fig, outPath);
    end

    % (CHANGE) Also save SVG alongside PNG
    [p,n,~] = fileparts(outPath);
    svgPath = fullfile(p, [n '.svg']);
    try
        exportgraphics(fig, svgPath, 'ContentType','vector');
    catch
        try
            print(fig, svgPath, '-dsvg');
        catch
            % do nothing
        end
    end
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit','ratBEHstruct2','ratBEHstruct'};
    for k=1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k})), beh = S.(cands{k}); return; end
    end
    f = fieldnames(S);
    for i = 1:numel(f)
        v = S.(f{i}); if isstruct(v) && numel(v)>1, beh = v; return; end
    end
    for i = 1:numel(f)
        v = S.(f{i}); if isstruct(v), beh = v; return; end
    end
end

function num = try_get_num_(C,s,ch,u)
    num = NaN;
    try
        v = C{s}{ch}{u};
        if isnumeric(v) && isscalar(v), num = double(v); end
    catch
        num = NaN;
    end
end

function M = cell2mat_or_empty_(C)
    if isempty(C), M = []; return; end
    M = cell2mat(reshape(C,1,1,[]));
    M = squeeze(M)';  % N x T
end

function draw_ms_scalebar_(t_ms, yvec, one_ms_len, fontSize)
    xl = [min(t_ms), max(t_ms)];
    yl = [min(yvec), max(yvec)];
    pad_ms = 1.3;
    x1 = xl(2) - pad_ms - one_ms_len;
    x2 = x1 + one_ms_len;
    y  = yl(1) + 0.10*(yl(2)-yl(1));
    plot([x1 x2], [y y], 'k', 'LineWidth', 2);
    text((x1+x2)/2, y - 0.08*(yl(2)-yl(1)), sprintf('%g ms', one_ms_len), ...
        'HorizontalAlignment','center','VerticalAlignment','top', ...
        'FontSize',fontSize,'Color','k');
    xlim(xl); ylim([yl(1)-0.05*(yl(2)-yl(1)), yl(2)+0.05*(yl(2)-yl(1))]);
end

function ct_clean = sanitize_celltype_(ct)
    ct_clean = ct;
    for s = 1:numel(ct_clean)
        if isempty(ct_clean{s}), continue; end
        for ch = 1:numel(ct_clean{s})
            if isempty(ct_clean{s}{ch}), continue; end
            for u = 1:numel(ct_clean{s}{ch})
                v = ct_clean{s}{ch}{u};
                if isempty(v), continue; end
                if ~isnumeric(v) || ~isscalar(v), ct_clean{s}{ch}{u} = 0; end
            end
        end
    end
end

function [nUncl, nMSN, nFSI, nTAN] = count_types_(ct)
    nUncl=0; nMSN=0; nFSI=0; nTAN=0;
    for s = 1:numel(ct)
        if isempty(ct{s}), continue; end
        for ch = 1:numel(ct{s})
            if isempty(ct{s}{ch}), continue; end
            for u = 1:numel(ct{s}{ch})
                v = ct{s}{ch}{u};
                if isempty(v), continue; end
                if ~isnumeric(v) || ~isscalar(v), v = 0; end
                switch v
                    case 1, nMSN = nMSN + 1;
                    case 2, nFSI = nFSI + 1;
                    case 3, nTAN = nTAN + 1;
                    otherwise, nUncl = nUncl + 1;
                end
            end
        end
    end
end

function make_stacked_histograms_(D_msn, D_fsi, xlabels_cell, titles_cell, colors2x, figDir, titleFS, labelFS, tickFS, axLW, figPos)
    f = figure('Color','w','Position',figPos);
    tl = tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');

    for i = 1:4
        ax = nexttile(tl,i); hold(ax,'on');

        msn = D_msn{i}; fsi = D_fsi{i};
        ALL = [msn(:); fsi(:)];
        ALL = ALL(isfinite(ALL));
        if isempty(ALL)
            title(ax,[titles_cell{i} ' (No data)'],'FontSize',titleFS,'FontWeight','bold');
            box(ax,'on'); continue;
        end

        nbins = 25;
        if min(ALL)==max(ALL), nbins = 5; end
        edges = linspace(min(ALL), max(ALL), nbins+1);

        C = zeros(2,numel(edges)-1);
        C(1,:) = histcounts(msn, edges);
        C(2,:) = histcounts(fsi, edges);

        ctrs = edges(1:end-1) + diff(edges)/2;
        b = bar(ax, ctrs, C', 'stacked', 'BarWidth', 1.0);
        b(1).FaceColor = colors2x(1,:); b(1).EdgeColor = 'none';
        b(2).FaceColor = colors2x(2,:); b(2).EdgeColor = 'none';

        xlabel(ax, xlabels_cell{i}, 'FontSize',labelFS);
        ylabel(ax, 'Count', 'FontSize',labelFS);
        title(ax, titles_cell{i}, 'FontWeight','bold','FontSize',titleFS);
        set(ax,'FontSize',tickFS,'TickDir','out','LineWidth',axLW,'Box','off','Layer','top');
        grid(ax,'on');

        if i==1
            lg = legend(ax, {'MSN','FSI'}, 'Location','northeast');
            set(lg,'FontSize',tickFS);
        end
    end

    sgtitle(tl,'Stacked Histograms by Cell Type','FontWeight','bold','FontSize',titleFS);

    % (CHANGE) also save SVG for stacked histogram figure
    export_or_save_(f, fullfile(figDir,'stacked_histograms_types.png'), 300);
end