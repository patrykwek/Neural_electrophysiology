%% ========= FIG3A_SESSION_ENSEMBLE_HEATMAP__A_THESIS_FINAL__GLOBALWARP__FIG2TRIALS_EXACT__CLASSIFIEDONLY.m =========
% ONE heatmap for sessions 1–38 (all on the same plot).
%
% Key properties (kept from your script):
%   1) TRIAL SELECTION matches FIG2-like behavioral rules (valid, RT>=MIN_RT, press/lick in fixedWin)
%   2) PER-UNIT trial inclusion additionally requires >=1 spike in [cue+fixedWin(1), cue+fixedWin(2)]
%   3) ONLY CLASSIFIED units are used: cell_type in {1,2,3} (MSN/FSI/TAN)
%   4) GLOBAL warp targets computed ONCE (behavior-only, no spike requirement)
%   5) For each unit: compute warped, per-trial z-scored SDF; use half-trial plot (even) mean

clear; clc;

%% ---- USER SETTINGS ----
matFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
cellTypeFile = '/Volumes/WD_BLACK/A_THESIS_FINAL/1_FIGURE_CLASSIFICATION/celltype_all_sessions_full.mat';

baseOut = '/Volumes/WD_BLACK/A_THESIS_FINAL';
outDir  = fullfile(baseOut, 'MATRIX');
if ~exist(outDir,'dir'), mkdir(outDir); end

preCueMs   = 500;
postLickMs = 3000;

fixedWin = [-1000, 6000];

MIN_RT_MS = 100;
minTrialsPerUnit = 10;
requireValid = true;

dt_ms        = 10;
gaussSigmaMs = 25;

climFixed = [-2 2];

xUnit = "s";
forceIntegerSecondsTicks = true;

% ---- FIGURE SIZE (make overall figure ~half as wide) ----
figPos = [120, 120, 650, 900];

cbFontSize = 22;

titleFontSize = 30;
labelFontSize = 30;
tickFontSize  = 30;

eventLineLW   = 5.0;
axesTickLW    = 4.0;
axesTickLen   = [0.02 0.02];

colCue   = [1 0 0];
colPress = [0.10 0.55 0.95];
colLick  = [0.15 0.70 0.20];

%% ---- LOAD BEHSTRUCT ----
assert(exist(matFile,'file')==2);
S = load(matFile);
beh = pickBehStruct_(S);
nSessions = numel(beh);
nSessions = min(nSessions, 36);
beh = beh(1:nSessions);

%% ---- LOAD CELL TYPES ----
CT = load(cellTypeFile);
cell_type = CT.cell_type;

%% ---- PRECOMPUTE SDF KERNEL ----
g = gaussianKernelUnitArea_(gaussSigmaMs, dt_ms);

%% ---- GLOBAL WARP TARGETS ----
all_pressLat_global = [];
all_lickLat_global  = [];

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    for k = 1:numel(trials)
        tr = trials(k);
        if requireValid && (~isfield(tr,'valid') || ~tr.valid), continue; end
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

Tpress = median(all_pressLat_global);
Tlick  = median(all_lickLat_global);

%% ---- TEMPLATE AXIS ----
winLeft  = -preCueMs;
winRight_template = Tlick + postLickMs;

tgrid_ms = winLeft:dt_ms:winRight_template;
nT = numel(tgrid_ms);

idx0T = timeToIdx_(0,      winLeft, dt_ms, nT);
idxPT = timeToIdx_(Tpress, winLeft, dt_ms, nT);
idxLT = timeToIdx_(Tlick,  winLeft, dt_ms, nT);

%% ---- BUILD SESSION MATRIX ----
PETH = nan(nT, nSessions);

for sIdx = 1:nSessions
    session = beh(sIdx);
    if ~isGoodTrialsStruct_(session), continue; end
    trials = session.trials;

    spikes_by_ch = getSpikesForSession_(session, S, sIdx);
    if isempty(spikes_by_ch), continue; end

    keepTrial = false(numel(trials),1);
    cueAbs   = nan(numel(trials),1);
    pressLat = nan(numel(trials),1);
    lickLat  = nan(numel(trials),1);

    for k = 1:numel(trials)
        tr = trials(k);
        if requireValid && (~isfield(tr,'valid') || ~tr.valid), continue; end
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

    unitCurves = [];

    for ch = 1:numel(spikes_by_ch)
        uc = spikes_by_ch{ch};
        if ~iscell(uc), continue; end

        for u = 1:numel(uc)
            ct = [];
            if sIdx <= numel(cell_type) && ...
               ch <= numel(cell_type{sIdx}) && ...
               u  <= numel(cell_type{sIdx}{ch})
                ct = cell_type{sIdx}{ch}{u};
            end
            if isempty(ct) || ~ismember(ct,[1 2 3]), continue; end

            spk_abs = double(uc{u}(:));
            if isempty(spk_abs), continue; end

            activeTrials = false(numel(idxKeep),1);
            for iTr = 1:numel(idxKeep)
                cue0 = cueAbsK(iTr);
                if any(spk_abs >= cue0+fixedWin(1) & spk_abs <= cue0+fixedWin(2))
                    activeTrials(iTr) = true;
                end
            end
            if nnz(activeTrials) < minTrialsPerUnit, continue; end

            cueU = cueAbsK(activeTrials);
            rtU  = pressLatK(activeTrials);
            lkU  = lickLatK(activeTrials);

            nTrU = numel(cueU);
            warpedZ = zeros(nTrU, nT);

            for iTr = 1:nTrU
                cue0   = cueU(iTr);
                tPress = rtU(iTr);
                tLick  = lkU(iTr);

                winRight_trial = tLick + postLickMs;
                tgrid_trial = winLeft:dt_ms:winRight_trial;
                nT_trial = numel(tgrid_trial);

                idx0 = timeToIdx_(0, winLeft, dt_ms, nT_trial);
                idxP = timeToIdx_(tPress, winLeft, dt_ms, nT_trial);
                idxL = timeToIdx_(tLick,  winLeft, dt_ms, nT_trial);

                spk_rel = spk_abs - cue0;
                spk_rel = spk_rel(spk_rel >= winLeft & spk_rel <= winRight_trial);

                counts = zeros(1,nT_trial);
                if ~isempty(spk_rel)
                    jj = round((spk_rel - winLeft)/dt_ms)+1;
                    jj = jj(jj>=1 & jj<=nT_trial);
                    for q=1:numel(jj)
                        counts(jj(q)) = counts(jj(q)) + 1;
                    end
                end

                sm = conv(counts,g,'same');
                y  = sm/(dt_ms/1000);

                pre = y(1:idx0);
                segA = y(idx0:idxP);
                wA = warpSegment_floorceilEqn_(segA, max(2,(idxPT-idx0T+1)));

                segB = y(idxP:idxL);
                wB = warpSegment_floorceilEqn_(segB, max(2,(idxLT-idxPT+1)));

                segC = y(idxL:end);
                kC_target = max(1,(nT-idxLT+1));
                if numel(segC) < kC_target
                    segC = [segC zeros(1,kC_target-numel(segC))];
                else
                    segC = segC(1:kC_target);
                end

                yy = [pre wA(2:end) wB(2:end) segC(2:end)];
                if numel(yy)<nT
                    yy = [yy zeros(1,nT-numel(yy))];
                else
                    yy = yy(1:nT);
                end

                warpedZ(iTr,:) = zscoreTrial_(yy);
            end

            idxPlot = 2:2:nTrU;
            if isempty(idxPlot), idxPlot=1:nTrU; end
            muPlot = mean(warpedZ(idxPlot,:),1);

            unitCurves(end+1,:) = muPlot; %#ok<AGROW>
        end
    end

    if isempty(unitCurves), continue; end

    ens  = mean(unitCurves,1);
    ensZ = zscoreTrial_(ens);

    PETH(:,sIdx) = ensZ(:);
end

%% ---- AXIS ----
if xUnit=="s"
    xPlot = tgrid_ms/1000;
    xlab  = 'Warped time from Cue (s)';
    cueLine=0; pressLine=Tpress/1000; lickLine=Tlick/1000;
else
    xPlot = tgrid_ms;
    xlab  = 'Warped time from Cue (ms)';
    cueLine=0; pressLine=Tpress; lickLine=Tlick;
end

%% ---- PLOT ----
fig = figure('Color','w','Position',figPos);
ax = axes(fig); hold(ax,'on');

Pimg = PETH';
hImg = imagesc(ax,xPlot,1:nSessions,Pimg); axis(ax,'tight');
set(ax,'YDir','reverse');
set(hImg,'Interpolation','nearest');

colormap(ax,flipud(gray(256)));
xlabel(ax,xlab,'FontSize',labelFontSize);
ylabel(ax,'Session','FontSize',labelFontSize);
title(ax,'Population activity per session','FontWeight','bold','FontSize',titleFontSize);

xline(ax,cueLine,'--','Color',colCue,'LineWidth',eventLineLW);
xline(ax,pressLine,'--','Color',colPress,'LineWidth',eventLineLW);
xline(ax,lickLine,'--','Color',colLick,'LineWidth',eventLineLW);

caxis(ax,climFixed);

majorY = 1:5:nSessions;
if majorY(end)~=nSessions
    majorY = unique([majorY nSessions]);
end
yticks(ax,majorY);
yticklabels(ax,string(majorY));
ax.YAxis.MinorTick='on';
ax.YAxis.MinorTickValues=setdiff(1:nSessions,majorY);

set(ax,'FontSize',tickFontSize,'TickDir','out','LineWidth',axesTickLW,...
    'TickLength',axesTickLen,'Box','off','Layer','top');

if xUnit=="s" && forceIntegerSecondsTicks
    xl=xlim(ax);
    ts=ceil(xl(1)); te=floor(xl(2));
    if te>=ts, xticks(ax,ts:te); end
end

%% ---- COLORBAR (longer scale, stronger ticks, same position) ----
cb = colorbar(ax,'eastoutside');
cb.Label.String='Z score';
cb.FontSize=cbFontSize;
cb.Label.FontSize=cbFontSize;

cb.LineWidth = 2;
cb.TickLength = 0.02;

axpos = ax.Position;
cbpos = cb.Position;

cbpos(3) = 0.040;                 % thin
cbpos(4) = 0.50 * axpos(4);       % longer scale
cbpos(2) = axpos(2) + (axpos(4)-cbpos(4))/2;
cbpos(1) = axpos(1) + axpos(3) + 0.07*axpos(3);
cbpos(1) = min(cbpos(1), 0.92);   % keep inside figure so label isn't clipped

cb.Position = cbpos;

%% ---- SAVE (PNG + SVG) ----
outbase = fullfile(outDir, sprintf('Sessions01_%02d_sessionEnsemble_heatmap_GLOBALWARP', nSessions));
set(fig, 'PaperPositionMode', 'auto');

% PNG (keep as exportgraphics)
exportgraphics(fig, [outbase '.png'], 'Resolution', 300, 'BackgroundColor', 'white');

% SVG (use saveas/print, since exportgraphics may not support svg on your MATLAB)
try
    saveas(fig, [outbase '.svg']);   % works in many setups
catch
    print(fig, [outbase '.svg'], '-dsvg');  % fallback
end

close(fig);

fprintf('Saved: %s.png\n', outbase);
fprintf('Saved: %s.svg\n', outbase);

%% ================= HELPERS =================
function z = zscoreTrial_(x)
x=x(:)'; mu=mean(x); sd=std(x);
if ~isfinite(sd)||sd<=0, sd=1; end
z=(x-mu)/sd;
end

function g = gaussianKernelUnitArea_(sigmaMs, dtMs)
halfWidth=ceil(5*sigmaMs/dtMs);
x=(-halfWidth:halfWidth)*dtMs;
g=exp(-0.5*(x./sigmaMs).^2);
g=g/sum(g);
end

function idx = timeToIdx_(t_ms, winLeft_ms, dt_ms, nT)
idx=round((t_ms-winLeft_ms)/dt_ms)+1;
idx=max(1,min(nT,idx));
end

function ywarp = warpSegment_floorceilEqn_(y,k)
y=y(:)'; n=numel(y);
if n<2||k<2
ywarp=zeros(1,k); if n>=1&&k>=1,ywarp(1)=y(1); end; return;
end
s=(k-1)/(n-1); ywarp=zeros(1,k);
for i=1:k
x=1+(i-1)/s; x0=floor(x); x1=ceil(x);
x0=max(1,min(n,x0)); x1=max(1,min(n,x1));
if x0==x1
ywarp(i)=y(x0);
else
w1=x-x0; w0=1-w1;
ywarp(i)=w0*y(x0)+w1*y(x1);
end
end
end

function beh = pickBehStruct_(S)
beh=[];
cands={'ratBEHstruct_unit','rat_BEHstruct_unit'};
for k=1:numel(cands)
if isfield(S,cands{k})&&isstruct(S.(cands{k})), beh=S.(cands{k}); return; end
end
f=fieldnames(S);
for i=1:numel(f)
v=S.(f{i});
if isstruct(v)&&numel(v)>1, beh=v; return; end
end
for i=1:numel(f)
v=S.(f{i});
if isstruct(v), beh=v; return; end
end
end

function tf = isGoodTrialsStruct_(session)
tf=false;
if ~isfield(session,'trials')||isempty(session.trials), return; end
if ~isstruct(session.trials), return; end
req={'valid','cue','press','lick'};
tf=all(isfield(session.trials,req));
end

function spikes = getSpikesForSession_(session,S,sessionIdx)
spikes=[];
if isfield(session,'spikes')&&~isempty(session.spikes)
spikes=session.spikes; return
end
if isfield(S,'spikes_session')&&numel(S.spikes_session)>=sessionIdx ...
&&~isempty(S.spikes_session{sessionIdx})
spikes=S.spikes_session{sessionIdx}; return
end
if isfield(S,'spikes_persession')&&numel(S.spikes_persession)>=sessionIdx ...
&&~isempty(S.spikes_persession{sessionIdx})
spikes=S.spikes_persession{sessionIdx}; return
end
if isfield(S,'spikes')&&iscell(S.spikes)&&~isempty(S.spikes)
spikes=S.spikes;
end
end