function build_trials_ALLSESS_FINAL_NO_LICK_BOUND()
%% Build strict cue->press->lick trials and add 'trials' to ALL sessions
% in /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_FINAL.mat
%
% CHANGE REQUESTED:
%   - REMOVE any lick "bound" rule (e.g., lick must occur before next cue/press).
%     Here we implement: LickCand = first lick AFTER press (no upper bound).
%
% TRIAL STRUCTURE BY SESSION:
%   - sessions 1–38:    cue-press-lick
%   - sessions 39–41:   cue1-press1-cue2-press2-lick
%   - sessions 42–120:  cue1-press1-cue2-press2-cue3-press3-lick
%
% COMMAND WINDOW:
%   - Print ALL rejection reasons per session.
%
% NEW ADDITION (ONLY ADDITION):
%   - Split output into two datasets (same top-level beh struct format):
%       * beh_1_38  : sessions 1–38
%       * beh_39_120: sessions 39–end (typically 120)
%     and save them as separate MAT files (in addition to the original outMat).
%
% OUTPUT:
%   - Full: /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_TESTFINAL_WITH_TRIALS.mat
%   - Split: /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_TESTFINAL_WITH_TRIALS__SESS01_38.mat
%   - Split: /Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_s1_TESTFINAL_WITH_TRIALS__SESS39_120.mat
%
% ADDITION REQUESTED NOW:
%   - Add trial fields for cue type and poke type (R/C/L), analogous to "correct"
%       * cueType         : from cuedNames (first character)
%       * pokeType        : from pokeNames (first character; e.g., "CL" -> "C")
%       * pokeNameRaw     : full original pokeNames cell string (e.g., "CL")
%       * pokeNameIsMulti : true if pokeNameRaw has multiple letters (numel>1)

clear; clc;

%% ---- USER SETTINGS ----
inMat  = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL.mat';
outMat = '/Volumes/WD_BLACK/A_THESIS_FINAL/rat_BEHstruct_unit_FINAL_WITH_TRIALS.mat';

% (NEW) split outputs
[outDir, outBase, ~] = fileparts(outMat);
outMat_1_36   = fullfile(outDir, [outBase '__SESS01_36.mat']);
outMat_37_151 = fullfile(outDir, [outBase '__SESS37_151.mat']);

% trial window padding (ms, same units as timestamps in your dataset)
winCue  = [-1000, 3000];  % used for t_start = cue + winCue(1)
winLick = [-1000, 3000];  % used for t_end   = lick + winLick(2)

% minimum reaction time (press must be at least this long after cue)
MIN_RT_MS = 100;

%% ---- LOAD / PICK STRUCT ----
assert(exist(inMat,'file')==2, 'File not found: %s', inMat);
S   = load(inMat);
beh = pickBehStruct_(S);
assert(~isempty(beh) && isstruct(beh),'Could not find behavior struct');

nSess = numel(beh);
fprintf('[INFO] Building trials for %d sessions... (MIN_RT_MS = %d)\n', nSess, MIN_RT_MS);

%% ---- MAIN ----
nBuilt = 0;
nSkippedMissingFields = 0;
nNoValidTrials = 0;

for sIdx = 1:nSess
    session = beh(sIdx);

    cueCells  = getfield_if_exists_(session,'cuedTimes');
    pokeTimesCellPerTrial = getfield_if_exists_(session,'pokeTimes');
    lickGlobal = getVecFromMaybeCell_( getfield_if_exists_(session,'LickOn') );

    % Names aligned 1:1 with cuedTimes (may not exist)
    cuedNames = getfield_if_exists_(session,'cuedNames');
    pokeNames = getfield_if_exists_(session,'pokeNames');

    % Determine cue/press count per trial by session index
    if sIdx <= 38
        nCuesPerTrial = 1;
    elseif sIdx <= 41
        nCuesPerTrial = 2;
    else
        nCuesPerTrial = 3;
    end

    % Require core inputs
    if isempty(cueCells) || isempty(pokeTimesCellPerTrial) || isempty(lickGlobal)
        fprintf('[sess %d] Missing cuedTimes/pokeTimes/LickOn -> no trials built.\n', sIdx);

        if nCuesPerTrial == 1
            beh(sIdx).trials = struct('cue',{},'press',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        elseif nCuesPerTrial == 2
            beh(sIdx).trials = struct('cue1',{},'press1',{},'cue2',{},'press2',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        else
            beh(sIdx).trials = struct('cue1',{},'press1',{},'cue2',{},'press2',{},'cue3',{},'press3',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        end

        nSkippedMissingFields = nSkippedMissingFields + 1;
        continue;
    end

    nTrialsDeclared = numel(cueCells);

    % Preallocate prototype struct with correct fields
    if nCuesPerTrial == 1
        trials = struct('cue',[],'press',[],'lick',[], ...
                        't_start',[],'t_end',[],'valid',false, ...
                        'correct',[], ...
                        'cueType','', 'pokeType','', 'pokeNameRaw','', 'pokeNameIsMulti',false);
    elseif nCuesPerTrial == 2
        trials = struct('cue1',[],'press1',[],'cue2',[],'press2',[],'lick',[], ...
                        't_start',[],'t_end',[],'valid',false, ...
                        'correct',[], ...
                        'cueType','', 'pokeType','', 'pokeNameRaw','', 'pokeNameIsMulti',false);
    else
        trials = struct('cue1',[],'press1',[],'cue2',[],'press2',[],'cue3',[],'press3',[],'lick',[], ...
                        't_start',[],'t_end',[],'valid',false, ...
                        'correct',[], ...
                        'cueType','', 'pokeType','', 'pokeNameRaw','', 'pokeNameIsMulti',false);
    end

    ti = 0;

    % ---- rejection diagnostics counters (per session) ----
    rej_nanCue            = 0;
    rej_incompleteBlock   = 0; % not enough cues left to form a full multi-cue trial
    rej_noPress           = 0;
    rej_pressTooFast      = 0; % (press-cue) < MIN_RT_MS
    rej_pressAfterNextCue = 0; % press not before within-trial next cue OR not before nextCue bound
    rej_noLickAfterPress  = 0;

    % For printing properly for multi-cue trials:
    nBlocksAttempted = 0;
    nBlocksRejected  = 0;

    % Iterate cues in blocks (1, 2, or 3 cues per "trial")
    t = 1;
    while t <= nTrialsDeclared

        % Ensure we have enough cues to form the block
        if (t + nCuesPerTrial - 1) > nTrialsDeclared
            rej_incompleteBlock = rej_incompleteBlock + 1;
            break;
        end

        nBlocksAttempted = nBlocksAttempted + 1;

        % ---- cue(s) ----
        cue1 = safeScalar_(cueCells, t);
        if isnan(cue1)
            rej_nanCue = rej_nanCue + 1;
            nBlocksRejected = nBlocksRejected + 1;
            t = t + nCuesPerTrial;
            continue;
        end

        if nCuesPerTrial >= 2
            cue2 = safeScalar_(cueCells, t+1);
            if isnan(cue2)
                rej_nanCue = rej_nanCue + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
        end

        if nCuesPerTrial >= 3
            cue3 = safeScalar_(cueCells, t+2);
            if isnan(cue3)
                rej_nanCue = rej_nanCue + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
        end

        % Next cue bound (used ONLY to ensure LAST press belongs to this trial)
        nextCue = inf;
        nextCueIdx = t + nCuesPerTrial;
        if nextCueIdx <= nTrialsDeclared
            nc = safeScalar_(cueCells, nextCueIdx);
            if ~isnan(nc), nextCue = nc; end
        end

        % ---- press(es) ----
        press1 = firstPressAfterCue_(pokeTimesCellPerTrial, t, cue1);
        if isnan(press1)
            rej_noPress = rej_noPress + 1;
            nBlocksRejected = nBlocksRejected + 1;
            t = t + nCuesPerTrial;
            continue;
        end
        if (press1 - cue1) < MIN_RT_MS
            rej_pressTooFast = rej_pressTooFast + 1;
            nBlocksRejected = nBlocksRejected + 1;
            t = t + nCuesPerTrial;
            continue;
        end
        if nCuesPerTrial >= 2
            if ~(press1 < cue2)
                rej_pressAfterNextCue = rej_pressAfterNextCue + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
        else
            if ~(press1 < nextCue)
                rej_pressAfterNextCue = rej_pressAfterNextCue + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
        end

        if nCuesPerTrial >= 2
            press2 = firstPressAfterCue_(pokeTimesCellPerTrial, t+1, cue2);
            if isnan(press2)
                rej_noPress = rej_noPress + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
            if (press2 - cue2) < MIN_RT_MS
                rej_pressTooFast = rej_pressTooFast + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
            if nCuesPerTrial >= 3
                if ~(press2 < cue3)
                    rej_pressAfterNextCue = rej_pressAfterNextCue + 1;
                    nBlocksRejected = nBlocksRejected + 1;
                    t = t + nCuesPerTrial;
                    continue;
                end
            else
                if ~(press2 < nextCue)
                    rej_pressAfterNextCue = rej_pressAfterNextCue + 1;
                    nBlocksRejected = nBlocksRejected + 1;
                    t = t + nCuesPerTrial;
                    continue;
                end
            end
        end

        if nCuesPerTrial >= 3
            press3 = firstPressAfterCue_(pokeTimesCellPerTrial, t+2, cue3);
            if isnan(press3)
                rej_noPress = rej_noPress + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
            if (press3 - cue3) < MIN_RT_MS
                rej_pressTooFast = rej_pressTooFast + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
            if ~(press3 < nextCue)
                rej_pressAfterNextCue = rej_pressAfterNextCue + 1;
                nBlocksRejected = nBlocksRejected + 1;
                t = t + nCuesPerTrial;
                continue;
            end
        end

        % ---- lick after FINAL press (NO UPPER BOUND) ----
        if nCuesPerTrial == 1
            lastPress = press1;
        elseif nCuesPerTrial == 2
            lastPress = press2;
        else
            lastPress = press3;
        end

        li = find(lickGlobal > lastPress, 1, 'first');
        if isempty(li)
            rej_noLickAfterPress = rej_noLickAfterPress + 1;
            nBlocksRejected = nBlocksRejected + 1;
            t = t + nCuesPerTrial;
            continue;
        end
        lickCand = lickGlobal(li);

        % correctness flag + cue/poke type fields from cuedNames/pokeNames (if available)
        corrFlag = NaN;

        cueType = '';
        pokeType = '';
        pokeNameRaw = '';
        pokeNameIsMulti = false;

        if ~isempty(cuedNames)
            cName = safeLabel_(cuedNames, t);
            cName = upper(strtrim(cName));
            if ~isempty(cName)
                cueType = cName(1);
                if ~ismember(cueType, ['R','C','L'])
                    cueType = '';
                end
            end
        end

        if ~isempty(pokeNames)
            pName = safeLabel_(pokeNames, t);
            pName = upper(strtrim(pName));
            pokeNameRaw = pName;
            if ~isempty(pName)
                pokeType = pName(1); % if e.g. "CL", take "C"
                if ~ismember(pokeType, ['R','C','L'])
                    pokeType = '';
                end
                pokeNameIsMulti = (numel(pName) > 1);
            end
        end

        if ~isempty(cuedNames) && ~isempty(pokeNames)
            cName2 = safeLabel_(cuedNames, t);
            pName2 = safeLabel_(pokeNames, t);
            if ~isempty(cName2) && ~isempty(pName2)
                corrFlag = double(strcmp(cName2, pName2)); % 1 if same, else 0
            end
        end

        % accept
        ti = ti + 1;

        if nCuesPerTrial == 1
            trials(ti).cue     = cue1;
            trials(ti).press   = press1;
            trials(ti).lick    = lickCand;
            trials(ti).t_start = cue1     + winCue(1);
            trials(ti).t_end   = lickCand + winLick(2);
            trials(ti).valid   = true;
            trials(ti).correct = corrFlag;

            trials(ti).cueType         = cueType;
            trials(ti).pokeType        = pokeType;
            trials(ti).pokeNameRaw     = pokeNameRaw;
            trials(ti).pokeNameIsMulti = pokeNameIsMulti;

        elseif nCuesPerTrial == 2
            trials(ti).cue1    = cue1;
            trials(ti).press1  = press1;
            trials(ti).cue2    = cue2;
            trials(ti).press2  = press2;
            trials(ti).lick    = lickCand;
            trials(ti).t_start = cue1     + winCue(1);
            trials(ti).t_end   = lickCand + winLick(2);
            trials(ti).valid   = true;
            trials(ti).correct = corrFlag;

            trials(ti).cueType         = cueType;
            trials(ti).pokeType        = pokeType;
            trials(ti).pokeNameRaw     = pokeNameRaw;
            trials(ti).pokeNameIsMulti = pokeNameIsMulti;

        else
            trials(ti).cue1    = cue1;
            trials(ti).press1  = press1;
            trials(ti).cue2    = cue2;
            trials(ti).press2  = press2;
            trials(ti).cue3    = cue3;
            trials(ti).press3  = press3;
            trials(ti).lick    = lickCand;
            trials(ti).t_start = cue1     + winCue(1);
            trials(ti).t_end   = lickCand + winLick(2);
            trials(ti).valid   = true;
            trials(ti).correct = corrFlag;

            trials(ti).cueType         = cueType;
            trials(ti).pokeType        = pokeType;
            trials(ti).pokeNameRaw     = pokeNameRaw;
            trials(ti).pokeNameIsMulti = pokeNameIsMulti;
        end

        t = t + nCuesPerTrial;
    end

    if ti == 0
        fprintf('[sess %d] No valid trials constructed. (declared %d cues)\n', sIdx, nTrialsDeclared);

        if nCuesPerTrial == 1
            beh(sIdx).trials = struct('cue',{},'press',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        elseif nCuesPerTrial == 2
            beh(sIdx).trials = struct('cue1',{},'press1',{},'cue2',{},'press2',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        else
            beh(sIdx).trials = struct('cue1',{},'press1',{},'cue2',{},'press2',{},'cue3',{},'press3',{},'lick',{}, ...
                                      't_start',{},'t_end',{},'valid',{}, ...
                                      'correct',{}, ...
                                      'cueType',{},'pokeType',{},'pokeNameRaw',{},'pokeNameIsMulti',{});
        end

        nNoValidTrials = nNoValidTrials + 1;

        fprintf('  Rejections (blocks attempted=%d, accepted=%d):\n', nBlocksAttempted, ti);
        fprintf('    Incomplete cue block at end:       %d\n', rej_incompleteBlock);
        fprintf('    NaN/missing cue:                   %d\n', rej_nanCue);
        fprintf('    No press after cue:                %d\n', rej_noPress);
        fprintf('    Press < %d ms after cue:           %d\n', MIN_RT_MS, rej_pressTooFast);
        fprintf('    Press not before next cue/bound:   %d\n', rej_pressAfterNextCue);
        fprintf('    No lick after final press:         %d\n', rej_noLickAfterPress);

    else
        beh(sIdx).trials = trials;
        nBuilt = nBuilt + 1;

        fprintf('[sess %d] Built %d trials (declared %d cues)\n', sIdx, ti, nTrialsDeclared);

        fprintf('  Rejections (blocks attempted=%d, accepted=%d, rejected=%d):\n', ...
            nBlocksAttempted, ti, nBlocksRejected);
        fprintf('    Incomplete cue block at end:       %d\n', rej_incompleteBlock);
        fprintf('    NaN/missing cue:                   %d\n', rej_nanCue);
        fprintf('    No press after cue:                %d\n', rej_noPress);
        fprintf('    Press < %d ms after cue:           %d\n', MIN_RT_MS, rej_pressTooFast);
        fprintf('    Press not before next cue/bound:   %d\n', rej_pressAfterNextCue);
        fprintf('    No lick after final press:         %d\n', rej_noLickAfterPress);

        sumRej = rej_incompleteBlock + rej_nanCue + rej_noPress + rej_pressTooFast + rej_pressAfterNextCue + rej_noLickAfterPress;
        if (sumRej ~= nBlocksRejected) && (nBlocksAttempted > 0)
            fprintf('    (Note) Uncounted rejects:          %d\n', nBlocksRejected - sumRej);
        end
    end
end

%% ---- SAVE (FULL DATASET) ----
S_out = S;
pref = getPreferredFieldName_(S_out);
S_out.(pref) = beh;
save(outMat, '-struct', 'S_out', '-v7.3');

%% ---- NEW: SPLIT INTO TWO DATASETS (NO OTHER CHANGES) ----
beh_1_36   = beh(1:min(36, numel(beh)));
if numel(beh) >= 37
    beh_37_151 = beh(37:end);
else
    beh_37_151 = beh([]); %#ok<NBRAK> % empty same type
end

S_out_1_36 = S_out;
S_out_1_36.(pref) = beh_1_36;
save(outMat_1_36, '-struct', 'S_out_1_36', '-v7.3');

S_out_37_151 = S_out;
S_out_37_151.(pref) = beh_37_151;
save(outMat_37_151, '-struct', 'S_out_37_151', '-v7.3');

fprintf('\nDONE.\n');
fprintf('Sessions w/ trials built: %d / %d\n', nBuilt, nSess);
fprintf('Sessions missing core fields: %d\n', nSkippedMissingFields);
fprintf('Sessions w/ 0 valid trials: %d\n', nNoValidTrials);
fprintf('Saved (full):  %s\n', outMat);
fprintf('Saved (1-38):  %s\n', outMat_1_36);
fprintf('Saved (39+):   %s\n', outMat_37_151);

end % function build_trials_ALLSESS_FINAL_NO_LICK_BOUND

%% ----------------- HELPERS -----------------
function pref = getPreferredFieldName_(S)
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k=1:numel(cands)
        if isfield(S,cands{k}), pref = cands{k}; return; end
    end
    f = fieldnames(S);
    for i=1:numel(f)
        if isstruct(S.(f{i})), pref = f{i}; return; end
    end
    error('No struct field to write back into.');
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k=1:numel(cands)
        if isfield(S,cands{k}) && isstruct(S.(cands{k}))
            beh = S.(cands{k}); return;
        end
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

function v = getfield_if_exists_(s, name)
    if isfield(s,name), v = s.(name); else, v = []; end
end

function vec = getVecFromMaybeCell_(x)
    vec = [];
    if isempty(x), return; end
    if isnumeric(x), vec = double(x(:)); return; end
    if iscell(x)
        tmp = nan(numel(x),1);
        for i=1:numel(x), tmp(i) = safeScalar_(x,i); end
        vec = double(tmp(~isnan(tmp))); return;
    end
end

function v = safeScalar_(cellA, idx)
    v = NaN;
    if ~iscell(cellA) || idx<1 || idx>numel(cellA), return; end
    xi = cellA{idx};
    if isempty(xi), return; end
    if isnumeric(xi)
        xi = double(xi(:)); if ~isempty(xi), v = xi(1); end
    elseif isstring(xi) || ischar(xi)
        vv = str2double(string(xi)); if ~isnan(vv), v = vv; end
    end
end

function lbl = safeLabel_(cellA, idx)
    lbl = '';
    if isempty(cellA) || ~iscell(cellA) || idx<1 || idx>numel(cellA), return; end
    xi = cellA{idx};
    if isempty(xi), return; end
    if isstring(xi)
        xi = char(xi);
    elseif isnumeric(xi)
        xi = char(string(xi(1)));
    end
    if ischar(xi)
        lbl = strtrim(xi);
    end
end

function press = firstPressAfterCue_(pokeTimesCellPerTrial, t, cue_t)
    press = NaN;
    if t<1 || t>numel(pokeTimesCellPerTrial), return; end
    pt = pokeTimesCellPerTrial{t};
    if isempty(pt), return; end
    pt = double(pt(:));
    j = find(pt > cue_t, 1, 'first');
    if ~isempty(j), press = pt(j); end
end