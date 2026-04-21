function keep_only_sessions_with_at_least_100_cuedNames()
%% Keep only sessions with at least 100 cuedNames in each behavior file.
%
% This script:
%   1) Loads each MAT file
%   2) Finds the main behavior struct (expected: ratBEHstruct)
%   3) Keeps only sessions where numel(cuedNames) >= 100
%   4) Saves a NEW MAT file containing only the filtered behavior struct,
%      using the same struct field name as the input file
%   5) Prints, for each file:
%        - how many sessions were kept
%        - total cuedNames across kept sessions
%
% INPUT FILES:
%   /Volumes/WD_BLACK/AB/B1Bay/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/D6Dahlia/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/F4Fig2/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/J3Jelly/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/J5Joy/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/J7Jasmine/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/L3Lychee/ratBEHstruct.mat
%   /Volumes/WD_BLACK/AB/T8Truffle/ratBEHstruct.mat
%
% OUTPUT DIRECTORY:
%   /Volumes/WD_BLACK/AB/AB
%
% OUTPUT FILES:
%   /Volumes/WD_BLACK/AB/AB/B1Bay_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/D6Dahlia_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/F4Fig2_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/J3Jelly_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/J5Joy_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/J7Jasmine_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/L3Lychee_ratBEHstruct_CUEDNAMES_MIN100.mat
%   /Volumes/WD_BLACK/AB/AB/T8Truffle_ratBEHstruct_CUEDNAMES_MIN100.mat

clear; clc;

%% ---------------- USER SETTINGS ----------------
inFiles = { ...
    '/Volumes/WD_BLACK/AB/B1Bay/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/D6Dahlia/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/F4Fig2/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/J3Jelly/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/J5Joy/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/J7Jasmine/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/L3Lychee/ratBEHstruct.mat', ...
    '/Volumes/WD_BLACK/AB/T8Truffle/ratBEHstruct.mat'};

outDir = '/Volumes/WD_BLACK/AB/AB';
MIN_CUEDNAMES_PER_SESSION = 100;

if exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

fprintf('[INFO] Starting filtering for %d files...\n\n', numel(inFiles));

summary_fileNames    = cell(numel(inFiles),1);
summary_sessionsKept = zeros(numel(inFiles),1);
summary_cuedTotal    = zeros(numel(inFiles),1);

%% ---------------- MAIN LOOP OVER FILES ----------------
for fIdx = 1:numel(inFiles)

    inMat = inFiles{fIdx};
    assert(exist(inMat,'file') == 2, 'File not found: %s', inMat);

    [inDir, inBase, inExt] = fileparts(inMat);
    [~, parentName] = fileparts(inDir);
    outBase = [parentName '_' inBase '_CUEDNAMES_MIN100'];
    outMat = fullfile(outDir, [outBase inExt]);

    fprintf('============================================================\n');
    fprintf('[FILE %d/%d] %s\n', fIdx, numel(inFiles), inMat);
    fprintf('Output: %s\n', outMat);

    S = load(inMat);
    beh = pickBehStruct_(S);
    assert(~isempty(beh) && isstruct(beh), 'Could not find behavior struct in %s', inMat);

    pref = getPreferredFieldName_(S);
    nSessOriginal = numel(beh);

    fprintf('[INFO] Found struct field: %s\n', pref);
    fprintf('[INFO] Sessions before filtering: %d\n', nSessOriginal);

    keepMask = false(1, nSessOriginal);
    nCuedPerSession = zeros(1, nSessOriginal);

    for sIdx = 1:nSessOriginal
        cuedNames = getfield_if_exists_(beh(sIdx), 'cuedNames');

        if isempty(cuedNames)
            nCuedPerSession(sIdx) = 0;
        else
            nCuedPerSession(sIdx) = numel(cuedNames);
        end

        keepMask(sIdx) = nCuedPerSession(sIdx) >= MIN_CUEDNAMES_PER_SESSION;
    end

    nSessKept = sum(keepMask);
    nSessRemoved = sum(~keepMask);

    beh_filtered = beh(keepMask);
    totalCuedNamesKept = sum(nCuedPerSession(keepMask));

    summary_fileNames{fIdx}    = outBase;
    summary_sessionsKept(fIdx) = nSessKept;
    summary_cuedTotal(fIdx)    = totalCuedNamesKept;

    % Save only the filtered behavior struct, using the same field name
    outStruct = struct();
    outStruct.(pref) = beh_filtered;
    save(outMat, '-struct', 'outStruct');

    fprintf('[FILE DONE] %s\n', inBase);
    fprintf('  Sessions before filtering:         %d\n', nSessOriginal);
    fprintf('  Sessions kept (>= %d cuedNames):   %d\n', MIN_CUEDNAMES_PER_SESSION, nSessKept);
    fprintf('  Sessions removed (< %d cuedNames): %d\n', MIN_CUEDNAMES_PER_SESSION, nSessRemoved);
    fprintf('  Total cuedNames kept:              %d\n', totalCuedNamesKept);
    fprintf('  Saved: %s\n\n', outMat);
end

fprintf('============================================================\n');
fprintf('ALL FILES DONE.\n');
fprintf('FINAL SUMMARY PER FILE\n');
fprintf('------------------------------------------------------------\n');
for fIdx = 1:numel(inFiles)
    fprintf('%s\n', summary_fileNames{fIdx});
    fprintf('  Sessions kept: %d\n', summary_sessionsKept(fIdx));
    fprintf('  Total cuedNames: %d\n', summary_cuedTotal(fIdx));
end
fprintf('============================================================\n');

end

%% ========================= HELPERS =========================

function pref = getPreferredFieldName_(S)
    cands = {'ratBEHstruct','ratBEHstruct_unit','rat_BEHstruct_unit'};
    for k = 1:numel(cands)
        if isfield(S, cands{k})
            pref = cands{k};
            return;
        end
    end

    f = fieldnames(S);
    for i = 1:numel(f)
        if isstruct(S.(f{i}))
            pref = f{i};
            return;
        end
    end

    error('No struct field found to write back into.');
end

function beh = pickBehStruct_(S)
    beh = [];
    cands = {'ratBEHstruct','ratBEHstruct_unit','rat_BEHstruct_unit'};

    for k = 1:numel(cands)
        if isfield(S, cands{k}) && isstruct(S.(cands{k}))
            beh = S.(cands{k});
            return;
        end
    end

    f = fieldnames(S);
    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v) && numel(v) > 1
            beh = v;
            return;
        end
    end

    for i = 1:numel(f)
        v = S.(f{i});
        if isstruct(v)
            beh = v;
            return;
        end
    end
end

function v = getfield_if_exists_(s, name)
    if isfield(s, name)
        v = s.(name);
    else
        v = [];
    end
end