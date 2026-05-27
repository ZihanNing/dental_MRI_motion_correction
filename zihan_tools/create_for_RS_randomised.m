%% create_for_RS_randomised.m
% Create anonymised/randomised image folders for radiologist scoring.
%
% Workflow:
% 1. Search /home/zn23/Data/ddMRI/ for case subfolders (folder names containing digits)
% 2. Create /home/zn23/Data/ddMRI/for_RS
% 3. For each case, find 3 images:
%       - No-MoCo       : filename contains '_Aq_MotCorr.nii'
%       - MoCo-fullFOV  : filename contains '_Di_MotCorr.nii'
%                         but NOT upperjaw/lowerjaw/fused
%       - MoCo-fused    : filename contains '_Di_fused_.nii'
% 4. Copy them into /for_RS/<caseID>/
% 5. Randomly rename them as A/B/C (preserving extension)
% 6. Save the mapping table for later recovery
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-11

% This is in use 24-Mar-2026 by Zihan

clear; clc;

rootDir = '/home/zn23/Data/ddMRI/';
outDir  = fullfile(rootDir, 'for_RS');

mappingXlsx = fullfile(outDir, 'RS_randomisation_key.xlsx');
mappingMat  = fullfile(outDir, 'RS_randomisation_key.mat');

rng('shuffle');   % randomise order each time script runs

% -------------------------------------------------------------------------
% Find case folders
% -------------------------------------------------------------------------
d = dir(rootDir);
isDir = [d.isdir] & ~startsWith({d.name}, '.');
names = {d(isDir).name};

hasDigit = cellfun(@(s) ~isempty(regexp(s, '\d', 'once')), names);
caseNamesAll = names(hasDigit);

if isempty(caseNamesAll)
    error('No case folders found under %s', rootDir);
end

% Sort by numeric part if available
caseNumsAll = nan(size(caseNamesAll));
for i = 1:numel(caseNamesAll)
    m = regexp(caseNamesAll{i}, '\d+', 'match');
    if ~isempty(m)
        caseNumsAll(i) = str2double(m{1});
    end
end
[~, ord] = sortrows([isnan(caseNumsAll(:)), caseNumsAll(:)], [1 2]);
caseNamesAll = caseNamesAll(ord);

% -------------------------------------------------------------------------
% Create output folder
% -------------------------------------------------------------------------
if ~isfolder(outDir)
    mkdir(outDir);
end

fprintf('Found %d case(s).\n', numel(caseNamesAll));
fprintf('Output folder: %s\n', outDir);

% -------------------------------------------------------------------------
% Prepare output mapping
% -------------------------------------------------------------------------
rows = {};

for c = 1:numel(caseNamesAll)
    caseID = caseNamesAll{c};
    caseDir = fullfile(rootDir, caseID);

    fprintf('\n[%s] Processing...\n', caseID);

    % Find the three source files
    p_noMoCo  = findReconNifti(caseDir, '_Aq_MotCorr.nii', 'NoMoCo');
    p_fullFOV = findReconNifti(caseDir, '_Di_MotCorr.nii', 'MoCoFullFOV');
    p_fused   = findReconNifti(caseDir, '_Di_fused_nearest_.nii', 'MoCoFused');

    if isempty(p_noMoCo)
        warning('[%s] No-MoCo image not found. Skipping case.', caseID);
        continue;
    end
    if isempty(p_fullFOV)
        warning('[%s] MoCo-fullFOV image not found. Skipping case.', caseID);
        continue;
    end
    if isempty(p_fused)
        warning('[%s] MoCo-fused image not found. Skipping case.', caseID);
        continue;
    end

    % Make case output folder
    caseOutDir = fullfile(outDir, caseID);
    if ~isfolder(caseOutDir)
        mkdir(caseOutDir);
    end

    % Random assignment of labels A/B/C
    srcPaths = {p_noMoCo, p_fullFOV, p_fused};
    srcCats  = {'No-MoCo', 'MoCo-fullFOV', 'MoCo-fused'};
    anonLabs = {'A', 'B', 'C'};

    perm = randperm(3);
    srcPaths = srcPaths(perm);
    srcCats  = srcCats(perm);

    % Copy and rename
    anonFiles = cell(1,3);
    for k = 1:3
        src = srcPaths{k};
        lab = anonLabs{k};

        extOut = getNiftiExtension(src);  % '.nii' or '.nii.gz'
        dst = fullfile(caseOutDir, [lab extOut]);

        copyfile(src, dst);
        anonFiles{k} = dst;

        fprintf('  %s -> %s (%s)\n', srcCats{k}, dst, lab);
    end

    % Save row
    row = struct();
    row.ID = string(caseID);

    % lab A/B/C corresponds to srcCats{1/2/3} after permutation
    row.A_category = string(srcCats{1});
    row.B_category = string(srcCats{2});
    row.C_category = string(srcCats{3});

    row.A_file = string(anonFiles{1});
    row.B_file = string(anonFiles{2});
    row.C_file = string(anonFiles{3});

    row.NoMoCo_source      = string(p_noMoCo);
    row.MoCoFullFOV_source = string(p_fullFOV);
    row.MoCoFused_source   = string(p_fused);

    rows{end+1} = row; %#ok<SAGROW>
end

if isempty(rows)
    error('No cases were processed successfully.');
end

% -------------------------------------------------------------------------
% Build table
% -------------------------------------------------------------------------
T_key = struct2table([rows{:}]', 'AsArray', true);

% Sort by numeric part of ID
idNums = nan(height(T_key), 1);
for i = 1:height(T_key)
    m = regexp(char(T_key.ID(i)), '\d+', 'match');
    if ~isempty(m)
        idNums(i) = str2double(m{1});
    end
end
[~, ord] = sortrows([isnan(idNums), idNums], [1 2]);
T_key = T_key(ord, :);

% Save full key
save(mappingMat, 'T_key');
writetable(T_key, mappingXlsx);

fprintf('\nDone.\n');
fprintf('Saved key MAT : %s\n', mappingMat);
fprintf('Saved key XLSX: %s\n', mappingXlsx);

% Optional: save a simplified sheet for checking
simpleXlsx = fullfile(outDir, 'RS_randomisation_key_simple.xlsx');
T_simple = T_key(:, {'ID','A_category','B_category','C_category'});
writetable(T_simple, simpleXlsx);
fprintf('Saved simple key: %s\n', simpleXlsx);

% =========================================================================
% Local functions
% =========================================================================
function imgPath = findReconNifti(caseDir, includeKey, modeName)
% Find recon file recursively.
% Rules:
% - NoMoCo: filename contains '_Aq_MotCorr.nii'
% - MoCoFullFOV: filename contains '_Di_MotCorr.nii' but excludes upper/lower/fused
% - MoCoFused: filename contains '_Di_fused_.nii'
% If multiple matches exist, choose the newest one.

    allHits = [dir(fullfile(caseDir, '**', '*.nii')); ...
               dir(fullfile(caseDir, '**', '*.nii.gz'))];

    if isempty(allHits)
        imgPath = '';
        return;
    end

    keep = false(numel(allHits), 1);
    for i = 1:numel(allHits)
        fn = allHits(i).name;
        fnLow = lower(fn);

        switch modeName
            case 'NoMoCo'
                keep(i) = contains(fn, includeKey);

            case 'MoCoFullFOV'
                keep(i) = contains(fn, includeKey) && ...
                          ~contains(fnLow, 'upperjaw') && ...
                          ~contains(fnLow, 'lowerjaw') && ...
                          ~contains(fnLow, 'fused');

            case 'MoCoFused'
                keep(i) = contains(fn, includeKey);

            otherwise
                error('Unknown modeName: %s', modeName);
        end
    end

    hits = allHits(keep);

    if isempty(hits)
        imgPath = '';
        return;
    end

    [~, idx] = max([hits.datenum]);
    imgPath = fullfile(hits(idx).folder, hits(idx).name);
end

function ext = getNiftiExtension(p)
% Return '.nii' or '.nii.gz'

    if endsWith(p, '.nii.gz', 'IgnoreCase', true)
        ext = '.nii.gz';
    elseif endsWith(p, '.nii', 'IgnoreCase', true)
        ext = '.nii';
    else
        error('Not a NIfTI file: %s', p);
    end
end