%% compute_NGS_on_teeth_mask_all_cases.m
% Compute NGS for 4 reconstructions within teeth masks (upper=1, lower=2)
% across all cases under /home/zn23/Data/ddMRI/, then save MAT + Excel.
%
% Assumptions:
% - Each case is a subfolder whose name contains at least one digit.
% - Each case folder contains exactly one teeth mask file matching:
%     *_msk_teeth_fullresol.nii.gz   (found recursively)
% - Each case folder contains 4 reconstruction NIfTI files, one for each recon:
%     NoMoCo_Aq, MoCo_FullFOV, MoCo_UpperJaw, MoCo_LowerJaw
%   The script finds them by searching filenames that CONTAIN the recon key.
%   You MUST edit reconSpecs.searchKey to match your actual filenames.
%
% NGS definition used here:
%   NGS = sum( |grad(I)|^2 within mask ) / sum( I^2 within mask )
% where |grad(I)|^2 = gx^2 + gy^2 + gz^2, computed by central differences.
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-03

clear; clc;

rootDir = '/home/zn23/Data/ddMRI/';

outMat  = fullfile(rootDir, 'NGS_teeth_all_cases.mat');
outXlsx = fullfile(rootDir, 'NGS_teeth_all_cases.xlsx');

% -------------------------------------------------------------------------
% Recon definitions (EDIT searchKey to match your filenames)
% -------------------------------------------------------------------------
reconSpecs = struct( ...
    'name',      {'NoMoCo_Aq',           'MoCo_FullFOV',       'MoCo_UpperJaw',             'MoCo_LowerJaw'}, ...
    'includeKey',{'_Aq_MotCorr.nii',     '_Di_MotCorr.nii',    '_Di_MotCorr_upperjaw_',     '_Di_MotCorr_lowerjaw_'} ...
);

% Image file extensions to consider
imgExts = {'.nii.gz', '.nii'};

% Teeth mask pattern (required)
maskPattern = '*_msk_teeth_fullresol.nii.gz';

% -------------------------------------------------------------------------
% 1) List case folders: subfolders with at least one digit in name
% -------------------------------------------------------------------------
d = dir(rootDir);
isDir = [d.isdir] & ~startsWith({d.name}, '.');
names = {d(isDir).name};

hasDigit = cellfun(@(s) ~isempty(regexp(s, '\d', 'once')), names);
caseNames = names(hasDigit);

if isempty(caseNames)
    error('No case folders (name contains digit) found under %s', rootDir);
end

% Sort cases by numeric part (robust)
caseNums = nan(size(caseNames));
for i = 1:numel(caseNames)
    m = regexp(caseNames{i}, '\d+', 'match');
    if ~isempty(m)
        caseNums(i) = str2double(m{1});
    end
end
[~, ord] = sortrows([isnan(caseNums(:)), caseNums(:)], [1 2]);
caseNames = caseNames(ord);

% -------------------------------------------------------------------------
% 2) Loop cases: find mask + recon images, compute NGS in upper/lower teeth
% -------------------------------------------------------------------------
rows = {};   % cell of structs
allVarNames = {};

for c = 1:numel(caseNames)
    caseID = caseNames{c};
    caseDir = fullfile(rootDir, caseID);

    fprintf('\n[%s] Processing...\n', caseID);

    % --- find teeth mask (recursive)
    mf = dir(fullfile(caseDir, '**', maskPattern));
    if isempty(mf)
        warning('[%s] Teeth mask not found (%s). Skipping.', caseID, maskPattern);
        continue;
    end
    % pick newest if multiple
    [~, midx] = max([mf.datenum]);
    maskPath = fullfile(mf(midx).folder, mf(midx).name);
    fprintf('[%s] Mask: %s\n', caseID, maskPath);

    % read mask
    try
        msk = readNiiAny(maskPath);
    catch ME
        warning('[%s] Failed to read mask: %s\n%s', caseID, maskPath, ME.message);
        continue;
    end

    upperMask = (msk == 1);
    lowerMask = (msk == 2);
    
    % Expand masks by 2 voxels in 3D
    se = strel('sphere', 2);   % 3D spherical structuring element, radius = 2 voxels

    upperMask = imdilate(upperMask, se);
    lowerMask = imdilate(lowerMask, se);
    

    if ~any(upperMask(:))
        warning('[%s] Upper teeth mask is empty (no voxels==1).', caseID);
    end
    if ~any(lowerMask(:))
        warning('[%s] Lower teeth mask is empty (no voxels==2).', caseID);
    end

    % Prepare row struct
    row = struct();
    row.ID = string(caseID);

    % --- for each recon: find image and compute NGS on upper/lower masks
    for r = 1:numel(reconSpecs)
        reconName = reconSpecs(r).name;
        key = reconSpecs(r).includeKey;

        imgPath = findReconNifti(caseDir, reconSpecs(r).includeKey, reconSpecs(r).name, imgExts);
        if isempty(imgPath)
            warning('[%s] Recon image not found for %s (key="%s").', caseID, reconName, key);
            row.(matlab.lang.makeValidName([reconName '_UpperTeeth'])) = NaN;
            row.(matlab.lang.makeValidName([reconName '_LowerTeeth'])) = NaN;
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_UpperTeeth']); %#ok<SAGROW>
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_LowerTeeth']); %#ok<SAGROW>
            continue;
        end

        fprintf('[%s] Recon %-12s: %s\n', caseID, reconName, imgPath);

        try
            I = readNiiAny(imgPath);
        catch ME
            warning('[%s] Failed to read image: %s\n%s', caseID, imgPath, ME.message);
            row.(matlab.lang.makeValidName([reconName '_UpperTeeth'])) = NaN;
            row.(matlab.lang.makeValidName([reconName '_LowerTeeth'])) = NaN;
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_UpperTeeth']); %#ok<SAGROW>
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_LowerTeeth']); %#ok<SAGROW>
            continue;
        end

        if ~isequal(size(I), size(msk))
            warning('[%s] Image/mask size mismatch for %s. image=%s mask=%s. Skipping recon.', ...
                caseID, reconName, mat2str(size(I)), mat2str(size(msk)));
            row.(matlab.lang.makeValidName([reconName '_UpperTeeth'])) = NaN;
            row.(matlab.lang.makeValidName([reconName '_LowerTeeth'])) = NaN;
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_UpperTeeth']); %#ok<SAGROW>
            allVarNames{end+1} = matlab.lang.makeValidName([reconName '_LowerTeeth']); %#ok<SAGROW>
            continue;
        end

        % Compute NGS
        ngsUpper = computeNGS(I, upperMask);
        ngsLower = computeNGS(I, lowerMask);

        fU = matlab.lang.makeValidName([reconName '_UpperTeeth']);
        fL = matlab.lang.makeValidName([reconName '_LowerTeeth']);
        row.(fU) = ngsUpper;
        row.(fL) = ngsLower;

        allVarNames{end+1} = fU; %#ok<SAGROW>
        allVarNames{end+1} = fL; %#ok<SAGROW>
    end

    rows{end+1} = row; %#ok<SAGROW>
end

if isempty(rows)
    error('No cases processed successfully.');
end

allVarNames = unique(allVarNames, 'stable');

% -------------------------------------------------------------------------
% 3) Build final table, sort by numeric ID, save MAT + Excel
% -------------------------------------------------------------------------
n = numel(rows);

T_teeth = table('Size', [n, 1 + numel(allVarNames)], ...
    'VariableTypes', [{'string'}, repmat({'double'}, 1, numel(allVarNames))], ...
    'VariableNames', [{'ID'}, allVarNames]);

for i = 1:n
    rs = rows{i};
    T_teeth.ID(i) = rs.ID;

    for v = 1:numel(allVarNames)
        fn = allVarNames{v};
        if isfield(rs, fn)
            T_teeth{i, 1+v} = rs.(fn);
        else
            T_teeth{i, 1+v} = NaN;
        end
    end
end

% Sort by numeric part of ID
idNums = nan(height(T_teeth), 1);
for i = 1:height(T_teeth)
    m = regexp(T_teeth.ID(i), '\d+', 'match');
    if ~isempty(m), idNums(i) = str2double(m{1}); end
end
[~, ord] = sortrows([isnan(idNums), idNums], [1 2]);
T_teeth = T_teeth(ord, :);

save(outMat, 'T_teeth', 'rootDir', 'reconSpecs', 'maskPattern');
writetable(T_teeth, outXlsx, 'FileType', 'spreadsheet');

fprintf('\nDone.\nSaved MAT : %s\nSaved XLSX: %s\n', outMat, outXlsx);

% =========================================================================
% Local functions
% =========================================================================
function vol = readNiiAny(p)
% Read .nii or .nii.gz using built-in niftiread.
% For .nii.gz, MATLAB supports it in newer versions; if not, gunzip to temp.

    p = char(p);
    [~, ~, ext] = fileparts(p);

    if strcmpi(ext, '.gz')
        % Handle .nii.gz
        tmpDir = tempname;
        mkdir(tmpDir);
        gunzip(p, tmpDir);
        niiFiles = dir(fullfile(tmpDir, '*.nii'));
        if isempty(niiFiles)
            error('gunzip succeeded but no .nii found in temp folder.');
        end
        niiPath = fullfile(niiFiles(1).folder, niiFiles(1).name);
        vol = niftiread(niiPath);
        % Cleanup
        try
            rmdir(tmpDir, 's');
        catch
        end
    else
        vol = niftiread(p);
    end

    vol = double(vol);
end

function imgPath = findReconNifti(caseDir, includeKey, reconName, imgExts)
% Find a recon NIfTI file under caseDir recursively whose filename contains includeKey.
% Special handling:
% - For MoCo_FullFOV, exclude files that contain 'upperjaw' or 'lowerjaw' to avoid mismatches.
% Picks the newest match if multiple.

    imgPath = '';
    hits = [];

    % Search both .nii and .nii.gz
    for e = 1:numel(imgExts)
        ext = imgExts{e};
        ff = dir(fullfile(caseDir, '**', ['*' ext]));
        hits = [hits; ff]; %#ok<AGROW>
    end
    if isempty(hits), return; end

    % Filter by includeKey in filename
    keep = false(numel(hits),1);
    for i = 1:numel(hits)
        fn = hits(i).name;
        keep(i) = contains(fn, includeKey);
    end
    hits = hits(keep);
    if isempty(hits), return; end

    % Extra exclusion for FullFOV: avoid matching upper/lower variants
    if strcmp(reconName, 'MoCo_FullFOV')
        keep2 = true(numel(hits),1);
        for i = 1:numel(hits)
            fn = lower(hits(i).name);
            if contains(fn, 'upperjaw') || contains(fn, 'lowerjaw')
                keep2(i) = false;
            end
        end
        hits = hits(keep2);
        if isempty(hits), return; end
    end

    % Pick newest
    [~, idx] = max([hits.datenum]);
    imgPath = fullfile(hits(idx).folder, hits(idx).name);
end

function ngs = computeNGS(I, mask)
% Compute NGS within mask:
%   NGS = sum( |grad(I)|^2 ) / sum( I^2 )
% Return NaN if mask empty or denominator == 0.

    if ~any(mask(:))
        ngs = NaN;
        return;
    end

    % Replace NaNs if any
    I(~isfinite(I)) = 0;

    % Central-difference gradients (simple and fast)
    gx = gradientCentral(I, 1);
    gy = gradientCentral(I, 2);
    gz = gradientCentral(I, 3);

    g2 = gx.^2 + gy.^2 + gz.^2;

    num = sum(g2(mask), 'omitnan');
    den = sum((I(mask)).^2, 'omitnan');

    if den <= 0 || ~isfinite(den)
        ngs = NaN;
    else
        ngs = num / den;
    end
end

function g = gradientCentral(I, dim)
% Central difference along dimension dim, with forward/backward at edges.

    sz = size(I);
    g = zeros(sz);

    idxAll = repmat({':'}, 1, ndims(I));

    % interior: (I(i+1)-I(i-1))/2
    idx1 = idxAll; idx2 = idxAll; idx3 = idxAll;
    idx1{dim} = 2:sz(dim)-1;
    idx2{dim} = 3:sz(dim);
    idx3{dim} = 1:sz(dim)-2;
    g(idx1{:}) = (I(idx2{:}) - I(idx3{:})) / 2;

    % first slice: forward difference
    idxF = idxAll; idxF2 = idxAll;
    idxF{dim}  = 1;
    idxF2{dim} = 2;
    g(idxF{:}) = I(idxF2{:}) - I(idxF{:});

    % last slice: backward difference
    idxL = idxAll; idxL2 = idxAll;
    idxL{dim}  = sz(dim);
    idxL2{dim} = sz(dim)-1;
    g(idxL{:}) = I(idxL{:}) - I(idxL2{:});
end