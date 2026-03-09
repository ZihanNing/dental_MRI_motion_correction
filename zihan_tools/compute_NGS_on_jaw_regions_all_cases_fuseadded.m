%% compute_NGS_on_jaw_regions_all_cases.m
% Compute NGS for multiple reconstructions within:
%   1) upper jaw region  (defined by idx1:idx2 along HF dimension)
%   2) lower jaw region  (defined by idx2:idx3 along HF dimension)
%   3) full FOV head mask (from *msk_head_fullresol.nii.gz)
% across selected cases under /home/zn23/Data/ddMRI/, then save MAT + CSV.
%
% NGS definition:
%   NGS = sum( |grad(I)|^2 within mask ) / sum( I^2 within mask )
% where |grad(I)|^2 = gx^2 + gy^2 + gz^2, computed by central differences.
%
% NOTE:
% - Upper/lower jaw regions are NOT computed from the small teeth mask itself.
% - Instead, teeth mask is only used to derive anatomical landmarks:
%       idx1 = superior boundary of upper jaw
%       idx2 = inter-arch boundary
%       idx3 = inferior boundary of lower jaw
% - Then larger slab masks are created using fillMaskBetween().
%
% Required external helper functions on MATLAB path:
%   compute_landmarks_from_teeth_mask
%   readLocationTxt
%   fillMaskBetween
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-09

% this computation is different from the batch processing, DO NOT USE FOR
% COMPARISON


clear; clc;

rootDir = '/home/zn23/Data/ddMRI/';

outMat = fullfile(rootDir, 'NGS_all_cases.mat');
outCsv = fullfile(rootDir, 'NGS_all_cases.csv');

selectedCases = [];   % [] means all

% -------------------------------------------------------------------------
% User settings
% -------------------------------------------------------------------------
hfDim = 1;   % Foot-Head / Superior-Inferior dimension in the image volume

% Recon definitions
reconSpecs = struct( ...
    'name',      {'NoMoCo_Aq',           'MoCo_FullFOV',       'MoCo_UpperJaw',             'MoCo_LowerJaw',             'MoCo_Fused'}, ...
    'includeKey',{'_Aq_MotCorr.nii',     '_Di_MotCorr.nii',    '_Di_MotCorr_upperjaw_',     '_Di_MotCorr_lowerjaw_',     'Di_fused'} ...
);

% Image file extensions to consider
imgExts = {'.nii.gz', '.nii'};

% Required mask patterns
teethMaskPattern = '*_msk_teeth_fullresol.nii.gz';
headMaskPattern  = '*msk_head_fullresol.nii.gz';

% -------------------------------------------------------------------------
% 1) List case folders: subfolders with at least one digit in name
% -------------------------------------------------------------------------
d = dir(rootDir);
isDir = [d.isdir] & ~startsWith({d.name}, '.');
names = {d(isDir).name};

hasDigit = cellfun(@(s) ~isempty(regexp(s, '\d', 'once')), names);
caseNamesAll = names(hasDigit);

if isempty(caseNamesAll)
    error('No case folders (name contains digit) found under %s', rootDir);
end

% Sort cases by numeric part
caseNumsAll = nan(size(caseNamesAll));
for i = 1:numel(caseNamesAll)
    m = regexp(caseNamesAll{i}, '\d+', 'match');
    if ~isempty(m)
        caseNumsAll(i) = str2double(m{1});
    end
end
[~, ord] = sortrows([isnan(caseNumsAll(:)), caseNumsAll(:)], [1 2]);
caseNamesAll = caseNamesAll(ord);
caseNumsAll  = caseNumsAll(ord);

% Optional subset selection
if isempty(selectedCases)
    caseNames = caseNamesAll;
else
    keep = ismember(caseNumsAll, selectedCases);
    caseNames = caseNamesAll(keep);

    if isempty(caseNames)
        error('None of the selectedCases were found under %s', rootDir);
    end
end

fprintf('Processing %d case(s).\n', numel(caseNames));
disp(caseNames(:));

% -------------------------------------------------------------------------
% 2) Loop cases: find masks + recon images, compute NGS
% -------------------------------------------------------------------------
rows = {};
allVarNames = {};

for c = 1:numel(caseNames)
    caseID = caseNames{c};
    caseDir = fullfile(rootDir, caseID);

    fprintf('\n[%s] Processing...\n', caseID);

    % --- find teeth mask (recursive)
    mfTeeth = dir(fullfile(caseDir, '**', teethMaskPattern));
    if isempty(mfTeeth)
        warning('[%s] Teeth mask not found (%s). Skipping case.', caseID, teethMaskPattern);
        continue;
    end
    [~, midxTeeth] = max([mfTeeth.datenum]);
    teethMaskPath = fullfile(mfTeeth(midxTeeth).folder, mfTeeth(midxTeeth).name);
    fprintf('[%s] Teeth mask: %s\n', caseID, teethMaskPath);

    % --- find head mask (recursive)
    mfHead = dir(fullfile(caseDir, '**', headMaskPattern));
    if isempty(mfHead)
        warning('[%s] Head mask not found (%s). FullFOV NGS will be NaN.', caseID, headMaskPattern);
        headMaskPath = '';
    else
        [~, midxHead] = max([mfHead.datenum]);
        headMaskPath = fullfile(mfHead(midxHead).folder, mfHead(midxHead).name);
        fprintf('[%s] Head mask : %s\n', caseID, headMaskPath);
    end

    % --- read teeth mask
    try
        mskTeeth = readNiiAny(teethMaskPath);
    catch ME
        warning('[%s] Failed to read teeth mask: %s\n%s', caseID, teethMaskPath, ME.message);
        continue;
    end

    % --- derive location.txt and jaw-region landmarks
    locFile = fullfile(caseDir, 'location.txt');
    try
        [idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid, locFile] = ...
            compute_landmarks_from_teeth_mask(teethMaskPath, locFile); %#ok<ASGLU>

        [idx1, idx2, idx3] = readLocationTxt(locFile);
    catch ME
        warning('[%s] Failed to compute/read jaw landmarks from teeth mask.\n%s', caseID, ME.message);
        continue;
    end

    fprintf('[%s] Jaw landmarks from %s: idx1=%d, idx2=%d, idx3=%d\n', ...
        caseID, locFile, idx1, idx2, idx3);

    % --- validate jaw landmark order
    if ~(isscalar(idx1) && isscalar(idx2) && isscalar(idx3) && ...
            isfinite(idx1) && isfinite(idx2) && isfinite(idx3))
        warning('[%s] Invalid jaw landmark indices. Skipping case.', caseID);
        continue;
    end

    idx1 = round(idx1);
    idx2 = round(idx2);
    idx3 = round(idx3);

    if ~(idx1 >= idx2 && idx2 >= idx3)
        warning('[%s] Unexpected landmark order: idx1=%d, idx2=%d, idx3=%d. Skipping case.', ...
            caseID, idx1, idx2, idx3);
        continue;
    end

    % --- read head mask if available
    headMask = [];
    if ~isempty(headMaskPath)
        try
            mskHead = readNiiAny(headMaskPath);
            headMask = (mskHead > 0);
            if ~any(headMask(:))
                warning('[%s] Head mask is empty.', caseID);
            end
        catch ME
            warning('[%s] Failed to read head mask: %s\n%s', caseID, headMaskPath, ME.message);
            headMask = [];
        end
    end

    % Prepare row struct
    row = struct();
    row.ID = string(caseID);
    row.idx1 = idx1;
    row.idx2 = idx2;
    row.idx3 = idx3;

    allVarNames{end+1} = 'idx1'; %#ok<SAGROW>
    allVarNames{end+1} = 'idx2'; %#ok<SAGROW>
    allVarNames{end+1} = 'idx3'; %#ok<SAGROW>

    % --- for each recon: find image and compute NGS on upper/lower/fullFOV
    for r = 1:numel(reconSpecs)
        reconName = reconSpecs(r).name;
        key = reconSpecs(r).includeKey;

        imgPath = findReconNifti(caseDir, key, reconName, imgExts);

        fU = matlab.lang.makeValidName([reconName '_UpperJaw']);
        fL = matlab.lang.makeValidName([reconName '_LowerJaw']);
        fF = matlab.lang.makeValidName([reconName '_FullFOV']);

        if isempty(imgPath)
            warning('[%s] Recon image not found for %s (key="%s").', caseID, reconName, key);
            row.(fU) = NaN;
            row.(fL) = NaN;
            row.(fF) = NaN;
            allVarNames{end+1} = fU; %#ok<SAGROW>
            allVarNames{end+1} = fL; %#ok<SAGROW>
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        fprintf('[%s] Recon %-12s: %s\n', caseID, reconName, imgPath);

        try
            I = readNiiAny(imgPath);
        catch ME
            warning('[%s] Failed to read image: %s\n%s', caseID, imgPath, ME.message);
            row.(fU) = NaN;
            row.(fL) = NaN;
            row.(fF) = NaN;
            allVarNames{end+1} = fU; %#ok<SAGROW>
            allVarNames{end+1} = fL; %#ok<SAGROW>
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        % Size checks against masks
        if ~isequal(size(I), size(mskTeeth))
            warning('[%s] Image/teeth-mask size mismatch for %s. image=%s teethMask=%s. Skipping recon.', ...
                caseID, reconName, mat2str(size(I)), mat2str(size(mskTeeth)));
            row.(fU) = NaN;
            row.(fL) = NaN;
            row.(fF) = NaN;
            allVarNames{end+1} = fU; %#ok<SAGROW>
            allVarNames{end+1} = fL; %#ok<SAGROW>
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        % Clamp indices to image size along HF dimension
        szI = size(I);
        idx1c = max(1, min(szI(hfDim), idx1));
        idx2c = max(1, min(szI(hfDim), idx2));
        idx3c = max(1, min(szI(hfDim), idx3));

        if ~(idx1c >= idx2c && idx2c >= idx3c)
            warning('[%s] Clamped landmark order invalid for %s. FullFOV only.', caseID, reconName);
            maskUpper = false(szI);
            maskLower = false(szI);
        else
            maskUpper = false(szI);
            maskLower = false(szI);

            maskUpper = fillMaskBetween(maskUpper, hfDim, idx2c, idx1c);
            maskLower = fillMaskBetween(maskLower, hfDim, idx3c, idx2c);
        end

        if ~any(maskUpper(:))
            warning('[%s] Upper jaw mask is empty for %s.', caseID, reconName);
        end
        if ~any(maskLower(:))
            warning('[%s] Lower jaw mask is empty for %s.', caseID, reconName);
        end

        % Optional: restrict jaw slabs to head mask
        % This avoids background voxels outside the head region contributing to NGS
        if ~isempty(headMask)
            if isequal(size(headMask), size(I))
                maskUpper = maskUpper & headMask;
                maskLower = maskLower & headMask;
            else
                warning('[%s] Image/head-mask size mismatch for %s. Jaw masks not intersected with head mask.', ...
                    caseID, reconName);
            end
        end

        % Compute NGS: upper/lower jaw regions
        ngsUpper = compute_ngs(I, maskUpper);
        ngsLower = compute_ngs(I, maskLower);

        % Compute NGS: full FOV (head mask)
        if isempty(headMask)
            ngsFull = NaN;
        elseif ~isequal(size(I), size(headMask))
            warning('[%s] Image/head-mask size mismatch for %s. image=%s headMask=%s. FullFOV set to NaN.', ...
                caseID, reconName, mat2str(size(I)), mat2str(size(headMask)));
            ngsFull = NaN;
        else
            ngsFull = compute_ngs(I, headMask);
        end

        row.(fU) = ngsUpper;
        row.(fL) = ngsLower;
        row.(fF) = ngsFull;

        allVarNames{end+1} = fU; %#ok<SAGROW>
        allVarNames{end+1} = fL; %#ok<SAGROW>
        allVarNames{end+1} = fF; %#ok<SAGROW>
    end

    rows{end+1} = row; %#ok<SAGROW>

    % Optional: save per-case MAT
    caseOutMat = fullfile(caseDir, [caseID '_NGS_all_regions.mat']);
    NGS_case = struct2table(row, 'AsArray', true);
    save(caseOutMat, 'row', 'NGS_case', 'teethMaskPath', 'headMaskPath', ...
        'reconSpecs', 'locFile', 'hfDim');
    fprintf('[%s] Saved per-case MAT: %s\n', caseID, caseOutMat);
end

if isempty(rows)
    error('No cases processed successfully.');
end

allVarNames = unique(allVarNames, 'stable');

% -------------------------------------------------------------------------
% 3) Build final table, sort by numeric ID
% -------------------------------------------------------------------------
n = numel(rows);

varTypes = cell(1, 1 + numel(allVarNames));
varTypes{1} = 'string';
for i = 2:numel(varTypes)
    varTypes{i} = 'double';
end

T_NGS = table('Size', [n, 1 + numel(allVarNames)], ...
    'VariableTypes', varTypes, ...
    'VariableNames', [{'ID'}, allVarNames]);

for i = 1:n
    rs = rows{i};
    T_NGS.ID(i) = rs.ID;

    for v = 1:numel(allVarNames)
        fn = allVarNames{v};
        if isfield(rs, fn)
            T_NGS{i, 1+v} = rs.(fn);
        else
            T_NGS{i, 1+v} = NaN;
        end
    end
end

% Sort by numeric part of ID
idNums = nan(height(T_NGS), 1);
for i = 1:height(T_NGS)
    m = regexp(char(T_NGS.ID(i)), '\d+', 'match');
    if ~isempty(m)
        idNums(i) = str2double(m{1});
    end
end
[~, ord] = sortrows([isnan(idNums), idNums], [1 2]);
T_NGS = T_NGS(ord, :);

% -------------------------------------------------------------------------
% 4) Save MAT + CSV
% -------------------------------------------------------------------------
save(outMat, 'T_NGS', 'rootDir', 'reconSpecs', 'teethMaskPattern', ...
    'headMaskPattern', 'selectedCases', 'hfDim');
writetable(T_NGS, outCsv);

fprintf('\nDone.\nSaved MAT : %s\nSaved CSV : %s\n', outMat, outCsv);

% =========================================================================
% Local functions
% =========================================================================
function vol = readNiiAny(p)
% Read .nii or .nii.gz using built-in niftiread.
% For .nii.gz, gunzip to temp if needed.

    p = char(p);
    [~, ~, ext] = fileparts(p);

    if strcmpi(ext, '.gz')
        tmpDir = tempname;
        mkdir(tmpDir);
        gunzip(p, tmpDir);
        niiFiles = dir(fullfile(tmpDir, '*.nii'));
        if isempty(niiFiles)
            error('gunzip succeeded but no .nii found in temp folder.');
        end
        niiPath = fullfile(niiFiles(1).folder, niiFiles(1).name);
        vol = niftiread(niiPath);
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
% - For MoCo_FullFOV, exclude files that contain upper/lower/fused variants.
% Picks the newest match if multiple.

    imgPath = '';
    hits = [];

    for e = 1:numel(imgExts)
        ext = imgExts{e};
        ff = dir(fullfile(caseDir, '**', ['*' ext]));
        hits = [hits; ff]; %#ok<AGROW>
    end
    if isempty(hits), return; end

    keep = false(numel(hits),1);
    for i = 1:numel(hits)
        fn = hits(i).name;
        keep(i) = contains(fn, includeKey);
    end
    hits = hits(keep);
    if isempty(hits), return; end

    if strcmp(reconName, 'MoCo_FullFOV')
        keep2 = true(numel(hits),1);
        for i = 1:numel(hits)
            fn = lower(hits(i).name);
            if contains(fn, 'upperjaw') || contains(fn, 'lowerjaw') || contains(fn, 'fused')
                keep2(i) = false;
            end
        end
        hits = hits(keep2);
        if isempty(hits), return; end
    end

    [~, idx] = max([hits.datenum]);
    imgPath = fullfile(hits(idx).folder, hits(idx).name);
end

function g = gradientCentral(I, dim)
% Central difference along dimension dim, with forward/backward at edges.

    sz = size(I);
    g = zeros(sz);

    if sz(dim) < 2
        return;
    end

    idxAll = repmat({':'}, 1, ndims(I));

    if sz(dim) > 2
        idx1 = idxAll; idx2 = idxAll; idx3 = idxAll;
        idx1{dim} = 2:sz(dim)-1;
        idx2{dim} = 3:sz(dim);
        idx3{dim} = 1:sz(dim)-2;
        g(idx1{:}) = (I(idx2{:}) - I(idx3{:})) / 2;
    end

    idxF = idxAll; idxF2 = idxAll;
    idxF{dim}  = 1;
    idxF2{dim} = min(2, sz(dim));
    g(idxF{:}) = I(idxF2{:}) - I(idxF{:});

    idxL = idxAll; idxL2 = idxAll;
    idxL{dim}  = sz(dim);
    idxL2{dim} = max(sz(dim)-1, 1);
    g(idxL{:}) = I(idxL{:}) - I(idxL2{:});
end
function [idx1, idx2, idx3] = readLocationTxt(locFile)
    vals = readmatrix(locFile, 'FileType', 'text');
    vals = vals(~isnan(vals));
    if numel(vals) < 3
        error('location.txt does not contain 3 valid indices: %s', locFile);
    end
    idx1 = round(vals(1));
    idx2 = round(vals(2));
    idx3 = round(vals(3));
end
function mask = fillMaskBetween(mask, dim, idxA, idxB)
% Fill mask between two indices (inclusive), regardless of order.
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-09

    lo = min(idxA, idxB);
    hi = max(idxA, idxB);

    sz = size(mask);
    lo = max(1, min(sz(dim), round(lo)));
    hi = max(1, min(sz(dim), round(hi)));

    idx = repmat({':'}, 1, ndims(mask));
    idx{dim} = lo:hi;
    mask(idx{:}) = true;
end