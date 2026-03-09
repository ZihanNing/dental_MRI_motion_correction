%% compute_NGS_fullFOV_all_cases.m
% Compute NGS in the FullFOV region only for five reconstructions:
%   1) NoMoCo_Aq
%   2) MoCo_FullFOV
%   3) MoCo_UpperJaw
%   4) MoCo_LowerJaw
%   5) MoCo_Fused
%
% The folder traversal follows the style of the previous all-case script:
%   - iterate case folders under rootDir
%   - find reconstructions recursively
%   - find head mask recursively
%   - find teeth mask recursively only to derive location.txt / idxHeadTop
%
% IMPORTANT:
% The FullFOV ROI is defined exactly as in script 2:
%   roi_full = headMask & hfInFull3
% where
%   fullHF_lo = 1;
%   fullHF_hi = idxHeadTop_txt;
%
% Thus, for the four reconstructions already computed in script 2, the
% FullFOV NGS values should match, assuming the same files are selected.
%
% This script uses external compute_ngs(I, mask, [vx vy vz]).
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-09

clear; clc;

rootDir = '/home/zn23/Data/ddMRI/';

outMat = fullfile(rootDir, 'NGS_fullFOV_all_cases.mat');
outCsv = fullfile(rootDir, 'NGS_fullFOV_all_cases.csv');

selectedCases = [];   % [] means all

% -------------------------------------------------------------------------
% Recon definitions
% -------------------------------------------------------------------------
reconSpecs = struct( ...
    'name',      {'NoMoCo_Aq',           'MoCo_FullFOV',       'MoCo_UpperJaw',             'MoCo_LowerJaw',             'MoCo_Fused'}, ...
    'includeKey',{'_Aq_MotCorr.nii',     '_Di_MotCorr.nii',    '_Di_MotCorr_upperjaw_',     '_Di_MotCorr_lowerjaw_',     'Di_fused'} ...
);

imgExts = {'.nii.gz', '.nii'};

teethMaskPattern = '*_msk_teeth_fullresol.nii.gz';
headMaskPattern  = '*msk_head_fullresol.nii.gz';

% -------------------------------------------------------------------------
% 1) List case folders
% -------------------------------------------------------------------------
d = dir(rootDir);
isDir = [d.isdir] & ~startsWith({d.name}, '.');
names = {d(isDir).name};

hasDigit = cellfun(@(s) ~isempty(regexp(s, '\d', 'once')), names);
caseNamesAll = names(hasDigit);

if isempty(caseNamesAll)
    error('No case folders (name contains digit) found under %s', rootDir);
end

% Sort by numeric part
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
% 2) Loop over cases
% -------------------------------------------------------------------------
rows = {};
allVarNames = {};

for c = 1:numel(caseNames)
    caseID = caseNames{c};
    caseDir = fullfile(rootDir, caseID);

    fprintf('\n[%s] Processing...\n', caseID);

    % --- find teeth mask (only for deriving location.txt / idxHeadTop)
    mfTeeth = dir(fullfile(caseDir, '**', teethMaskPattern));
    if isempty(mfTeeth)
        warning('[%s] Teeth mask not found (%s). Skipping case.', caseID, teethMaskPattern);
        continue;
    end
    [~, midxTeeth] = max([mfTeeth.datenum]);
    teethMaskPath = fullfile(mfTeeth(midxTeeth).folder, mfTeeth(midxTeeth).name);
    fprintf('[%s] Teeth mask: %s\n', caseID, teethMaskPath);

    % --- compute/read location.txt so that idxHeadTop is available
    locFile = fullfile(caseDir, 'location.txt');
    try
        [idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid, locFile] = ...
            compute_landmarks_from_teeth_mask(teethMaskPath, locFile); %#ok<NASGU,ASGLU>
        vals = readLocationTxt(locFile);
        idxHeadTop_txt = round(vals(1));
    catch ME
        warning('[%s] Failed to compute/read location.txt from teeth mask.\n%s', caseID, ME.message);
        continue;
    end

    % --- find head mask
    mfHead = dir(fullfile(caseDir, '**', headMaskPattern));
    if isempty(mfHead)
        warning('[%s] Head mask not found (%s). Skipping case.', caseID, headMaskPattern);
        continue;
    end
    [~, midxHead] = max([mfHead.datenum]);
    headMaskPath = fullfile(mfHead(midxHead).folder, mfHead(midxHead).name);
    fprintf('[%s] Head mask : %s\n', caseID, headMaskPath);

    try
        headInfo = niftiinfo(headMaskPath);
        headMask = niftiread(headInfo);
        headMask = (headMask > 0);
    catch ME
        warning('[%s] Failed to read head mask: %s\n%s', caseID, headMaskPath, ME.message);
        continue;
    end

    if ~any(headMask(:))
        warning('[%s] Head mask is empty. Skipping case.', caseID);
        continue;
    end

    % Prepare row
    row = struct();
    row.ID = string(caseID);
    row.idxHeadTop = idxHeadTop_txt;

    allVarNames{end+1} = 'idxHeadTop'; %#ok<SAGROW>

    % --- compute FullFOV NGS for each recon
    for r = 1:numel(reconSpecs)
        reconName = reconSpecs(r).name;
        key = reconSpecs(r).includeKey;

        fF = matlab.lang.makeValidName([reconName '_FullFOV']);

        imgPath = findReconNifti(caseDir, key, reconName, imgExts);
        if isempty(imgPath)
            warning('[%s] Recon image not found for %s (key="%s").', caseID, reconName, key);
            row.(fF) = NaN;
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        fprintf('[%s] Recon %-12s: %s\n', caseID, reconName, imgPath);

        try
            infoI = niftiinfo(imgPath);
            I = single(niftiread(infoI));
        catch ME
            warning('[%s] Failed to read recon image: %s\n%s', caseID, imgPath, ME.message);
            row.(fF) = NaN;
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        sz = size(I);
        if ~isequal(size(headMask), sz)
            warning('[%s] Image/head-mask size mismatch for %s. image=%s headMask=%s. FullFOV set to NaN.', ...
                caseID, reconName, mat2str(size(I)), mat2str(size(headMask)));
            row.(fF) = NaN;
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        % -----------------------------------------------------------------
        % FullFOV ROI definition: EXACTLY as script 2
        % -----------------------------------------------------------------
        T = infoI.Transform.T;
        A = T(1:3,1:3);

        worldAxisOfVoxel = zeros(1,3); % 1=X(RL), 2=Y(AP), 3=Z(HF)
        for vdim = 1:3
            [~, w] = max(abs(A(:,vdim)));
            worldAxisOfVoxel(vdim) = w;
        end

        dimHF = find(worldAxisOfVoxel == 3, 1);
        if isempty(dimHF)
            warning('[%s] Failed to infer HF dimension for %s. FullFOV set to NaN.', caseID, reconName);
            row.(fF) = NaN;
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        HFmax = sz(dimHF);
        idxHeadTop_txt_clamped = max(1, min(HFmax, idxHeadTop_txt));

        fullHF_lo = 1;
        fullHF_hi = idxHeadTop_txt_clamped;

        hfIdx = 1:HFmax;
        hfInFull = (hfIdx >= fullHF_lo & hfIdx <= fullHF_hi);

        shape = ones(1,3);
        shape(dimHF) = HFmax;
        hfInFull3 = reshape(hfInFull, shape);

        roi_full = headMask & hfInFull3;

        if ~any(roi_full(:))
            warning('[%s] FullFOV ROI is empty for %s. FullFOV set to NaN.', caseID, reconName);
            row.(fF) = NaN;
            allVarNames{end+1} = fF; %#ok<SAGROW>
            continue;
        end

        vx = infoI.PixelDimensions(1);
        vy = infoI.PixelDimensions(2);
        vz = infoI.PixelDimensions(3);

        ngsFull = compute_ngs(I, roi_full, [vx vy vz]);

        row.(fF) = ngsFull;
        allVarNames{end+1} = fF; %#ok<SAGROW>
    end

    rows{end+1} = row; %#ok<SAGROW>

    % Optional per-case save
    caseOutMat = fullfile(caseDir, [caseID '_NGS_fullFOV.mat']);
    NGS_case = struct2table(row, 'AsArray', true);
    save(caseOutMat, 'row', 'NGS_case', 'teethMaskPath', 'headMaskPath', 'locFile', 'reconSpecs');
    fprintf('[%s] Saved per-case MAT: %s\n', caseID, caseOutMat);
end

if isempty(rows)
    error('No cases processed successfully.');
end

allVarNames = unique(allVarNames, 'stable');

% -------------------------------------------------------------------------
% 3) Build final table
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
save(outMat, 'T_NGS', 'rootDir', 'reconSpecs', 'teethMaskPattern', 'headMaskPattern', 'selectedCases');
writetable(T_NGS, outCsv);

fprintf('\nDone.\nSaved MAT : %s\nSaved CSV : %s\n', outMat, outCsv);
fprintf('\nVariables in final table:\n');
disp(T_NGS.Properties.VariableNames');

% =========================================================================
% Local functions
% =========================================================================
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

function vals = readLocationTxt(locFile)
% Read numeric values from location.txt

    vals = readmatrix(locFile, 'FileType', 'text');
    vals = vals(~isnan(vals));

    if numel(vals) < 3
        error('location.txt does not contain at least 3 valid indices: %s', locFile);
    end
end