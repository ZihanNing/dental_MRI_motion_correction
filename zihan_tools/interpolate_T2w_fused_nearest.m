function interpolate_T2w_fused_nearest(caseList)
% interpolate_T2w_fused_nearest
%
% Go through case folders under /home/zn23/Data/ddMRI/, find NIfTI files
% whose names contain 'Di_fused_nearest_.nii', check whether they are
% T2wSPACE images, and for T2wSPACE images interpolate the smallest matrix
% dimension to double its size, then save with '_interp' appended.
%
% Usage:
%   interpolate_T2w_fused_nearest
%   interpolate_T2w_fused_nearest([1 2 5 10])
%   interpolate_T2w_fused_nearest({'1','2','10'})
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 27-Mar-2026

    basePath = '/home/zn23/Data/ddMRI/';

    if nargin < 1 || isempty(caseList)
        caseFolders = get_numeric_subfolders(basePath);
    else
        caseFolders = normalize_case_list(caseList);
    end

    fprintf('Found %d case(s) to process.\n', numel(caseFolders));

    for iCase = 1:numel(caseFolders)
        caseID = caseFolders{iCase};
        casePath = fullfile(basePath, caseID);

        if ~isfolder(casePath)
            warning('Case folder does not exist: %s. Skipping.', casePath);
            continue;
        end

        fprintf('\n========================================\n');
        fprintf('Processing case: %s\n', caseID);

        niiFile = find_target_nifti(casePath);

        if isempty(niiFile)
            fprintf('No file containing ''Di_fused_nearest_.nii'' found. Skipping.\n');
            continue;
        end

        fprintf('Found file: %s\n', niiFile.name);

        if ~contains(niiFile.name, 'T2wSPACE', 'IgnoreCase', true)
            fprintf('not T2w image, and no interpolations\n');
            continue;
        end

        inputFile = fullfile(casePath, niiFile.name);

        [~, nameOnly, ext] = fileparts(niiFile.name);
        outputFile = fullfile(casePath, [nameOnly '_interp' ext]);

        if exist(outputFile, 'file')
            fprintf('Output already exists, skipping: %s\n', outputFile);
            continue;
        end

        try
            info = niftiinfo(inputFile);
            img  = niftiread(inputFile);

            imgSize = size(img);
            if numel(imgSize) < 3
                warning('Image is not at least 3D. Skipping: %s', inputFile);
                continue;
            end

            spatialSize = imgSize(1:3);
            [~, interpDim] = min(spatialSize);

            fprintf('Original size: [%s]\n', num2str(imgSize));
            fprintf('Interpolating dimension %d from %d to %d\n', ...
                interpDim, spatialSize(interpDim), spatialSize(interpDim) * 2);

            targetSize = imgSize;
            targetSize(interpDim) = targetSize(interpDim) * 2;

            imgInterp = interpolate_along_smallest_dim(img, targetSize);

            save_interp_nifti(inputFile, imgInterp, interpDim, outputFile);
            fprintf('Saved interpolated image to: %s\n', outputFile);

        catch ME
            warning('Failed for case %s: %s', caseID, ME.message);
        end
    end
end


function caseFolders = get_numeric_subfolders(basePath)
    d = dir(basePath);
    isSubDir = [d.isdir];
    names = {d(isSubDir).name};
    names = names(~ismember(names, {'.', '..'}));

    isNumericName = cellfun(@(x) ~isempty(regexp(x, '^\d+$', 'once')), names);
    caseFolders = names(isNumericName);

    caseNums = cellfun(@str2double, caseFolders);
    [~, idx] = sort(caseNums);
    caseFolders = caseFolders(idx);
end


function caseFolders = normalize_case_list(caseList)
    if isnumeric(caseList)
        caseFolders = arrayfun(@num2str, caseList, 'UniformOutput', false);
    elseif iscell(caseList)
        caseFolders = cell(size(caseList));
        for k = 1:numel(caseList)
            if isnumeric(caseList{k})
                caseFolders{k} = num2str(caseList{k});
            else
                caseFolders{k} = char(string(caseList{k}));
            end
        end
    elseif isstring(caseList)
        caseFolders = cellstr(caseList);
    elseif ischar(caseList)
        caseFolders = {caseList};
    else
        error('Unsupported input type for caseList.');
    end
end


function niiFile = find_target_nifti(casePath)
    d = dir(fullfile(casePath, '*.nii'));
    matchIdx = find(contains({d.name}, 'Di_MotCorr.nii', 'IgnoreCase', true), 1, 'first');

    if isempty(matchIdx)
        niiFile = [];
    else
        niiFile = d(matchIdx);
    end
end

function imgInterp = interpolate_along_smallest_dim(img, targetSize)
    if ndims(img) == 3
        imgInterp = imresize3(single(img), targetSize(1:3), 'lanczos3');
    elseif ndims(img) == 4
        nVol = size(img, 4);
        imgInterp = zeros(targetSize, 'single');
        for t = 1:nVol
            imgInterp(:,:,:,t) = imresize3(single(img(:,:,:,t)), targetSize(1:3), 'lanczos3');
        end
    else
        error('Images with ndims > 4 are not supported in this script.');
    end
end


function save_interp_nifti(inputFile, imgInterp, interpDim, outputFile)
% Save interpolated NIfTI in MATLAB only, while preserving geometry.
%
% Strategy:
% 1. Write imgInterp once with default header -> guaranteed consistent size/datatype
% 2. Read that temp header back
% 3. Copy/adjust geometry fields from original header
% 4. Rewrite final output

    % ----- files -----
    [outFolder, outName, ~] = fileparts(outputFile);
    tmpFile = fullfile(outFolder, [outName '_tmp_for_header_fix.nii']);

    % ----- original info -----
    origInfo = niftiinfo(inputFile);

    % ----- step 1: write temp file with default consistent header -----
    niftiwrite(single(imgInterp), tmpFile, 'Compressed', false);

    % ----- step 2: read temp header back -----
    newInfo = niftiinfo(tmpFile);

    % ----- step 3: copy non-size-related metadata from original -----
    % keep size/datatype from newInfo, because they already match imgInterp
    newInfo.Description           = origInfo.Description;
    newInfo.SpaceUnits            = origInfo.SpaceUnits;
    newInfo.TimeUnits             = origInfo.TimeUnits;
    newInfo.AdditiveOffset        = origInfo.AdditiveOffset;
    newInfo.MultiplicativeScaling = origInfo.MultiplicativeScaling;
    newInfo.TimeOffset            = origInfo.TimeOffset;
    newInfo.SliceCode             = origInfo.SliceCode;
    newInfo.FrequencyDimension    = origInfo.FrequencyDimension;
    newInfo.PhaseDimension        = origInfo.PhaseDimension;
    newInfo.SpatialDimension      = origInfo.SpatialDimension;
    newInfo.TransformName         = origInfo.TransformName;
    newInfo.Qfactor               = origInfo.Qfactor;

    % ----- pixel size -----
    newInfo.PixelDimensions = origInfo.PixelDimensions;
    newInfo.PixelDimensions(interpDim) = origInfo.PixelDimensions(interpDim) / 2;

    % ----- transform -----
    % Start from original transform and scale the interpolated axis step by 1/2.
    %
    % For niftiinfo, Transform is returned in a MATLAB transform object.
    % In practice, editing Transform.T together with raw.srow_* and raw.pixdim
    % is the most reliable MATLAB-only approach.
    T = origInfo.Transform.T;

    % For MATLAB image array dims:
    %   dim 1 -> rows    -> Y direction in voxel grid
    %   dim 2 -> cols    -> X direction in voxel grid
    %   dim 3 -> slices  -> Z direction in voxel grid
    %
    % In the affine matrix used by NIfTI, the first three rows of one column
    % describe the physical step for one voxel axis.
    %
    % We therefore halve the step vector of the interpolated axis.
    switch interpDim
        case 1
            T(2,1:3) = T(2,1:3) / 2;
        case 2
            T(1,1:3) = T(1,1:3) / 2;
        case 3
            T(3,1:3) = T(3,1:3) / 2;
        otherwise
            error('interpDim must be 1, 2, or 3.');
    end

    newInfo.Transform.T = T;

    % ----- raw header -----
    % Copy original raw header first, then overwrite fields that must match
    % the new image.
    raw = origInfo.raw;

    % datatype from temp image / single precision
    raw.datatype = newInfo.raw.datatype;
    raw.bitpix   = newInfo.raw.bitpix;

    % image dimensions must match the new image
    imgSize = size(imgInterp);
    raw.dim = newInfo.raw.dim;  % start from valid temp header
    raw.dim(1) = ndims(imgInterp);
    raw.dim(2) = imgSize(1);
    raw.dim(3) = imgSize(2);
    raw.dim(4) = imgSize(3);
    if numel(imgSize) >= 4
        raw.dim(5) = imgSize(4);
    else
        raw.dim(5) = 1;
    end

    % pixel dimensions
    raw.pixdim = newInfo.raw.pixdim;  % start from valid temp header
    raw.pixdim(2) = newInfo.PixelDimensions(1);
    raw.pixdim(3) = newInfo.PixelDimensions(2);
    raw.pixdim(4) = newInfo.PixelDimensions(3);
    raw.pixdim(1) = origInfo.raw.pixdim(1); % qfac

    % Copy original spatial transform codes if present
    if isfield(origInfo.raw, 'qform_code'); raw.qform_code = origInfo.raw.qform_code; end
    if isfield(origInfo.raw, 'sform_code'); raw.sform_code = origInfo.raw.sform_code; end

    % Update sform rows from the edited transform
    if isfield(raw, 'srow_x'); raw.srow_x = T(1, :); end
    if isfield(raw, 'srow_y'); raw.srow_y = T(2, :); end
    if isfield(raw, 'srow_z'); raw.srow_z = T(3, :); end

    % qform quaternion fields are harder to recompute reliably in plain MATLAB.
    % To avoid inconsistency, keep qform code only if original had none;
    % otherwise prefer sform by zeroing qform_code.
    if isfield(raw, 'qform_code')
        raw.qform_code = 0;
    end

    newInfo.raw = raw;

    % ----- step 4: write final output -----
    niftiwrite(single(imgInterp), outputFile, newInfo, 'Compressed', false);

    % remove temp file
    if exist(tmpFile, 'file')
        delete(tmpFile);
    end
end