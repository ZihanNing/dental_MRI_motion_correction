%% Fuse MoCo-upper and MoCo-lower into one image via lower-jaw registration
% This script:
% 1) finds I_U and I_L in each numeric case folder
% 2) reads idx_1, idx_2, idx_3 from location.txt
% 3) builds upper/lower masks along the head-foot dimension
% 4) registers I_U to I_L in the lower-jaw region to get T_{U->L}
% 5) inverts the transform to get T_{L->U}
% 6) warps I_L into I_U space
% 7) builds feathering weights w(x) from distance transforms
% 8) computes the fused image
% 9) saves the fused image as a NIfTI in the same case folder
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-09

% This is in use 24-Mar-2026 by Zihan

clear; clc;

%% User settings
rootDir = '/home/zn23/Data/ddMRI/';
hfDim = 1;                 % head-foot dimension in the NIfTI array
overlapHalfWidth = 5;      % voxels; controls overlap band for feathering
saveAsSingle = true;       % save fused image as single precision
verbose = true;

%% Find numeric case folders
d = dir(rootDir);
isCase = [d.isdir] & ~startsWith({d.name}, '.') & ...
         cellfun(@(x) ~isempty(regexp(x, '^\d+$', 'once')), {d.name});
caseDirs = d(isCase);

fprintf('Found %d numeric case folders under %s\n', numel(caseDirs), rootDir);

for iCase = 1:numel(caseDirs)
    caseName = caseDirs(iCase).name;
    caseFolder = fullfile(rootDir, caseName);

    fprintf('\n==============================\n');
    fprintf('Processing case: %s\n', caseName);
    fprintf('Folder: %s\n', caseFolder);

    try
        %% ----------------------------------------------------------------
        % 1. Find I_U and I_L
        %% ----------------------------------------------------------------
        upperFile = findSingleFile(caseFolder, '*_Di_MotCorr_upperjaw_*.nii*');
        lowerFile = findSingleFile(caseFolder, '*_Di_MotCorr_lowerjaw_*.nii*');

        if isempty(upperFile)
            warning('Case %s: MoCo-upper file not found. Skipping.', caseName);
            continue;
        end
        if isempty(lowerFile)
            warning('Case %s: MoCo-lower file not found. Skipping.', caseName);
            continue;
        end

        if verbose
            fprintf('I_U: %s\n', upperFile);
            fprintf('I_L: %s\n', lowerFile);
        end

        %% ----------------------------------------------------------------
        % 2. Read images and location.txt
        %% ----------------------------------------------------------------
        locFile = fullfile(caseFolder, 'location_.txt');

        % If location.txt does not exist, try to generate it from teeth mask
        if ~exist(locFile, 'file')
            fprintf('Case %s: location.txt not found. Trying to generate from teeth mask...\n', caseName);

            teethMaskFile = findLatestTeethMask(caseFolder);

            if isempty(teethMaskFile)
                warning('Case %s: neither location.txt nor *msk_teeth_fullresol.nii.gz found. Skipping.', caseName);
                continue;
            end

            try
                [idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid, locFile] = ...
                    compute_landmarks_from_teeth_mask(teethMaskFile, locFile);

                fprintf('Landmarks saved to: %s\n', locFile);
                fprintf('HF indices: head=%d, lips=%d, chin=%d (RL mid=%d)\n', ...
                    idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid);

            catch ME
                warning('Case %s: failed to generate location.txt from teeth mask: %s', ...
                    caseName, ME.message);
                continue;
            end
        end

        infoU = niftiinfo(upperFile);
        IU = double(niftiread(infoU));

        infoL = niftiinfo(lowerFile);
        IL = double(niftiread(infoL));

        if ~isequal(size(IU), size(IL))
            warning('Case %s: I_U and I_L have different sizes. Skipping.', caseName);
            continue;
        end

        [idx1, idx2, idx3] = readLocationTxt(locFile);

        %% ----------------------------------------------------------------
        % 3. Build upper/lower hard masks and expanded masks
        %% ----------------------------------------------------------------
        volSize = size(IU);

        maskUpper = false(volSize);
        maskLower = false(volSize);

        maskUpper = fillMaskBetween(maskUpper, hfDim, idx1, idx2);
        maskLower = fillMaskBetween(maskLower, hfDim, idx2, idx3);

        % Expanded masks for feathering
        maskUpperExp = dilateAlongDim(maskUpper, hfDim, overlapHalfWidth);
        maskLowerExp = dilateAlongDim(maskLower, hfDim, overlapHalfWidth);

        % Keep fusion local: outside the union, default to I_U
        unionMask = maskUpperExp | maskLowerExp;
        overlapMask = maskUpperExp & maskLowerExp;

        %% ----------------------------------------------------------------
        % 4. Register I_U to I_L in the upper-jaw region to get T_{L->U}
        %
        % We use masked images for upper-jaw registration:
        %   T_{L->U} = argmin D_{R_U}( I_U, I_L ∘ T )
        %
        % Here:
        %   fixed  = I_U(upper region)
        %   moving = I_L(upper region)
        %% ----------------------------------------------------------------
        % Register I_U to I_L in the upper-jaw region to get T_{U->L}
        fixedReg  = IL .* maskUpper;
        movingReg = IU .* maskUpper;

        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = 300;
        optimizer.MinimumStepLength = 1e-5;
        optimizer.RelaxationFactor = 0.5;

        fixedRegSm  = imgaussfilt3(fixedReg, 1.0);
        movingRegSm = imgaussfilt3(movingReg, 1.0);

        Rfixed  = imref3d(size(fixedRegSm));
        Rmoving = imref3d(size(movingRegSm));

        tform_U_to_L = imregtform(movingRegSm, Rmoving, fixedRegSm, Rfixed, ...
                                  'rigid', optimizer, metric);

        % Optional: inverse for record only
        Tinv = inv(tform_U_to_L.T);
        Tinv(1:3,4) = 0;
        Tinv(4,4)   = 1;
        tform_L_to_U = affine3d(Tinv);

        if verbose
            fprintf('Rigid registration finished for case %s.\n', caseName);
            fprintf('T_{L->U}:\n');
            disp(tform_L_to_U.T);
            fprintf('T_{U->L} = inv(T_{U->L}):\n');
            disp(tform_U_to_L.T);
        end

        %% ----------------------------------------------------------------
        % 5. Warp I_U into I_L space
        %% ----------------------------------------------------------------
        % Warp I_U into I_L space using T_{U->L}
        IU_hat = imwarp(IU, tform_U_to_L, 'OutputView', imref3d(size(IL)), ...
                'Interp', 'nearest', 'FillValues', 0);

        %% ----------------------------------------------------------------
        % 6. Build w(x) using distance-transform feathering
        %    but only blend inside the overlap region
        %    and use sharper weights to reduce blur
        %% ----------------------------------------------------------------
        upperOnly = maskUpperExp & ~maskLowerExp;
        lowerOnly = maskLowerExp & ~maskUpperExp;
        overlap   = maskUpperExp & maskLowerExp;

        dU = bwdist(~maskUpperExp);
        dL = bwdist(~maskLowerExp);

        % Option 3: make transition steeper
        p = 2;   % try p = 2, 3, or 4
        w = dU.^p ./ (dU.^p + dL.^p + eps);

        %% ----------------------------------------------------------------
        % 7. Compute fused image
        %    Option 2: blend only in overlap
        %% ----------------------------------------------------------------
        % Final image in I_L space
        Ifused = IL;   % default: use I_L everywhere

        % Upper-only region: directly use warped upper image
        Ifused(upperOnly) = IU_hat(upperOnly);

        % Overlap region: weighted blending
        Ifused(overlap) = w(overlap) .* IU_hat(overlap) + ...
                          (1 - w(overlap)) .* IL(overlap);

        %% ----------------------------------------------------------------
        % 8. Save fused image
        %% ----------------------------------------------------------------
        outFile = buildOutputName(caseFolder, upperFile);

        outInfo = infoU;
        if saveAsSingle
            IfusedToSave = single(Ifused);
        else
            IfusedToSave = Ifused;
        end

        % niftiwrite can reuse header-like info from niftiinfo
        niftiwrite(IfusedToSave, outFile, outInfo, 'Compressed', false);

        fprintf('Saved fused image: %s\n', outFile);

        % Save transform matrices as well for traceability
        save(fullfile(caseFolder, 'Di_fused_transform.mat'), ...
             'tform_U_to_L', 'tform_L_to_U', ...
             'idx1', 'idx2', 'idx3', ...
             'hfDim', 'overlapHalfWidth');

    catch ME
        warning('Case %s failed: %s', caseName, ME.message);
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

fprintf('\nAll done.\n');


%% ========================= Local functions ==============================

function filePath = findSingleFile(folder, pattern)
    files = dir(fullfile(folder, pattern));
    if isempty(files)
        filePath = '';
        return;
    end
    if numel(files) > 1
        % pick the first but warn
        warning('Multiple files found for pattern %s in %s. Using the first one.', pattern, folder);
    end
    filePath = fullfile(folder, files(1).name);
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

function mask = fillMaskBetween(mask, dim, a, b)
    lo = min(a, b);
    hi = max(a, b);

    sz = size(mask);
    lo = max(1, min(sz(dim), lo));
    hi = max(1, min(sz(dim), hi));

    idx = repmat({':'}, 1, ndims(mask));
    idx{dim} = lo:hi;
    mask(idx{:}) = true;
end

function out = dilateAlongDim(mask, dim, radius)
    if radius <= 0
        out = mask;
        return;
    end

    kernelSize = ones(1, ndims(mask));
    kernelSize(dim) = 2 * radius + 1;
    kernel = ones(kernelSize);

    out = convn(double(mask), kernel, 'same') > 0;
end

function outFile = buildOutputName(caseFolder, upperFile)
    [~, baseName, ~] = fileparts(upperFile);

    % If input was .nii.gz, fileparts removes only .gz, so handle that
    if endsWith(baseName, '.nii', 'IgnoreCase', true)
        baseName = erase(baseName, '.nii');
    end

    % Replace upperjaw token if present; otherwise append
    if contains(baseName, '_Di_MotCorr_upperjaw_')
        outBase = strrep(baseName, '_Di_MotCorr_upperjaw_', '_Di_fused_nearest_');
    else
        outBase = [baseName, '_Di_fused_nearest_'];
    end

    outFile = fullfile(caseFolder, [outBase, '.nii']);
end
function teethMaskPath = findLatestTeethMask(caseFolder)
    f = dir(fullfile(caseFolder, '*msk_teeth_fullresol.nii.gz'));

    if isempty(f)
        teethMaskPath = '';
        return;
    end

    if numel(f) > 1
        [~, idx] = max([f.datenum]);
        warning('Multiple teeth masks found in %s, using most recent: %s', ...
                caseFolder, f(idx).name);
        teethMaskPath = fullfile(caseFolder, f(idx).name);
    else
        teethMaskPath = fullfile(caseFolder, f(1).name);
    end
end