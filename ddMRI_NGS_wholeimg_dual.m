%% batch_compute_NGS.m
% Batch compute NGS for cases, each with 4 image variants (no mouth_moco).
% Sequence can be selected via seqTag ('T2wSPACE' or 'PDwSPACE').
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2025-12-11

clear; clc;

% -------------------------------------------------------------------------
% Configuration
% -------------------------------------------------------------------------
baseDir = '/home/zn23/matlab/3D-moco-siemens-working/Studies-deploy/ddMRI_new';
% Choose which sequence to process:
%   'T2wSPACE'  or  'PDwSPACE'
seqTag  = 'PDwSPACE';   % <--- change to 'PDwSPACE' for PDw sequence

% Output CSV will include the seqTag in the filename
outCSV  = fullfile(baseDir, ['NGS_results_' seqTag '.csv']);

% Define all cases (folder names are 1–8)
nCases   = 4;
caseList = arrayfun(@num2str, 1:nCases, 'UniformOutput', false);

% Define image types and corresponding filename patterns (no mouth_moco)
imgTypes = {'Aq', 'full_FOV_moco', 'upperjaw_moco', 'lowerjaw_moco'};
patterns = {'_Aq_MotCorr.nii', ...
            'Di_MotCorr.nii', ...
            'Di_MotCorr_upperjaw_.nii', ...
            'Di_MotCorr_lowerjaw_.nii'};

% -------------------------------------------------------------------------
% Prepare result storage
% -------------------------------------------------------------------------
NGS_table = array2table(nan(nCases, numel(imgTypes)), ...
    'VariableNames', imgTypes);
NGS_table.Case = (1:nCases)';

% Move Case column to first position
NGS_table = movevars(NGS_table, 'Case', 'Before', 1);

% -------------------------------------------------------------------------
% Loop over all cases
% -------------------------------------------------------------------------
for iCase = 1:nCases
    caseDir = fullfile(baseDir, caseList{iCase}, 'An-Ve');
    fprintf('Processing Case %d (%s), sequence tag: %s\n', ...
            iCase, caseDir, seqTag);

    for j = 1:numel(imgTypes)
        % Find the file matching both seqTag and pattern
        searchPattern = ['*' seqTag '*' patterns{j}];
        files = dir(fullfile(caseDir, searchPattern));
        if isempty(files)
            warning('Missing file for case %d (%s): pattern %s', ...
                    iCase, seqTag, searchPattern);
            continue;
        end

        imgPath = fullfile(caseDir, files(1).name);

        % Read NIfTI image
        try
            img = niftiread(imgPath);
        catch
            warning('Failed to read image: %s', imgPath);
            continue;
        end

        % Use only half FOV in 2nd dimension (as in your original script)
        img = img(:, 1:ceil(size(img,2)/2), :);
        
        switch seqTag
            case 'T2wSPACE'
                ins_thres = 0.35;
            case 'PDwSPACE'
                ins_thres = 0.15;
            otherwise
                ins_thres = 0.5;
        end

        % Compute NGS using provided function
        try
            NGS_val = computeNGS(img,ins_thres);
        catch ME
            warning('computeNGS failed for case %d (%s, %s): %s', ...
                    iCase, seqTag, imgTypes{j}, ME.message);
            continue;
        end

        % Store in table
        NGS_table{iCase, imgTypes{j}} = NGS_val;
    end
end

% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
writetable(NGS_table, outCSV);
fprintf('\nAll done. Results for %s saved to:\n%s\n', seqTag, outCSV);


% ========================================================================
% Helper functions
% ========================================================================
function ngs = computeNGS(I,ins_thres)
% Normalized Gradient-Squared (3D):
%   ngs = sum(|∇I|^2 over foreground) / sum(I^2 over foreground)
% Foreground mask from adaptive threshold on |I| to suppress background.
% ins_thres is a parameter to adjust the threshold moduling for T
    I = double(I);
    I = I - median(I(:), 'omitnan');       % robust centering
    Iabs = abs(I);

    % Foreground mask (adaptive): Otsu on rescaled |I|
    Imax = max(Iabs(:));
    if Imax <= 0
        ngs = NaN; return;
    end
    Iu = Iabs / Imax;
    try
        T = ins_thres*graythresh(clamp01(Iu));
        mask = Iu > max(T, 0.05);          % floor at 5% to avoid ultra-conservative mask
    catch
        mask = Iu > 0.05;
    end
    mask = imfill3(mask);                  % fill small holes (custom below)

    % Gradients (central differences via gradient)
    [gx, gy, gz] = gradient(I);
    num = sum((gx(mask).^2 + gy(mask).^2 + gz(mask).^2), 'omitnan');
    den = sum((I(mask)).^2, 'omitnan') + eps;

    ngs = num / den;
end

function B = clamp01(A)
    B = min(max(A,0),1);
end

function BWf = imfill3(BW)
% Simple 3D morphological hole-filling & cleanup for logical masks
    if ~islogical(BW), BW = BW > 0; else, BW = BW; end
    % Remove small specks
    BW = bwareaopen(BW, 100);
    % 3D closing to connect small gaps
    se = strel('sphere', 1);
    BW = imclose(BW, se);
    % Fill slicewise to be robust if 3D fill not available
    BWf = false(size(BW));
    for k = 1:size(BW,3)
        BWf(:,:,k) = imfill(BW(:,:,k), 'holes');
    end
end
