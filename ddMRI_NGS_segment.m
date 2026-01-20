%% batch_compute_NGS_segmented.m
% Batch compute NGS for 8 cases, each with 5 image types and 3 anatomical segments.
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 20-Oct-2025

clear; clc;

% -------------------------------------------------------------------------
% Configuration
% -------------------------------------------------------------------------
baseDir = '/data/gadgetron/matlab_study/ddMRI';
outCSV  = fullfile(baseDir, 'NGS_results_segmented.csv');
nCases  = 8;

imgTypes = {'Aq', 'full_FOV_moco', 'mouth_moco', 'upperjaw_moco', 'lowerjaw_moco'};
patterns = {'_Aq_MotCorr.nii', ...
             'Di_MotCorr.nii', ...
             'Di_MotCorr_mouth_.nii', ...
             'Di_MotCorr_upperjaw_.nii', ...
             'Di_MotCorr_lowerjaw_.nii'};

segNames = {'mouth', 'upperjaw', 'lowerjaw'};

% -------------------------------------------------------------------------
% Prepare output table
% Columns: Case | Aq_mouth | full_FOV_moco_mouth | ... | Aq_upperjaw | ... | Aq_lowerjaw | ...
% -------------------------------------------------------------------------
colNames = {};
for s = 1:numel(segNames)
    for t = 1:numel(imgTypes)
        colNames{end+1} = sprintf('%s_%s', imgTypes{t}, segNames{s}); %#ok<SAGROW>
    end
end
NGS_table = array2table(nan(nCases, numel(colNames)), 'VariableNames', colNames);
NGS_table.Case = (1:nCases)';
NGS_table = movevars(NGS_table, 'Case', 'Before', 1);

% -------------------------------------------------------------------------
% Loop through all cases
% -------------------------------------------------------------------------
for iCase = 1:nCases
    caseDir = fullfile(baseDir, num2str(iCase));
    anveDir = fullfile(caseDir, 'An-Ve');
    fprintf('Processing Case %d: %s\n', iCase, caseDir);

    %% ---- Read location.txt ----
    locFile = fullfile(caseDir, 'location.txt');
    if ~isfile(locFile)
        warning('location.txt missing for case %d', iCase);
        continue;
    end
    loc = readmatrix(locFile);
    if numel(loc) < 3
        warning('Invalid location.txt in case %d (expect 3 values)', iCase);
        continue;
    end
    loc = sort(loc, 'descend');  % ensure idx(1)>idx(2)>idx(3)
    idx1 = loc(1); idx2 = loc(2); idx3 = loc(3);

    % Define segment ranges along first dimension
    segRanges.mouth     = idx3:idx1;
    segRanges.upperjaw  = idx2:idx1;
    segRanges.lowerjaw  = idx3:idx2;

    %% ---- Loop through all image types ----
    for j = 1:numel(imgTypes)
        files = dir(fullfile(anveDir, ['*' patterns{j}]));
        if isempty(files)
            warning('Missing file for case %d: pattern %s', iCase, patterns{j});
            continue;
        end

        imgPath = fullfile(anveDir, files(1).name);

        try
            img = niftiread(imgPath);
        catch ME
            warning('Failed to read %s: %s', imgPath, ME.message);
            continue;
        end

        % ---- Compute NGS for each segment ----
        for s = 1:numel(segNames)
            seg = segNames{s};
            range1 = segRanges.(seg);

            if max(range1) > size(img,1)
                warning('Segment range exceeds image size for case %d (%s)', iCase, seg);
                continue;
            end

            img_seg = img(range1,1:ceil(size(img,2)/2),:);

            try
                NGS_val = computeNGS(img_seg);
            catch ME
                warning('computeNGS failed (case %d, %s, %s): %s', iCase, imgTypes{j}, seg, ME.message);
                continue;
            end

            colName = sprintf('%s_%s', imgTypes{j}, seg);
            NGS_table{iCase, colName} = NGS_val;
        end
    end
end

% -------------------------------------------------------------------------
% Save results
% -------------------------------------------------------------------------
writetable(NGS_table, outCSV);
fprintf('\n✅ All done. Results saved to:\n%s\n', outCSV);


function ngs = computeNGS(I)
% Normalized Gradient-Squared (3D):
%   ngs = sum(|∇I|^2 over foreground) / sum(I^2 over foreground)
% Foreground mask from adaptive threshold on |I| to suppress background.
    I = double(I);
    I = I - median(I(:), 'omitnan');       % robust centering
    Iabs = abs(I);

    % Foreground mask (adaptive): Otsu on rescaled |I|
    Imax = max(Iabs(:));
    if Imax <= 0
        ngs = NaN; return;
    end
    Iu = Iabs / Imax;
%     try
%         T = graythresh(clamp01(Iu));
%         mask = Iu > max(T, 0.05);          % floor at 5% to avoid ultra-conservative mask
%     catch
%         mask = Iu > 0.05;
%     end
    mask = Iu > mean(Iu(:))*0.5;
    mask = bwareaopen(mask, 700);  
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




