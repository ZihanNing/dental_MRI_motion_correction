%% batch_compute_NGS_segmented.m
% Batch compute NGS for 8 cases, each with 4 image types (no mouth_moco)
% and 3 anatomical segments, for a selected sequence type.
%
% Sequence is selected via seqTag ('T2wSPACE' or 'PDwSPACE').
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 11-Dec-2025

clear; clc;

% -------------------------------------------------------------------------
% Configuration
% -------------------------------------------------------------------------
baseDir = '/data/gadgetron/matlab_study/ddMRI_new';
nCases  = 8;

% Choose which sequence to process: 'T2wSPACE' or 'PDwSPACE'
seqTag  = 'T2wSPACE';   % <--- change to 'PDwSPACE' for PDw sequence

% Output CSV will include seqTag in its name
outCSV  = fullfile(baseDir, ['NGS_results_segmented_' seqTag '.csv']);

% Image types (mouth_moco removed) and patterns
imgTypes = {'Aq', 'full_FOV_moco', 'upperjaw_moco', 'lowerjaw_moco'};
patterns = {'_Aq_MotCorr.nii', ...
            'Di_MotCorr.nii', ...
            'Di_MotCorr_upperjaw_.nii', ...
            'Di_MotCorr_lowerjaw_.nii'};

% Segments (unchanged)
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
    fprintf('Processing Case %d (%s), sequence tag: %s\n', iCase, caseDir, seqTag);

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
        % Only pick files that match both seqTag and the pattern
        searchPattern = ['*' seqTag '*' patterns{j}];
        files = dir(fullfile(anveDir, searchPattern));
        if isempty(files)
            warning('Missing file for case %d (%s): pattern %s', ...
                    iCase, seqTag, searchPattern);
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

            img_seg = img(range1, 1:ceil(size(img,2)/2), :);

            try
                NGS_val = computeNGS(img_seg);
            catch ME
                warning('computeNGS failed (case %d, %s, %s, %s): %s', ...
                        iCase, seqTag, imgTypes{j}, seg, ME.message);
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
fprintf('\nAll done. Segmented results for %s saved to:\n%s\n', seqTag, outCSV);


% ========================================================================
% Helper functions
% ========================================================================
function ngs = computeNGS(I)
% Normalized Gradient-Squared (3D):
%   ngs = sum(|∇I|^2 over foreground) / sum(I^2 over foreground)
% Foreground mask from adaptive threshold on |I| to suppress background.
    I = double(I);
    I = I - median(I(:), 'omitnan');       % robust centering
    Iabs = abs(I);

    % Foreground mask (simple heuristic based on mean intensity)
    Imax = max(Iabs(:));
    if Imax <= 0
        ngs = NaN; return;
    end
    Iu = Iabs / Imax;
    mask = Iu > mean(Iu(:)) * 0.5;
    mask = bwareaopen(mask, 700);
    mask = imfill3(mask);                  % fill small holes (custom below)
    mask(:,size(mask,2)-round(size(mask,2)*1/4):end,:)=0;

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
