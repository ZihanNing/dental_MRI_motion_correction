function compute_motion_trace_correlation_and_split_outlierremoval()
% Compare MoCo-upper and MoCo-lower motion traces to identify
% lower-jaw-involved motion.
%
% Logic:
% 1) load T from MoCo-fullFOV / MoCo-upper / MoCo-lower
% 2) align shot number
% 3) detect shared outlier time points across reconstructions within a case
% 4) replace outlier rows by nearest valid row
% 5) compute Pearson correlation only between MoCo-upper and MoCo-lower
% 6) convert to difference d = 1-r
% 7) use only WorstDiffUL for classification:
%    if WorstDiffUL is large, classify as LowerJawInvolvedMotion
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 17-Mar-2026


clc;
clear;
close all;

%% =========================
%  User settings
%  =========================
rootDir = '/home/zn23/Data/ddMRI/';

% Threshold for WorstDiffUL
thr_worst_diff_ul = 0.9;   % equivalent to corr <= 0.30 when d = 1-r

% Outlier detection threshold (robust z-score on shot-to-shot change)
outlier_z_thr = 15.0;

% Dimension names
dimNames = {'Tra1','Tra2','Tra3','Rot1','Rot2','Rot3'};

% Output files
outMat   = fullfile(rootDir, 'motion_trace_ULsplit_results.mat');
outXlsx  = fullfile(rootDir, 'motion_trace_ULsplit_results.xlsx');

%% =========================
%  Find case folders
%  =========================
d = dir(rootDir);
isSub = [d.isdir];
subNames = {d(isSub).name};

isCase = false(size(subNames));
for i = 1:numel(subNames)
    isCase(i) = ~isempty(regexp(subNames{i}, '^\d+$', 'once'));
end
caseNames = subNames(isCase);

caseNums = cellfun(@str2double, caseNames);
[~, sortIdx] = sort(caseNums);
caseNames = caseNames(sortIdx);

nCases = numel(caseNames);
fprintf('Found %d case folders.\n', nCases);

%% =========================
%  Prepare result containers
%  =========================
results = struct();

CaseID = strings(nCases,1);

Corr_UL_Tra1 = nan(nCases,1);
Corr_UL_Tra2 = nan(nCases,1);
Corr_UL_Tra3 = nan(nCases,1);
Corr_UL_Rot1 = nan(nCases,1);
Corr_UL_Rot2 = nan(nCases,1);
Corr_UL_Rot3 = nan(nCases,1);

Diff_UL_Tra1 = nan(nCases,1);
Diff_UL_Tra2 = nan(nCases,1);
Diff_UL_Tra3 = nan(nCases,1);
Diff_UL_Rot1 = nan(nCases,1);
Diff_UL_Rot2 = nan(nCases,1);
Diff_UL_Rot3 = nan(nCases,1);

Mean_rUL = nan(nCases,1);
Mean_diffUL = nan(nCases,1);

WorstDiffUL = nan(nCases,1);
SecondWorstDiffUL = nan(nCases,1);
Worst_Dim_Name = strings(nCases,1);
SecondWorst_Dim_Name = strings(nCases,1);

NumFlaggedDims = nan(nCases,1);
FlaggedDimNames = strings(nCases,1);

NumOutlierShots = nan(nCases,1);
OutlierShotIdx = strings(nCases,1);

MotionPattern = strings(nCases,1);
Status = strings(nCases,1);

%% =========================
%  Main loop
%  =========================
for icase = 1:nCases
    caseID = caseNames{icase};
    caseFolder = fullfile(rootDir, caseID);
    CaseID(icase) = string(caseID);

    fprintf('\n==============================\n');
    fprintf('Processing case %s\n', caseID);

    % ---- Find files ----
    fileFull  = find_single_file(caseFolder, '*_Tr_MotCorr.mat');
    fileUpper = find_single_file(caseFolder, '*_Tr_MotCorr_upperjaw_*.mat');
    fileLower = find_single_file(caseFolder, '*_Tr_MotCorr_lowerjaw_*.mat');

    if isempty(fileFull) || isempty(fileUpper) || isempty(fileLower)
        warning('Case %s: missing one or more required files. Skipping.', caseID);
        Status(icase) = "Missing file(s)";
        MotionPattern(icase) = "Unknown";
        continue;
    end

    try
        % ---- Load motion traces ----
        T_full  = load_motion_trace(fileFull);
        T_upper = load_motion_trace(fileUpper);
        T_lower = load_motion_trace(fileLower);

        if size(T_full,2) ~= 6 || size(T_upper,2) ~= 6 || size(T_lower,2) ~= 6
            warning('Case %s: one or more traces are not [N x 6]. Skipping.', caseID);
            Status(icase) = "Invalid trace shape";
            MotionPattern(icase) = "Unknown";
            continue;
        end

        % ---- Align shot numbers ----
        nShot = min([size(T_full,1), size(T_upper,1), size(T_lower,1)]);
        if nShot < 3
            warning('Case %s: too few shots (%d). Skipping.', caseID, nShot);
            Status(icase) = "Too few shots";
            MotionPattern(icase) = "Unknown";
            continue;
        end

        T_full  = T_full(1:nShot, :);
        T_upper = T_upper(1:nShot, :);
        T_lower = T_lower(1:nShot, :);

        % ---- Detect shared outlier shot indices ----
        outlierIdx = detect_shared_outlier_rows(T_full, T_upper, T_lower, outlier_z_thr);

        % ---- Replace outlier rows with nearest valid row ----
        T_full_clean  = replace_outlier_rows_with_nearest(T_full, outlierIdx);
        T_upper_clean = replace_outlier_rows_with_nearest(T_upper, outlierIdx);
        T_lower_clean = replace_outlier_rows_with_nearest(T_lower, outlierIdx);

        % ---- Compute UL correlations only ----
        rUL = nan(1,6);
        diffUL_dim = nan(1,6);
        flagDim = false(1,6);

        for idim = 1:6
            xU = T_upper_clean(:, idim);
            xL = T_lower_clean(:, idim);

            rUL(idim) = safe_corr(xU, xL);
            diffUL_dim(idim) = 1 - rUL(idim);

            flagDim(idim) = (diffUL_dim(idim) >= thr_worst_diff_ul);
        end

        mean_rUL = mean(rUL, 'omitnan');
        mean_diffUL = mean(diffUL_dim, 'omitnan');

        % ---- Dimension ranking ----
        [sortedDiffUL, sortIdx] = sort(diffUL_dim, 'descend', 'MissingPlacement', 'last');

        worstDiffUL = sortedDiffUL(1);
        worst_dim_name = string(dimNames{sortIdx(1)});

        if numel(sortedDiffUL) >= 2
            secondWorstDiffUL = sortedDiffUL(2);
            secondWorst_dim_name = string(dimNames{sortIdx(2)});
        else
            secondWorstDiffUL = NaN;
            secondWorst_dim_name = "";
        end

        nFlaggedDims = sum(flagDim);

        if nFlaggedDims > 0
            flaggedNames = string(dimNames(flagDim));
            flaggedDimNamesStr = strjoin(flaggedNames, ',');
        else
            flaggedDimNamesStr = "";
        end

        if any(outlierIdx)
            outlierIdxStr = strjoin(string(find(outlierIdx)), ',');
        else
            outlierIdxStr = "";
        end

        % ---- Classification: use only WorstDiffUL ----
        if isfinite(worstDiffUL) && (worstDiffUL >= thr_worst_diff_ul)
            patternLabel = "LowerJawInvolvedMotion";
        else
            patternLabel = "MostlySharedMotion";
        end

        % ---- Save into table variables ----
        Corr_UL_Tra1(icase) = rUL(1);
        Corr_UL_Tra2(icase) = rUL(2);
        Corr_UL_Tra3(icase) = rUL(3);
        Corr_UL_Rot1(icase) = rUL(4);
        Corr_UL_Rot2(icase) = rUL(5);
        Corr_UL_Rot3(icase) = rUL(6);

        Diff_UL_Tra1(icase) = diffUL_dim(1);
        Diff_UL_Tra2(icase) = diffUL_dim(2);
        Diff_UL_Tra3(icase) = diffUL_dim(3);
        Diff_UL_Rot1(icase) = diffUL_dim(4);
        Diff_UL_Rot2(icase) = diffUL_dim(5);
        Diff_UL_Rot3(icase) = diffUL_dim(6);

        Mean_rUL(icase) = mean_rUL;
        Mean_diffUL(icase) = mean_diffUL;

        WorstDiffUL(icase) = worstDiffUL;
        SecondWorstDiffUL(icase) = secondWorstDiffUL;
        Worst_Dim_Name(icase) = worst_dim_name;
        SecondWorst_Dim_Name(icase) = secondWorst_dim_name;

        NumFlaggedDims(icase) = nFlaggedDims;
        FlaggedDimNames(icase) = flaggedDimNamesStr;

        NumOutlierShots(icase) = sum(outlierIdx);
        OutlierShotIdx(icase) = outlierIdxStr;

        MotionPattern(icase) = patternLabel;
        Status(icase) = "OK";

        % ---- Save detailed results ----
        results(icase).caseID = caseID;
        results(icase).fileFull = fileFull;
        results(icase).fileUpper = fileUpper;
        results(icase).fileLower = fileLower;
        results(icase).nShotUsed = nShot;

        results(icase).T_full_raw = T_full;
        results(icase).T_upper_raw = T_upper;
        results(icase).T_lower_raw = T_lower;

        results(icase).outlierIdx = outlierIdx;
        results(icase).T_full_clean = T_full_clean;
        results(icase).T_upper_clean = T_upper_clean;
        results(icase).T_lower_clean = T_lower_clean;

        results(icase).rUL = rUL;
        results(icase).diffUL_dim = diffUL_dim;
        results(icase).flagDim = flagDim;
        results(icase).mean_rUL = mean_rUL;
        results(icase).mean_diffUL = mean_diffUL;
        results(icase).worstDiffUL = worstDiffUL;
        results(icase).secondWorstDiffUL = secondWorstDiffUL;
        results(icase).worst_dim_name = worst_dim_name;
        results(icase).secondWorst_dim_name = secondWorst_dim_name;
        results(icase).nFlaggedDims = nFlaggedDims;
        results(icase).flaggedDimNames = flaggedDimNamesStr;
        results(icase).motionPattern = patternLabel;

        fprintf('Case %s done.\n', caseID);
        fprintf('  mean r(U,L)        = %.4f\n', mean_rUL);
        fprintf('  mean diff(U,L)     = %.4f\n', mean_diffUL);
        fprintf('  worst diff(U,L)    = %.4f (%s)\n', worstDiffUL, worst_dim_name);
        fprintf('  2nd worst diff(U,L)= %.4f (%s)\n', secondWorstDiffUL, secondWorst_dim_name);
        fprintf('  flagged dims       = %d [%s]\n', nFlaggedDims, flaggedDimNamesStr);
        fprintf('  outlier shots      = %d [%s]\n', sum(outlierIdx), outlierIdxStr);
        fprintf('  Classified as: %s\n', patternLabel);

    catch ME
        warning('Case %s failed: %s', caseID, ME.message);
        Status(icase) = "Failed";
        MotionPattern(icase) = "Unknown";

        results(icase).caseID = caseID;
        results(icase).error = ME.message;
    end
end

%% =========================
%  Build summary table
%  =========================
T_summary = table( ...
    CaseID, ...
    Corr_UL_Tra1, Corr_UL_Tra2, Corr_UL_Tra3, ...
    Corr_UL_Rot1, Corr_UL_Rot2, Corr_UL_Rot3, ...
    Diff_UL_Tra1, Diff_UL_Tra2, Diff_UL_Tra3, ...
    Diff_UL_Rot1, Diff_UL_Rot2, Diff_UL_Rot3, ...
    Mean_rUL, Mean_diffUL, ...
    WorstDiffUL, SecondWorstDiffUL, ...
    Worst_Dim_Name, SecondWorst_Dim_Name, ...
    NumFlaggedDims, FlaggedDimNames, ...
    NumOutlierShots, OutlierShotIdx, ...
    MotionPattern, Status);

%% =========================
%  Save outputs
%  =========================
save(outMat, 'results', 'T_summary', ...
    'thr_worst_diff_ul', 'outlier_z_thr', 'dimNames');

writetable(T_summary, outXlsx, 'Sheet', 'Summary');

% Detailed sheet
detailedCell = cell(nCases, 1 + 12);
header = [{'CaseID'}, ...
          strcat('rUL_', dimNames), ...
          strcat('diffUL_', dimNames)];

for i = 1:nCases
    detailedCell{i,1} = char(CaseID(i));
    if isfield(results(i), 'rUL')
        detailedCell(i, 2:7)   = num2cell(results(i).rUL);
        detailedCell(i, 8:13)  = num2cell(results(i).diffUL_dim);
    else
        detailedCell(i, 2:13) = {nan};
    end
end

T_detail = cell2table(detailedCell, 'VariableNames', header);
writetable(T_detail, outXlsx, 'Sheet', 'UL_Corr');

fprintf('\n========================================\n');
fprintf('Finished.\n');
fprintf('MAT file : %s\n', outMat);
fprintf('Excel    : %s\n', outXlsx);
fprintf('========================================\n');

end

%% =========================================================
function filePath = find_single_file(folderPath, pattern)
dd = dir(fullfile(folderPath, pattern));

if isempty(dd)
    filePath = '';
    return;
end

[~, idx] = sort({dd.name});
dd = dd(idx);

filePath = fullfile(folderPath, dd(1).name);

if numel(dd) > 1
    warning('Multiple files found for pattern "%s" in %s. Using: %s', ...
        pattern, folderPath, dd(1).name);
end
end

%% =========================================================
function Ttrace = load_motion_trace(matFile)
S = load(matFile);

if ~isfield(S, 'T')
    error('File %s does not contain variable T.', matFile);
end

Ttrace = squeeze(S.T);

if ~ismatrix(Ttrace)
    error('After squeeze(T), result is not a 2D matrix in file: %s', matFile);
end

if size(Ttrace,2) ~= 6 && size(Ttrace,1) == 6
    Ttrace = Ttrace.';
end
end

%% =========================================================
function r = safe_corr(x, y)
x = x(:);
y = y(:);

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 3
    r = nan;
    return;
end

if std(x) < eps || std(y) < eps
    r = 0;
    return;
end

C = corrcoef(x, y);
r = C(1,2);

r = max(min(r, 0.999999), -0.999999);
end

%% =========================================================
function outlierIdx = detect_shared_outlier_rows(T_full, T_upper, T_lower, zThr)
% Detect shot indices with unusually large shot-to-shot change
% that appear across all reconstructions within the case.
%
% Strategy:
% 1) compute row-wise motion-change magnitude between consecutive shots
% 2) use robust z-score based on MAD
% 3) if either side of a shot is extreme in all three reconstructions,
%    mark that shot as an outlier row

nShot = size(T_full,1);

scoreF = row_change_score(T_full);
scoreU = row_change_score(T_upper);
scoreL = row_change_score(T_lower);

badF = scoreF > zThr;
badU = scoreU > zThr;
badL = scoreL > zThr;

outlierIdx = badF & badU & badL;

% avoid wiping out too many rows
if sum(outlierIdx) >= nShot
    outlierIdx(:) = false;
end
end

%% =========================================================
function score = row_change_score(T)
% For each row, quantify how abnormal the shot-to-shot change is.
% Score is based on the sum of absolute changes across 6 dimensions,
% then converted to robust z-score.

nShot = size(T,1);

if nShot < 2
    score = zeros(nShot,1);
    return;
end

dT = diff(T, 1, 1);                    % (nShot-1) x 6
changeMag = sum(abs(dT), 2);           % (nShot-1) x 1

% robust z-score on change magnitude
medv = median(changeMag, 'omitnan');
madv = mad(changeMag, 1);

if madv < eps
    z = zeros(size(changeMag));
else
    z = abs(changeMag - medv) ./ (1.4826 * madv);
end

% assign edge-change z-scores back to rows
score = zeros(nShot,1);
score(1) = z(1);
score(end) = z(end);
for i = 2:nShot-1
    score(i) = max(z(i-1), z(i));
end
end

%% =========================================================
function Tclean = replace_outlier_rows_with_nearest(T, outlierIdx)
% Replace each outlier row by the nearest non-outlier row.

Tclean = T;
nShot = size(T,1);

validIdx = find(~outlierIdx);

if isempty(validIdx)
    return;
end

badIdx = find(outlierIdx);

for i = 1:numel(badIdx)
    k = badIdx(i);

    [~, j] = min(abs(validIdx - k));
    nearestIdx = validIdx(j);

    Tclean(k,:) = T(nearestIdx,:);
end
end