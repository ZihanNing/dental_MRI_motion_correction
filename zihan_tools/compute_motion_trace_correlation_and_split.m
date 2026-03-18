function compute_motion_trace_correlation_and_split()
% Compute pairwise motion-trace correlations among three reconstructions
% and split cases into motion-pattern groups using a more general rule:
% regional motion is flagged when at least one regional trace differs from
% the MoCo-fullFOV trace, and the upper-vs-lower difference is even larger.
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

% Thresholds for classification
% diff = 1 - correlation

% Mean-based rule
thr_diff_from_full = 0.46;
thr_margin_ul      = 0.08;

% Dimension-wise rule (more sensitive)
thr_diff_from_full_dim = 0.40;
thr_margin_ul_dim      = 0.05;
minFlaggedDims         = 2;    % require at least 2 dimensions to support regional motion
% Dimension names for the 6 motion parameters
dimNames = {'Tra1','Tra2','Tra3','Rot1','Rot2','Rot3'};

% Output files
outMat   = fullfile(rootDir, 'motion_trace_correlation_results_v2.mat');
outXlsx  = fullfile(rootDir, 'motion_trace_correlation_results_v2.xlsx');

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
[caseNums, sortIdx] = sort(caseNums);
caseNames = caseNames(sortIdx);

nCases = numel(caseNames);
fprintf('Found %d case folders.\n', nCases);

%% =========================
%  Prepare result containers
%  =========================
results = struct();

CaseID = strings(nCases,1);

AvgCorr_Tra1 = nan(nCases,1);
AvgCorr_Tra2 = nan(nCases,1);
AvgCorr_Tra3 = nan(nCases,1);
AvgCorr_Rot1 = nan(nCases,1);
AvgCorr_Rot2 = nan(nCases,1);
AvgCorr_Rot3 = nan(nCases,1);

Mean_rFU = nan(nCases,1);
Mean_rFL = nan(nCases,1);
Mean_rUL = nan(nCases,1);

Diff_FU = nan(nCases,1);
Diff_FL = nan(nCases,1);
Diff_UL = nan(nCases,1);

MaxDiffFromFull = nan(nCases,1);
ULminusMaxFullDiff = nan(nCases,1);

Worst_ULadv = nan(nCases,1);
SecondWorst_ULadv = nan(nCases,1);
Worst_Dim_Name = strings(nCases,1);
SecondWorst_Dim_Name = strings(nCases,1);
NumFlaggedDims = nan(nCases,1);
FlaggedDimNames = strings(nCases,1);

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

        % ---- Compute per-dimension correlations ----
        rFU = nan(1,6);
        rFL = nan(1,6);
        rUL = nan(1,6);
        avgCorr3 = nan(1,6);
        
        diffFU_dim = nan(1,6);
        diffFL_dim = nan(1,6);
        diffUL_dim = nan(1,6);
        maxDiffFromFull_dim = nan(1,6);
        ULadv_dim = nan(1,6);
        flagDim = false(1,6);

        for idim = 1:6
            xF = T_full(:, idim);
            xU = T_upper(:, idim);
            xL = T_lower(:, idim);

            rFU(idim) = safe_corr(xF, xU);
            rFL(idim) = safe_corr(xF, xL);
            rUL(idim) = safe_corr(xU, xL);

            avgCorr3(idim) = fisher_mean_corr([rFU(idim), rFL(idim), rUL(idim)]);

            % Dimension-wise differences
            diffFU_dim(idim) = 1 - rFU(idim);
            diffFL_dim(idim) = 1 - rFL(idim);
            diffUL_dim(idim) = 1 - rUL(idim);

            maxDiffFromFull_dim(idim) = max(diffFU_dim(idim), diffFL_dim(idim));
            ULadv_dim(idim) = diffUL_dim(idim) - maxDiffFromFull_dim(idim);

            flagDim(idim) = (maxDiffFromFull_dim(idim) >= thr_diff_from_full_dim) && ...
                            (ULadv_dim(idim) >= thr_margin_ul_dim);
        end

        % ---- Mean correlations across 6 dimensions ----
        mean_rFU = mean(rFU, 'omitnan');
        mean_rFL = mean(rFL, 'omitnan');
        mean_rUL = mean(rUL, 'omitnan');

        % ---- Convert correlation to difference ----
        diffFU = 1 - mean_rFU;
        diffFL = 1 - mean_rFL;
        diffUL = 1 - mean_rUL;

        maxDiffFromFull = max(diffFU, diffFL);
        ulMinusMaxFull  = diffUL - maxDiffFromFull;
        
        % ---- Dimension-wise ranking ----
        [sortedULadv, sortIdxULadv] = sort(ULadv_dim, 'descend', 'MissingPlacement', 'last');

        worst_ULadv = sortedULadv(1);
        worst_dim_name = string(dimNames{sortIdxULadv(1)});

        if numel(sortedULadv) >= 2
            secondWorst_ULadv = sortedULadv(2);
            secondWorst_dim_name = string(dimNames{sortIdxULadv(2)});
        else
            secondWorst_ULadv = NaN;
            secondWorst_dim_name = "";
        end

        nFlaggedDims = sum(flagDim);

        if nFlaggedDims > 0
            flaggedNames = string(dimNames(flagDim));
            flaggedDimNamesStr = strjoin(flaggedNames, ',');
        else
            flaggedDimNamesStr = "";
        end

        % ---- Updated classification rule ----
        % Mean-based rule
        passMeanRule = (maxDiffFromFull >= thr_diff_from_full) && ...
                       (diffUL >= maxDiffFromFull + thr_margin_ul);

        % Dimension-wise rule
        passDimRule = (nFlaggedDims >= minFlaggedDims);

        if passMeanRule || passDimRule
            patternLabel = "IndependentRegionalMotion";
        else
            patternLabel = "MostlySharedMotion";
        end

        % ---- Save into table variables ----
        AvgCorr_Tra1(icase) = avgCorr3(1);
        AvgCorr_Tra2(icase) = avgCorr3(2);
        AvgCorr_Tra3(icase) = avgCorr3(3);
        AvgCorr_Rot1(icase) = avgCorr3(4);
        AvgCorr_Rot2(icase) = avgCorr3(5);
        AvgCorr_Rot3(icase) = avgCorr3(6);

        Mean_rFU(icase) = mean_rFU;
        Mean_rFL(icase) = mean_rFL;
        Mean_rUL(icase) = mean_rUL;

        Diff_FU(icase) = diffFU;
        Diff_FL(icase) = diffFL;
        Diff_UL(icase) = diffUL;

        MaxDiffFromFull(icase) = maxDiffFromFull;
        ULminusMaxFullDiff(icase) = ulMinusMaxFull;
        
        Worst_ULadv(icase) = worst_ULadv;
        SecondWorst_ULadv(icase) = secondWorst_ULadv;
        Worst_Dim_Name(icase) = worst_dim_name;
        SecondWorst_Dim_Name(icase) = secondWorst_dim_name;
        NumFlaggedDims(icase) = nFlaggedDims;
        FlaggedDimNames(icase) = flaggedDimNamesStr;

        MotionPattern(icase) = patternLabel;
        Status(icase) = "OK";

        % ---- Save detailed results ----
        results(icase).caseID = caseID;
        results(icase).fileFull = fileFull;
        results(icase).fileUpper = fileUpper;
        results(icase).fileLower = fileLower;
        results(icase).nShotUsed = nShot;

        results(icase).rFU = rFU;
        results(icase).rFL = rFL;
        results(icase).rUL = rUL;
        results(icase).avgCorr3 = avgCorr3;

        results(icase).mean_rFU = mean_rFU;
        results(icase).mean_rFL = mean_rFL;
        results(icase).mean_rUL = mean_rUL;

        results(icase).diffFU = diffFU;
        results(icase).diffFL = diffFL;
        results(icase).diffUL = diffUL;
        results(icase).maxDiffFromFull = maxDiffFromFull;
        results(icase).ulMinusMaxFull = ulMinusMaxFull;

        results(icase).motionPattern = patternLabel;
        
        results(icase).diffFU_dim = diffFU_dim;
        results(icase).diffFL_dim = diffFL_dim;
        results(icase).diffUL_dim = diffUL_dim;
        results(icase).maxDiffFromFull_dim = maxDiffFromFull_dim;
        results(icase).ULadv_dim = ULadv_dim;
        results(icase).flagDim = flagDim;
        results(icase).nFlaggedDims = nFlaggedDims;
        results(icase).flaggedDimNames = flaggedDimNamesStr;
        results(icase).worst_ULadv = worst_ULadv;
        results(icase).secondWorst_ULadv = secondWorst_ULadv;
        results(icase).worst_dim_name = worst_dim_name;
        results(icase).secondWorst_dim_name = secondWorst_dim_name;
        results(icase).passMeanRule = passMeanRule;
        results(icase).passDimRule = passDimRule;

        fprintf('Case %s done.\n', caseID);
        fprintf('  mean r(F,U) = %.4f\n', mean_rFU);
        fprintf('  mean r(F,L) = %.4f\n', mean_rFL);
        fprintf('  mean r(U,L) = %.4f\n', mean_rUL);
        fprintf('  diff(F,U)   = %.4f\n', diffFU);
        fprintf('  diff(F,L)   = %.4f\n', diffFL);
        fprintf('  diff(U,L)   = %.4f\n', diffUL);
        fprintf('  maxDiffFromFull   = %.4f\n', maxDiffFromFull);
        fprintf('  ULminusMaxFull    = %.4f\n', ulMinusMaxFull);
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
    AvgCorr_Tra1, AvgCorr_Tra2, AvgCorr_Tra3, ...
    AvgCorr_Rot1, AvgCorr_Rot2, AvgCorr_Rot3, ...
    Mean_rFU, Mean_rFL, Mean_rUL, ...
    Diff_FU, Diff_FL, Diff_UL, ...
    MaxDiffFromFull, ULminusMaxFullDiff, ...
    Worst_ULadv, SecondWorst_ULadv, ...
    Worst_Dim_Name, SecondWorst_Dim_Name, ...
    NumFlaggedDims, FlaggedDimNames, ...
    MotionPattern, Status);

%% =========================
%  Save outputs
%  =========================
save(outMat, 'results', 'T_summary', ...
    'thr_diff_from_full', 'thr_margin_ul', 'dimNames');

writetable(T_summary, outXlsx, 'Sheet', 'Summary');

% Detailed sheet
detailedCell = cell(nCases, 1 + 18);
header = [{'CaseID'}, ...
          strcat('rFU_', dimNames), ...
          strcat('rFL_', dimNames), ...
          strcat('rUL_', dimNames)];

for i = 1:nCases
    detailedCell{i,1} = char(CaseID(i));
    if isfield(results(i), 'rFU')
        detailedCell(i, 2:7)   = num2cell(results(i).rFU);
        detailedCell(i, 8:13)  = num2cell(results(i).rFL);
        detailedCell(i, 14:19) = num2cell(results(i).rUL);
    else
        detailedCell(i, 2:19) = {nan};
    end
end

T_detail = cell2table(detailedCell, 'VariableNames', header);
writetable(T_detail, outXlsx, 'Sheet', 'PairwiseCorr');

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
function rMean = fisher_mean_corr(rVals)
rVals = rVals(isfinite(rVals));

if isempty(rVals)
    rMean = nan;
    return;
end

rVals = max(min(rVals, 0.999999), -0.999999);

z = atanh(rVals);
rMean = tanh(mean(z));
end