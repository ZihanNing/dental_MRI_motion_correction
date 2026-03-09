%% conclude_all_NGS_in_table.m
% Collect NGS_table from all cases under a root folder and export to MAT + Excel.
% Computation based on the head mask
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-03

clear; clc;

rootDir = '/home/zn23/Data/ddMRI/';

outMat  = fullfile(rootDir, 'NGS_all_cases.mat');
outXlsx = fullfile(rootDir, 'NGS_all_cases.xlsx');

% -------------------------------------------------------------------------
% 1) List case folders
% -------------------------------------------------------------------------
d = dir(rootDir);
isCase = [d.isdir] & ~startsWith({d.name}, '.');
caseFolders = d(isCase);

if isempty(caseFolders)
    error('No case folders found under: %s', rootDir);
end

% Sort by numeric value if possible
caseNames = {caseFolders.name};
caseNums = nan(size(caseNames));
for i = 1:numel(caseNames)
    tmp = str2double(caseNames{i});
    if ~isnan(tmp), caseNums(i) = tmp; end
end
if all(~isnan(caseNums))
    [~, idxSort] = sort(caseNums);
else
    [~, idxSort] = sort(lower(caseNames));
end
caseFolders = caseFolders(idxSort);

% -------------------------------------------------------------------------
% 2) Loop and collect each case as a struct (stored in a cell array)
% -------------------------------------------------------------------------
allRows = {};      % cell array of structs, one per case
allVarNames = {};  % union of flattened variable names

for c = 1:numel(caseFolders)
    caseID = caseFolders(c).name;
    caseDir = fullfile(rootDir, caseID);

    % Find NGS mat file (recursive)
    f = dir(fullfile(caseDir, '**', '*_NGS.mat'));
    if isempty(f)
        warning('[%s] No *_NGS.mat found. Skipping.', caseID);
        continue;
    end

    % If multiple matches, pick the newest file
    [~, newestIdx] = max([f.datenum]);
    matPath = fullfile(f(newestIdx).folder, f(newestIdx).name);

    S = load(matPath);

    if ~isfield(S, 'NGS_table')
        warning('[%s] %s does not contain variable "NGS_table". Skipping.', caseID, matPath);
        continue;
    end

    NGS_table = S.NGS_table;

    if ~istable(NGS_table)
        warning('[%s] NGS_table is not a table in %s. Skipping.', caseID, matPath);
        continue;
    end

    if isempty(NGS_table.Properties.RowNames) || isempty(NGS_table.Properties.VariableNames)
        warning('[%s] NGS_table missing RowNames/VariableNames in %s. Skipping.', caseID, matPath);
        continue;
    end

    rn = NGS_table.Properties.RowNames;
    cn = NGS_table.Properties.VariableNames;

    % Convert to numeric matrix
    try
        V = table2array(NGS_table);
    catch
        warning('[%s] Failed to convert NGS_table to array in %s. Skipping.', caseID, matPath);
        continue;
    end

    if ~isnumeric(V)
        warning('[%s] NGS_table values are not numeric in %s. Skipping.', caseID, matPath);
        continue;
    end

    % Flatten into a struct
    rowStruct = struct();
    rowStruct.ID = string(caseID);

    for i = 1:numel(rn)
        for j = 1:numel(cn)
            fname = matlab.lang.makeValidName(sprintf('%s_%s', rn{i}, cn{j}));
            rowStruct.(fname) = V(i, j);
            allVarNames{end+1} = fname; %#ok<SAGROW>
        end
    end

    allRows{end+1} = rowStruct; %#ok<SAGROW>
    fprintf('[%s] Loaded: %s\n', caseID, matPath);
end

if isempty(allRows)
    error('No valid cases were loaded. Please check folder structure and *_NGS*.mat files.');
end

allVarNames = unique(allVarNames, 'stable');

% -------------------------------------------------------------------------
% 3) Build final table (missing values -> NaN)
% -------------------------------------------------------------------------
nCases = numel(allRows);

T = table('Size', [nCases, 1 + numel(allVarNames)], ...
          'VariableTypes', [{ 'string' }, repmat({'double'}, 1, numel(allVarNames))], ...
          'VariableNames', [{'ID'}, allVarNames]);

for k = 1:nCases
    rs = allRows{k};
    T.ID(k) = rs.ID;

    for v = 1:numel(allVarNames)
        fn = allVarNames{v};
        if isfield(rs, fn)
            T{k, 1+v} = rs.(fn);
        else
            T{k, 1+v} = NaN;
        end
    end
end

% Sort rows by numeric ID if possible
nRows = height(T);
idNums = nan(nRows,1);

for i = 1:nRows
    % Extract numeric part from ID
    numStr = regexp(T.ID(i), '\d+', 'match');
    if ~isempty(numStr)
        idNums(i) = str2double(numStr{1});
    else
        idNums(i) = NaN;
    end
end

% Sort: numeric first (ascending), then NaNs at bottom
[~, ord] = sortrows([isnan(idNums), idNums], [1 2]);

T = T(ord,:);

% -------------------------------------------------------------------------
% 4) Save MAT + Excel
% -------------------------------------------------------------------------
save(outMat, 'T', 'rootDir');
writetable(T, outXlsx, 'FileType', 'spreadsheet');

fprintf('\nDone.\nSaved MAT : %s\nSaved XLSX: %s\n', outMat, outXlsx);