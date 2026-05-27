%% compute_radiologist_icc.m
% Compute inter-rater ICC between two radiologists for three reconstructions
%
% The script reads the first row of the Excel file as headers exactly as written,
% then performs robust matching by normalising header strings
% (ignoring spaces, hyphens, underscores, etc.).
%
% ICC computed: ICC(2,1)
% two-way random-effects, absolute agreement, single measurement
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 04-Apr-2026

clear; clc;

%% Input file
xlsxFile = '/home/zn23/Data/ddMRI/allRS.xlsx';

%% Read raw Excel content
[~, ~, raw] = xlsread(xlsxFile);

if isempty(raw) || size(raw,1) < 2
    error('Excel file is empty or does not contain enough rows.');
end

headers = raw(1, :);
dataRaw = raw(2:end, :);

disp('Headers read from Excel:')
disp(headers')

%% Define requested column pairs
pairs = {
    'No-MoCo',       'Rad1 - No-MoCo',       'Rad2 - No-MoCo';
    'MoCo-fullFOV',  'Rad1 - MoCo-fullFOV',  'Rad2 - MoCo-fullFOV';
    'MoCo-fused',    'Rad1 - MoCo-fused',    'Rad2 - MoCo-fused';
    };

%% Prepare results
results = table('Size', [size(pairs,1), 6], ...
    'VariableTypes', {'string','double','double','double','double','double'}, ...
    'VariableNames', {'Reconstruction','N','ICC','F','df1','df2'});

fprintf('=============================================\n');
fprintf('Inter-rater ICC results (ICC(2,1))\n');
fprintf('File: %s\n', xlsxFile);
fprintf('=============================================\n\n');

%% Compute ICC
for i = 1:size(pairs,1)
    reconName = pairs{i,1};
    header1   = pairs{i,2};
    header2   = pairs{i,3};

    idx1 = find_header_index(headers, header1);
    idx2 = find_header_index(headers, header2);

    if isempty(idx1)
        error('Cannot find column header: %s', header1);
    end
    if isempty(idx2)
        error('Cannot find column header: %s', header2);
    end

    fprintf('Matched "%s" --> "%s"\n', header1, headers{idx1});
    fprintf('Matched "%s" --> "%s"\n', header2, headers{idx2});

    x1_raw = dataRaw(:, idx1);
    x2_raw = dataRaw(:, idx2);

    x1 = cellfun(@convert_to_numeric, x1_raw);
    x2 = cellfun(@convert_to_numeric, x2_raw);

    validMask = ~(isnan(x1) | isnan(x2));
    X = [x1(validMask), x2(validMask)];

    if size(X,1) < 2
        warning('Not enough valid paired scores for %s. Skipping.', reconName);
        results.Reconstruction(i) = string(reconName);
        results.N(i)   = size(X,1);
        results.ICC(i) = NaN;
        results.F(i)   = NaN;
        results.df1(i) = NaN;
        results.df2(i) = NaN;
        fprintf('\n');
        continue;
    end

    stats = compute_icc_2_1(X);

    results.Reconstruction(i) = string(reconName);
    results.N(i)   = stats.n;
    results.ICC(i) = stats.ICC;
    results.F(i)   = stats.F;
    results.df1(i) = stats.df1;
    results.df2(i) = stats.df2;

    fprintf('Reconstruction: %s\n', reconName);
    fprintf('  N   = %d\n', stats.n);
    fprintf('  ICC = %.4f\n', stats.ICC);
    fprintf('  F   = %.4f\n', stats.F);
    fprintf('  df1 = %d\n', stats.df1);
    fprintf('  df2 = %d\n\n', stats.df2);
end

%% Show and save results
disp(results);

outFile = '/home/zn23/Data/ddMRI/allRS_ICC_results.csv';
writetable(results, outFile);
fprintf('Results saved to: %s\n', outFile);


%% ===== local functions =====

function idx = find_header_index(headers, targetHeader)
% Find header index using normalised matching:
% ignore spaces, hyphens, underscores and non-alphanumeric characters.

    targetNorm = normalize_header(targetHeader);
    idx = [];

    for ii = 1:numel(headers)
        thisHeader = headers{ii};

        if ischar(thisHeader) || isstring(thisHeader)
            thisNorm = normalize_header(char(thisHeader));
            if strcmp(thisNorm, targetNorm)
                idx = ii;
                return;
            end
        end
    end
end

function s = normalize_header(s)
% Convert header to lower case and remove all non-alphanumeric characters

    s = lower(char(s));
    s = regexprep(s, '[^a-z0-9]', '');
end

function val = convert_to_numeric(x)
% Convert Excel cell content to numeric value

    if isempty(x)
        val = NaN;
    elseif isnumeric(x)
        val = x;
    elseif ischar(x)
        val = str2double(x);
    elseif isstring(x)
        val = str2double(char(x));
    else
        val = NaN;
    end
end

function stats = compute_icc_2_1(X)
% Compute ICC(2,1): two-way random-effects, absolute agreement, single measurement

    [n, k] = size(X);

    if k ~= 2
        error('This script currently expects exactly 2 raters.');
    end
    if n < 2
        error('Not enough valid cases to compute ICC.');
    end

    rowMeans  = mean(X, 2);
    colMeans  = mean(X, 1);
    grandMean = mean(X(:));

    SSR = k * sum((rowMeans - grandMean).^2);
    SSC = n * sum((colMeans - grandMean).^2);
    SSE = sum(sum((X - rowMeans(:,ones(1,k)) - colMeans(ones(n,1),:) + grandMean).^2));

    dfR = n - 1;
    dfC = k - 1;
    dfE = (n - 1) * (k - 1);

    MSR = SSR / dfR;
    MSC = SSC / dfC;
    MSE = SSE / dfE;

    ICC = (MSR - MSE) / (MSR + (k - 1)*MSE + k*(MSC - MSE)/n);
    F = MSR / MSE;

    stats = struct();
    stats.ICC = ICC;
    stats.F   = F;
    stats.df1 = dfR;
    stats.df2 = dfE;
    stats.n   = n;
    stats.k   = k;
    stats.MSR = MSR;
    stats.MSC = MSC;
    stats.MSE = MSE;
end