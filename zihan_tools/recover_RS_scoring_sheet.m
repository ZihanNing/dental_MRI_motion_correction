%% recover_RS_scoring_sheet_v2.m
% Recover radiologist scoring sheet from anonymised A/B/C order
% back to the true reconstruction categories.
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-25

clear; clc;

% -------------------------------------------------------------------------
% Paths
% -------------------------------------------------------------------------
rootDir   = '/home/zn23/Data/ddMRI/';
forRSDir  = fullfile(rootDir, 'for_RS');

keyXlsx   = fullfile(forRSDir, 'RS_randomisation_key.xlsx');
scoreXlsx = fullfile(rootDir, 'scoring_sheet_260327_JD.xlsx');

outXlsx   = fullfile(rootDir, 'scoring_sheet_JD_recovered.xlsx');
outMat    = fullfile(rootDir, 'scoring_sheet_JD_recovered.mat');

% -------------------------------------------------------------------------
% Check files
% -------------------------------------------------------------------------
if ~isfile(keyXlsx)
    error('Randomisation key not found:\n%s', keyXlsx);
end
if ~isfile(scoreXlsx)
    error('Scoring sheet not found:\n%s', scoreXlsx);
end

% -------------------------------------------------------------------------
% Read key
% -------------------------------------------------------------------------
T_key = readtable(keyXlsx, 'TextType', 'string');
T_key.ID = strtrim(string(T_key.ID));

fprintf('Loaded key   : %s\n', keyXlsx);
fprintf('Loaded score : %s\n', scoreXlsx);

% -------------------------------------------------------------------------
% Read score sheet as raw cell array
% This is more robust for Excel sheets with merged multi-row headers.
% -------------------------------------------------------------------------
raw = readcell(scoreXlsx);

if size(raw,2) < 10
    error('Scoring sheet seems to have fewer than 10 columns.');
end

% Expected layout:
% row 1: merged upper/lower header or blank
% row 2: actual column names
% row 3+: data
%
% We will detect the first valid data row automatically:
% valid data row = first column numeric, second column string like P01-01/HV01-01

startRow = [];
for r = 1:size(raw,1)
    c1 = raw{r,1};
    c2 = raw{r,2};

    isIDNumeric = isnumeric(c1) && isfinite(c1);
    isUnqIDText = (ischar(c2) || isstring(c2)) && ~isempty(strtrim(string(c2))) ...
        && ~contains(lower(strtrim(string(c2))), "unqid");

    if isIDNumeric && isUnqIDText
        startRow = r;
        break;
    end
end

if isempty(startRow)
    error('Could not detect the first data row in scoring sheet.');
end

fprintf('Detected first data row at Excel row %d\n', startRow);

% -------------------------------------------------------------------------
% Recover rows
% -------------------------------------------------------------------------
rows = {};

for r = startRow:size(raw,1)

    % Stop if first column empty
    if isempty(raw{r,1}) || (isstring(raw{r,1}) && strlength(raw{r,1})==0)
        continue;
    end

    % Column mapping
    % 1: ID
    % 2: UnqID
    % 3: Seq Type
    % 4-6: upper A/B/C
    % 7-9: lower A/B/C
    % 10: Best Image

    idCell = raw{r,1};
    unqID  = string(raw{r,2});
    seqType = string(raw{r,3});

    % skip non-data rows
    if ~(isnumeric(idCell) && isfinite(idCell))
        continue;
    end

    scoreID = string(idCell);   % this should match T_key.ID

    idxKey = find(T_key.ID == scoreID, 1);
    if isempty(idxKey)
        warning('No randomisation key found for scoring-sheet ID "%s" (row %d). Skipping.', scoreID, r);
        continue;
    end

    upperA = raw{r,4};
    upperB = raw{r,5};
    upperC = raw{r,6};

    lowerA = raw{r,7};
    lowerB = raw{r,8};
    lowerC = raw{r,9};

    bestRaw = "";
    if size(raw,2) >= 10 && ~isempty(raw{r,10})
        bestRaw = lower(strtrim(string(raw{r,10})));
    end

    mapA = string(T_key.A_category(idxKey));
    mapB = string(T_key.B_category(idxKey));
    mapC = string(T_key.C_category(idxKey));

    [upper_NoMoCo, upper_MoCoFullFOV, upper_MoCoFused] = ...
        assignRecoveredScores(mapA, mapB, mapC, upperA, upperB, upperC);

    [lower_NoMoCo, lower_MoCoFullFOV, lower_MoCoFused] = ...
        assignRecoveredScores(mapA, mapB, mapC, lowerA, lowerB, lowerC);

    switch bestRaw
        case "a"
            bestImageName = mapA;
        case "b"
            bestImageName = mapB;
        case "c"
            bestImageName = mapC;
        otherwise
            bestImageName = missing;
    end

    S = struct();
    S.ID                = idCell;
    S.UnqID             = unqID;
    S.SeqType           = seqType;

    S.upper_NoMoCo      = upper_NoMoCo;
    S.upper_MoCoFullFOV = upper_MoCoFullFOV;
    S.upper_MoCoFused   = upper_MoCoFused;

    S.lower_NoMoCo      = lower_NoMoCo;
    S.lower_MoCoFullFOV = lower_MoCoFullFOV;
    S.lower_MoCoFused   = lower_MoCoFused;

    S.BestImageName     = bestImageName;

    rows{end+1} = S; %#ok<SAGROW>
end

if isempty(rows)
    error('No rows were successfully recovered.');
end

T_out = struct2table([rows{:}]', 'AsArray', true);

% optional sorting by numeric ID
[~, ord] = sort(T_out.ID);
T_out = T_out(ord,:);

writetable(T_out, outXlsx);
save(outMat, 'T_out');

fprintf('\nRecovered scoring sheet saved to:\n%s\n', outXlsx);
fprintf('Recovered MAT saved to:\n%s\n', outMat);

% -------------------------------------------------------------------------
% Local function
% -------------------------------------------------------------------------
function [scoreNoMoCo, scoreFullFOV, scoreFused] = ...
    assignRecoveredScores(catA, catB, catC, scoreA, scoreB, scoreC)

    scoreNoMoCo = missing;
    scoreFullFOV = missing;
    scoreFused = missing;

    cats   = {char(catA), char(catB), char(catC)};
    scores = {scoreA, scoreB, scoreC};

    for ii = 1:3
        switch cats{ii}
            case 'No-MoCo'
                scoreNoMoCo = scores{ii};
            case 'MoCo-fullFOV'
                scoreFullFOV = scores{ii};
            case 'MoCo-fused'
                scoreFused = scores{ii};
            otherwise
                warning('Unknown category: %s', cats{ii});
        end
    end
end