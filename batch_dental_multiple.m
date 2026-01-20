%% batchRecon_dental_dual.m
% Batch processing script for ddMRI_new cases
%
% For each case folder (e.g. ddMRI_new/1/):
%   1) Find ALL .dat files
%   2) Process them ONE BY ONE
%   3) Move to the next case
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

clear; clc;

% ---- USER SETTINGS ----
rootFolder  = '/data/gadgetron/matlab_study/dental_paper_all';
studiesFile = fullfile('./Studies-deploy', 'studies.m');
numCases    = 17;
caseList    = 3:17;   % subset if needed

% >>> NEW: sequence selection <<<
% Leave EMPTY {} to reconstruct ALL sequences (default behaviour)
% Otherwise choose one or more of:
%   {'MPRAGE'}, {'PDwSPACE'}, {'T2wSPACE'}, or combinations
seqSelect = {'MPRAGE','T2wSPACE'};

for caseIdx = caseList

    clc;
    fprintf('\n=========================================\n');
    fprintf('Processing Case %d\n', caseIdx);
    fprintf('=========================================\n');

    % ---- Case folder ----
    caseFolder = [fullfile(rootFolder, num2str(caseIdx)) '/'];

    if ~isfolder(caseFolder)
        warning('Case folder not found: %s. Skipping...', caseFolder);
        continue;
    end

    % ---- 1) Find ALL .dat files in this case ----
    datFilesAll = dir(fullfile(caseFolder, '*.dat'));

    % ---- Optional sequence filtering ----
    if isempty(seqSelect)
        datFiles = datFilesAll;  % default: keep all
    else
        keepMask = false(size(datFilesAll));
        for k = 1:numel(datFilesAll)
            namek = datFilesAll(k).name;
            for s = 1:numel(seqSelect)
                if contains(namek, seqSelect{s}, 'IgnoreCase', true)
                    keepMask(k) = true;
                    break;
                end
            end
        end
        datFiles = datFilesAll(keepMask);
    end


    if isempty(datFiles)
        warning('No .dat files found in %s. Skipping case...', caseFolder);
        continue;
    end

    fprintf('Found %d .dat files in Case %d\n', numel(datFiles), caseIdx);

    % ---- 2) Process each .dat file independently ----
    for fIdx = 1:numel(datFiles)

        [~, baseName, ~] = fileparts(datFiles(fIdx).name);

        % ---- Detect sequence type (for logging only) ----
        seqType = 'Unknown';
        if contains(baseName, 'PDwSPACE', 'IgnoreCase', true)
            seqType = 'PDwSPACE';
        elseif contains(baseName, 'T2wSPACE', 'IgnoreCase', true)
            seqType = 'T2wSPACE';
        elseif contains(baseName, 'MPRAGE', 'IgnoreCase', true)
            seqType = 'MPRAGE';
        end

        fprintf('\n--- Case %d: file %d / %d ---\n', ...
                caseIdx, fIdx, numel(datFiles));
        fprintf('File name: %s\n', baseName);
        fprintf('Seq type : %s\n', seqType);

        % ---- Generate studies.m for THIS file only ----
        fid = fopen(studiesFile, 'wt');
        if fid == -1
            error('Cannot open %s for writing.', studiesFile);
        end

        fprintf(fid, '%% Auto-generated studies.m\n');
        fprintf(fid, '%% Case %d, file %s\n\n', caseIdx, baseName);

        fprintf(fid, 'pathData{1}=''%s'';\n', caseFolder);
        fprintf(fid, 'pathIn{1}=pathData{1};\n');
        fprintf(fid, 'pathRef{1}=strcat(pathData{1},''/Re-Se'');\n');
        fprintf(fid, 'pathRemote{1}=pathData{1};\n\n');

        fprintf(fid, '%% ACQUISITIONS DATA\n');
        fprintf(fid, 'fileIn{1}{1}=''%s'';\n\n', baseName);

        fprintf(fid, '%% REFERENCE DATA\n');
        fprintf(fid, 'refIn  = fillCell(fileIn, '''');\n');
        fprintf(fid, 'refBIn = fillCell(fileIn, '''');\n\n');

        fprintf(fid, '%% B0 and B1\n');
        fprintf(fid, 'B0In = fillCell(fileIn, '''');\n');
        fprintf(fid, 'B1In = fillCell(fileIn, '''');\n\n');

        fprintf(fid, '%% PILOT TONE FLAG\n');
        fprintf(fid, 'isPT     = fillCell(fileIn, 0);\n');
        fprintf(fid, 'facFOVTh = fillCell(fileIn, -.45);\n\n');

        fprintf(fid, '%% DATA CONVERSION SPECIFICATIONS\n');
        fprintf(fid, 'noiseFile      = fillCell(fileIn, refIn{1}{1});\n');
        fprintf(fid, 'supportReadout = fillCell(fileIn, []);\n');
        fprintf(fid, 'resRec         = fillCell(fileIn, []);\n');
        fprintf(fid, 'RDesired       = fillCell(fileIn, []);\n');

        fclose(fid);

        % ---- Save loop state (deployRecon_dental may clear) ----
        stateFile = 'batchRecon_dental_state.mat';
        currDir   = pwd;

        save(stateFile, ...
             'rootFolder', 'studiesFile', 'numCases', 'caseList', ...
             'caseIdx', 'fIdx', 'currDir', ...
             'datFiles', 'caseFolder', 'seqSelect');


        % ---- Run reconstruction ----
        try
            fprintf('Running deployRecon_dental...\n');
            deployRecon_dental;
            fprintf('Reconstruction completed.\n');
        catch ME
            fprintf(2, ...
                'Error in Case %d, file %s:\n%s\n', ...
                caseIdx, baseName, ME.message);
        end

        % ---- Restore state ----
        if exist('batchRecon_dental_state.mat', 'file')
            load('batchRecon_dental_state.mat', ...
                 'rootFolder', 'studiesFile', 'numCases', 'caseList', ...
                 'caseIdx', 'fIdx', 'currDir', ...
                 'datFiles', 'caseFolder', 'seqSelect');

            try
                cd(currDir);
            catch
                % If directory no longer exists, just ignore
            end
            delete('batchRecon_dental_state.mat');
        else
            warning('State file %s not found. Loop variables may be lost.', stateFile);
        end

    end % loop over dat files

end % loop over cases

fprintf('\nAll requested cases processed.\n');
