%% batchRecon_dental_dual.m
% Batch processing script for ddMRI_new cases
% Each case folder (e.g. ddMRI_new/1/) may contain:
%   - location.txt
%   - PDwSPACE *.dat   (e.g. ..._PDwSPACE_...)
%   - T2wSPACE *.dat   (e.g. ..._T2wSPACE_...)
%
% For each .dat file in a case folder, this script:
%   1) Auto-generates ./Studies-deploy/studies.m
%      containing ONLY that single acquisition
%   2) Calls deployRecon_dental
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2025-12-08

clear; clc;

% ---- USER SETTINGS ----
rootFolder  = '/data/gadgetron/matlab_study/ddMRI_new';   % parent folder of numbered cases
studiesFile = fullfile('./Studies-deploy', 'studies.m');                              % studies.m used by deployRecon_dental
numCases    = 8;                                                                      % total number of cases (e.g. 1..15)
caseList    = [5 6 7 8];                                                             % modify if you want a subset, e.g. [1 3 5]

for caseIdx = caseList

    clc;
    fprintf('\n=========================================\n');
    fprintf('Processing Case %d of %d\n', caseIdx, numCases);
    fprintf('=========================================\n');

    % Ensure trailing slash in caseFolder (some downstream code may rely on it)
    caseFolder = [fullfile(rootFolder, num2str(caseIdx)) '/'];

    if ~isfolder(caseFolder)
        warning('Case folder not found: %s. Skipping...', caseFolder);
        continue;
    end

    % Find all .dat files in this case folder
    datFiles = dir(fullfile(caseFolder, '*.dat'));
    if isempty(datFiles)
        warning('No .dat files found in %s. Skipping case...', caseFolder);
        continue;
    end

    % Loop over each .dat file and reconstruct independently
    for fIdx = 1:numel(datFiles)

        % Recompute paths in case they were modified anywhere
        rootFolder  = '/data/gadgetron/matlab_study/ddMRI_new';
        caseFolder  = [fullfile(rootFolder, num2str(caseIdx)) '/'];
        datFiles    = dir(fullfile(caseFolder, '*.dat'));

        [~, baseName, ~] = fileparts(datFiles(fIdx).name);

        % Optional: detect sequence type from name (for logging only)
        seqType = 'Unknown';
        if contains(baseName, 'PDwSPACE', 'IgnoreCase', true)
            seqType = 'PDwSPACE';
        elseif contains(baseName, 'T2wSPACE', 'IgnoreCase', true)
            seqType = 'T2wSPACE';
        elseif contains(baseName, 'MPRAGE', 'IgnoreCase', true)
            seqType = 'MPRAGE';
        end

        fprintf('\n--- Case %d: reconstructing file %d of %d ---\n', ...
                caseIdx, fIdx, numel(datFiles));
        fprintf('File name: %s\n', baseName);
        fprintf('Seq type : %s\n', seqType);

        % ---- Generate new studies.m for THIS single .dat file ----
        fid = fopen(studiesFile, 'wt');
        if fid == -1
            error('Cannot open %s for writing.', studiesFile);
        end

        fprintf(fid, '%% Auto-generated studies.m for Case %d, file %s\n', caseIdx, baseName);
        fprintf(fid, 'pathData{1}=''%s'';\n', caseFolder);
        fprintf(fid, 'pathIn{1}=pathData{1};\n');
        fprintf(fid, 'pathRef{1}=strcat(pathData{1},''/Re-Se'');\n');
        fprintf(fid, 'pathRemote{1}=pathData{1};\n\n');

        fprintf(fid, '%% ACQUISITIONS DATA\n');
        % Only one acquisition for this reconstruction:
        fprintf(fid, 'fileIn{1}{1}=''%s'';\n', baseName);
        fprintf(fid, '%%Keep adding files if needed\n\n');

        fprintf(fid, '%% REFERENCE DATA\n');
        fprintf(fid, 'refIn = fillCell(fileIn, '''');\n');
        fprintf(fid, 'refBIn = fillCell(fileIn, '''');\n\n');

        fprintf(fid, '%% B0 and B1\n');
        fprintf(fid, 'B0In = fillCell(fileIn, '''');\n');
        fprintf(fid, 'B1In = fillCell(fileIn, '''');\n\n');

        fprintf(fid, '%% PILOT TONE FLAG - leave untouched\n');
        fprintf(fid, 'isPT = fillCell(fileIn, 0);\n');
        fprintf(fid, 'facFOVTh = fillCell(fileIn, -.45);\n\n');

        fprintf(fid, '%% DATA CONVERSION SPECIFICATIONS - leave untouched\n');
        fprintf(fid, 'noiseFile = fillCell(fileIn,refIn{1}{1});\n');
        fprintf(fid, 'supportReadout = fillCell(fileIn,[]);\n');
        fprintf(fid, 'resRec = fillCell(fileIn,[]);\n');
        fprintf(fid, 'RDesired = fillCell(fileIn,[]);\n');

        fclose(fid);

        % ---- Save loop state before calling deployRecon_dental ----
        % (because deployRecon_dental starts with "clear")
        stateFile = 'batchRecon_dental_state.mat';
        currDir   = pwd;  % also save current directory in case deployRecon_dental does cd

        save(stateFile, ...
             'rootFolder', 'studiesFile', 'numCases', 'caseList', ...
             'caseIdx', 'fIdx', 'currDir');

        % ---- Run reconstruction for this single file ----
        try
            fprintf('Running deployRecon_dental for %s ...\n', baseName);
            deployRecon_dental;   % this may contain "clear"
            fprintf('Reconstruction completed!\n');
        catch ME
            fprintf(2, 'Error in Case %d, file %s: %s\n', caseIdx, baseName, ME.message);
        end

        % ---- Restore loop state after deployRecon_dental returns ----
        % (even if it has executed "clear" in its script body)
        if exist('batchRecon_dental_state.mat', 'file')
            load('batchRecon_dental_state.mat', ...
                 'rootFolder', 'studiesFile', 'numCases', 'caseList', ...
                 'caseIdx', 'fIdx', 'currDir');
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
