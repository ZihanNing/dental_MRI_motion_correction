%% batchRecon_dental.m
% Batch processing script for 16 cases in ./Studies-deploy/neonates_mp2rage/
% Auto-generates studies.m for each case and runs deployRecon_mp2rage.m
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2025-09-25

clear; clc;

% ---- USER SETTINGS ----
% rootFolder   = './Studies-deploy/neonates_mp2rage/';
% numCases     = 15;
% studiesFile  = fullfile('./Studies-deploy', 'studies.m');

for caseIdx = 9
    
    clc
    rootFolder   = '/data/gadgetron/matlab_study/ddMRI';
    studiesFile  = fullfile('./Studies-deploy', 'studies.m');

    fprintf('\n=========================================\n');
    fprintf('Processing Case %d of %d\n', caseIdx, 15);
    fprintf('=========================================\n');

    % Locate case folder (ensure trailing slash)
    caseFolder = [fullfile(rootFolder, num2str(caseIdx)) '/'];

    % Find .dat files in this folder
    datFiles = dir(fullfile(caseFolder, '*.dat'));
    if isempty(datFiles)
        warning('No .dat file found in %s, skipping...', caseFolder);
        continue;
    end

    % Prepare list of file names (without .dat extension)
    fileNames = cell(1, numel(datFiles));
    for k = 1:numel(datFiles)
        [~, baseName, ~] = fileparts(datFiles(k).name);
        fileNames{k} = baseName;
    end

    % ---- Generate new studies.m ----
    fid = fopen(studiesFile, 'wt');
    if fid == -1
        error('Cannot open studies.m for writing.');
    end

    fprintf(fid, '%% Auto-generated studies.m for Case %d\n', caseIdx);
    fprintf(fid, 'pathData{1}=''%s'';\n', caseFolder);
    fprintf(fid, 'pathIn{1}=pathData{1};\n');
    fprintf(fid, 'pathRef{1}=strcat(pathData{1},''/Re-Se'');\n');
    fprintf(fid, 'pathRemote{1}=pathData{1};\n\n');

    fprintf(fid, '%% ACQUISITIONS DATA\n');
    for k = 1:numel(fileNames)
        fprintf(fid, 'fileIn{1}{%d}=''%s'';\n', k, fileNames{k});
    end
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

    % ---- Verification ----
    fprintf('Generated new studies.m:\n');
    type(studiesFile);
    
    pause(5);

    % ---- Run reconstruction ----
    try
        deployRecon_dental;
    catch ME
        fprintf(2, 'Error in Case %d: %s\n', caseIdx, ME.message);
    end

end

fprintf('\nAll cases processed.\n');