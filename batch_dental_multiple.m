%% batchRecon_dental_multiple.m
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
addpath(genpath(pwd))

% ---- USER SETTINGS ----
rootFolder  = [pwd,'/Studies-deploy'];
studiesFile = fullfile('./Studies-deploy', 'studies.m');
numCases    = 1;
caseList    = 14;   % subset if needed

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
         
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% SENSE reconstruction 
        % prepare the image for recon
        % predict coil sensitivity map by ESPIRIT
        % run SENSE and saved the image 
        %       in $caseFolder/An-Aq
        %       file name:
        %       [erase(datFiles(fIdx).name,".dat"),'_Aq_womsk.nii']
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Call python for preprocessing for segmentation
        % the python script are within /Python sub-folder
        % preprocessing to generate the cases for nnUNet segmentation
        % the steps contains
        %   crop and downsampling 
        %       by ddMRI_cropds_recover.py and generate
        %       '_Aq_womsk_cropds.nii.gz' image and a corresponding json
        %       file
        %       crop to FOV: "RL": 144.0, "AP": 160.0, "SI": 150.0
        %       downsample to 2mm iso
        %   after crop and downsample, norm the intenstiy to [0 1000] via
        %   robust normalization by normalize_intensity_robust.py
        %   and generate '_Aq_womsk_cropds_norm.nii.gz' (outputNii_norm)
        %   which ready to be the testing case for nnUNet seg prediction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -------- Python preprocessing: crop + downsample --------
        % Python info
        PY_CROPDS = fullfile(currDir, 'Python', 'ddMRI_cropds_recover.py');
        PY_NORM = fullfile(currDir, 'Python', 'normalize_intensity_robust.py');

        % input files: SENSE reconed image
        inputNii = fullfile(caseFolder, 'An-Aq', ...
            [erase(datFiles(fIdx).name, '.dat'), '_Aq_womsk.nii']);

        % output files
        pathPrep = fullfile(caseFolder, 'An-Aq/seg'); if ~exist( pathPrep,'dir');mkdir(pathPrep);end
        outputNii_cropds = fullfile(pathPrep, ...
            [erase(datFiles(fIdx).name, '.dat'), '_Aq_womsk_cropds.nii.gz']);
        outputNii_norm = fullfile(pathPrep, ...
            [erase(datFiles(fIdx).name, '.dat'), '_Aq_womsk_cropds_norm.nii.gz']);

        % --- Call python for crop and downsampling
        CONDA = '/home/zn23/anaconda3/bin/conda';
        ENVNAME = 'nnunetv2';

        cmd = sprintf(['"%s" run -n %s python "%s" "%s" "%s" --mode generic'], ...
            CONDA, ENVNAME, PY_CROPDS, inputNii, outputNii_cropds);
        [status, cmdout] = system(cmd);
        
         % --- Call python for normalization
        cmd = sprintf(['"%s" run -n %s python "%s" "%s" "%s" --plow 0.5 --phigh 99.5'], ...
            CONDA, ENVNAME, PY_NORM, outputNii_cropds, outputNii_norm);
        [status, cmdout] = system(cmd);

        if status ~= 0
            fprintf(2, 'Python preprocessing failed:\n%s\n', cmdout);
            error('ddMRI_cropds_recover.py failed');
        else
            fprintf('Python preprocessing finished:\n%s\n', cmdout);
        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Use nnUNet for segmentation prediction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Input file: SENSE reconed image after preprocessing (cropds +
        % norm)
        inputNiiNorm = outputNii_norm;  

        % nnUNet temp folders
        nnUNetInDir  = fullfile(caseFolder, 'An-Aq', 'nnunet_in');
        nnUNetOutDir = fullfile(caseFolder, 'An-Aq', 'nnunet_out');

        if ~exist(nnUNetInDir, 'dir'), mkdir(nnUNetInDir); end
        if ~exist(nnUNetOutDir, 'dir'), mkdir(nnUNetOutDir); end
        
        % nnUNet requires *_0000.nii.gz
        nnUNetInputFile = fullfile(nnUNetInDir, 'ddMRI_001_0000.nii.gz');
        copyfile(inputNiiNorm, nnUNetInputFile);
        
        % call trained network for segmentation
        NNUNET_BASE = '/home/zn23/nnUNet';
        
        % segment of upper and lower teeth (trained network 2)
        cmd = sprintf([ ...
            '"%s" run -n %s ' ...
            'bash -lc '' ' ...
            'export nnUNet_raw="%s/nnUNet_raw"; ' ...
            'export nnUNet_preprocessed="%s/nnUNet_preprocessed"; ' ...
            'export nnUNet_results="%s/nnUNet_results"; ' ...
            'nnUNetv2_predict ' ...
            '-i "%s" -o "%s" ' ...
            '-d 2 -c 3d_fullres -f 0 1 2 3 4 ' ...
            ''''], ...
            CONDA, ENVNAME, ...
            NNUNET_BASE, NNUNET_BASE, NNUNET_BASE, ...
            nnUNetInDir, nnUNetOutDir);
        [status, cmdout] = system(cmd);

        if status ~= 0
            fprintf(2, 'nnUNet prediction failed:\n%s\n', cmdout);
            error('nnUNetv2_predict failed');
        else
            fprintf('nnUNet prediction finished:\n%s\n', cmdout);
        end
        
        % copy the result back and free temp
        nnUNetOutputFile = fullfile(nnUNetOutDir, 'ddMRI_001.nii.gz');
        nnUNetMaskteeth = fullfile(pathPrep, ...
            [erase(datFiles(fIdx).name, '.dat'), '_msk_teeth.nii.gz']);
        copyfile(nnUNetOutputFile, nnUNetMaskteeth);
        rmdir(nnUNetInDir, 's'); 
        rmdir(nnUNetOutDir, 's'); 




        

    end % loop over dat files

end % loop over cases

fprintf('\nAll requested cases processed.\n');
