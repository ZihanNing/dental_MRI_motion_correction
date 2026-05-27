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

% This is in use 24-Mar-2026 by Zihan

clear; clc;
addpath(genpath(pwd))

%% ---- Start logging ----
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

logFile = fullfile(pwd, ...
    sprintf('log_PIPE_%s.txt', timestamp));

diary(logFile);
diary on;


if exist(logFile, 'file')
    delete(logFile);   % optional: remove old log
end

diary(logFile);
diary on;

fprintf('==== PIPE started: %s ====\n', datestr(now));

%% ---- USER SETTINGS ----
rootFolder  = './Studies-deploy';
studiesFile = fullfile('./Studies-deploy', 'studies.m');
numCases    = 1;
caseList    = [1];   % subset if needed, e.g., [15]

% >>> NEW: sequence selection <<<
% Leave EMPTY {} to reconstruct ALL sequences (default behaviour)
% Otherwise choose one or more of:
%   {'MPRAGE'}, {'PDwSPACE'}, {'T2wSPACE'}, or combinations
seqSelect = {'MPRAGE','T2wSPACE','PDwSPACE'};

for caseIdx = caseList

    clc;
    fprintf('\n=========================================\n');
    fprintf('Processing Case %d\n', caseIdx);
    fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
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
            fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
            deployRecon_dental_SENSE;
            fprintf('Reconstruction completed.\n');
            fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
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
        %%% Use nnUNet for segmentation prediction (upper and lower teeth)
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
            '-d 3 -c 3d_fullres -f 0 1 2 3 4 ' ...
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

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Use nnUNet for segmentation prediction (head seg)
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
            '-d 4 -c 3d_fullres -f 0 1 2 3 4 ' ...
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
            [erase(datFiles(fIdx).name, '.dat'), '_msk_head.nii.gz']);
        copyfile(nnUNetOutputFile, nnUNetMaskteeth);
        rmdir(nnUNetInDir, 's'); 
        rmdir(nnUNetOutDir, 's'); 
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Recover the masks (teeth & head) to original resol & FOV
        % generate large masks:
        % - head: /An-Aq/seg/...msk_head_fullresol.nii.gz
        % - teeth: /An-Aq/seg/...msk_teeth_fullresol.nii.gz
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Required variables
        % pathPrep = '...'; % folder with masks + json (set in above
        % section)
        % CONDA   = '/home/zn23/anaconda3/bin/conda';
        % ENVNAME = 'nnunetv2';
        
        PY_RECOVER = fullfile(currDir, 'Python', 'restore_mask_from_json.py');

        % find JSON (use most recent if multiple)
        jsonList = dir(fullfile(pathPrep, '*.json'));
        if isempty(jsonList)
            error('No json found in: %s', pathPrep);
        end
        [~, idxJ] = max([jsonList.datenum]);
        jsonFile = fullfile(pathPrep, jsonList(idxJ).name);

        % restore head
        headList = dir(fullfile(pathPrep, '*msk_head.nii.gz'));
        if isempty(headList)
            warning('No head mask (*msk_head.nii.gz) found in: %s', pathPrep);
        else
            [~, idxH] = max([headList.datenum]);
            headMask = fullfile(pathPrep, headList(idxH).name);
            outHead  = strrep(headMask, '.nii.gz', '_fullresol.nii.gz');

            cmd = sprintf(['"%s" run -n %s python "%s" --mask "%s" --json "%s" --out "%s"'], ...
                CONDA, ENVNAME, PY_RECOVER, headMask, jsonFile, outHead);

            fprintf('\n[restore head]\nmask: %s\njson: %s\nout : %s\n', headMask, jsonFile, outHead);
            [status, cmdout] = system(cmd);
            if status ~= 0
                error('Head restore failed:\n%s', cmdout);
            else
                fprintf('Head restored OK: %s\n', outHead);
            end
        end

        % restore teeth 
        teethList = dir(fullfile(pathPrep, '*msk_teeth.nii.gz'));
        if isempty(teethList)
            warning('No teeth mask (*msk_teeth.nii.gz) found in: %s', pathPrep);
        else
            [~, idxT] = max([teethList.datenum]);
            teethMask = fullfile(pathPrep, teethList(idxT).name);
            outTeeth  = strrep(teethMask, '.nii.gz', '_fullresol.nii.gz');

            cmd = sprintf(['"%s" run -n %s python "%s" --mask "%s" --json "%s" --out "%s"'], ...
                CONDA, ENVNAME, PY_RECOVER, teethMask, jsonFile, outTeeth);

            fprintf('\n[restore teeth]\nmask: %s\njson: %s\nout : %s\n', teethMask, jsonFile, outTeeth);
            [status, cmdout] = system(cmd);
            if status ~= 0
                error('Teeth restore failed:\n%s', cmdout);
            else
                fprintf('Teeth restored OK: %s\n', outTeeth);
            end
        end

        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define pixel indexes (landmarks) based on large teeth masks
        %%% (upper and lower teeth) for MoCo
        % generate location.txt in the path of the case as:
        % pixel index of top of the head
        % pixel index of the middle of the lips/teeth
        % pixel index of the bottom of the chin/lower teeth
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pathPrep = '...'; % already set previously

        % teethMaskPath should be the full path to ...msk_teeth_fullresol.nii.gz
        f = dir(fullfile(pathPrep, '*msk_teeth_fullresol.nii.gz'));
        if isempty(f)
            error('Cannot find *msk_teeth_fullresol.nii.gz in: %s', pathPrep);
        end
        if numel(f) > 1
            [~, idx] = max([f.datenum]);
            warning('Multiple matches found, using most recent: %s', f(idx).name);
            teethMaskPath = fullfile(pathPrep, f(idx).name);
        else
            teethMaskPath = fullfile(pathPrep, f(1).name);
        end
        
        locPath = fullfile(caseFolder, 'location.txt');

        [idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid, locPath] = ...
            compute_landmarks_from_teeth_mask(teethMaskPath, locPath);

        fprintf('Landmarks saved to: %s\n', locPath);
        fprintf('HF indices: head=%d, lips=%d, chin=%d (RL mid=%d)\n', ...
            idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid);
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Call the MoCo Recon
        % run moco reconstruction with AlignedSENSE (constrained in ROI
        % range in HF) for:
        % - MoCo in fullFOV
        % - MoCo in upper teeth/jaw only
        % - MoCo in lower teeth/jaw only
        %
        % The result will be saved in An-Ve folders
        % Automatically check to avoid repeated recon is ON - see skipRecon
        % related content
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %check whether we have run this before
        skipRecon = 0;
        anVeDir = fullfile(caseFolder, 'An-Ve');
        if ~isfolder(anVeDir)
            warning('Folder not found: %s', anVeDir);
        else
            % seqName: datFiles.name without trailing ".dat"
            seqName = datFiles(fIdx).name;
            if endsWith(seqName, '.dat', 'IgnoreCase', true)
                seqName = extractBefore(seqName, strlength(seqName) - strlength(".dat") + 1);
            else
                warning('datFiles.name does not end with .dat: %s', datFiles.name);
            end

            f1 = fullfile(anVeDir, [seqName, '_Di_MotCorr_lowerjaw_.nii']);
            f2 = fullfile(anVeDir, [seqName, '_Di_MotCorr_upperjaw_.nii']);
            f3 = fullfile(anVeDir, [seqName, '_Di_MotCorr.nii']);

            if isfile(f1) && isfile(f2) && isfile(f3)
                skipRecon = 1;
                fprintf('[PIPE] Found all An-Ve outputs. Skip following recon. \n');
            end
        end
        
        if ~skipRecon
            % ---- Save loop state (deployRecon_dental may clear) ----
            stateFile = 'batchRecon_dental_state.mat';
            currDir   = pwd;

            save(stateFile, ...
                 'rootFolder', 'studiesFile', 'numCases', 'caseList', ...
                 'caseIdx', 'fIdx', 'currDir', ...
                 'datFiles', 'caseFolder', 'seqSelect');

            % ---- Run reconstruction: fullFOV, upper & lower teeth ----

            try
                fprintf('Running deployRecon_dental (MoCo)...\n');
                fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
                deployRecon_dental_MoCo;
                fprintf('MoCo Reconstruction completed.\n');
                fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))

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
        end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Image fusion
        % Fuse the per-sequence upper-jaw and lower-jaw MoCo outputs and
        % save the fused image using the naming/output logic in
        % zihan_tools/image_fusion.m.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        seqName = datFiles(fIdx).name;
        if endsWith(seqName, '.dat', 'IgnoreCase', true)
            seqName = extractBefore(seqName, strlength(seqName) - strlength(".dat") + 1);
        end

        try
            fprintf('Running image fusion...\n');
            fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
            fusedFile = image_fusion(caseFolder, seqName, 'verbose', true);
            fprintf('Image fusion completed.\n');
            fprintf('Fused image saved to: %s\n', fusedFile);
            fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
        catch ME
            fprintf(2, ...
                'Image fusion failed in Case %d, file %s:\n%s\n', ...
                caseIdx, baseName, ME.message);
        end

        fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
%     end % loop over dat files

end % loop over cases

fprintf('Time: %s \n', datestr(now, 'yyyy/mm/dd HH:MM:SS'))
fprintf('\nAll requested cases processed.\n')

fprintf('==== PIPE finished: %s ====\n', datestr(now));
diary off
