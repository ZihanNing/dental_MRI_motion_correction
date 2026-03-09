%% ZN: only used for batch processing

clearvars -except varargin;    % clear function workspace
clear global;                  % clear any globals that may hold pathData etc.
clear persistent;              % reset persistent variables in other functions
clear functions;               % clear all loaded function workspaces from memory

close all;

%% Main computation
addpath(genpath('brackenier-tools'));
addpath(genpath('dhcp-repo-release07'));
addpath(genpath('zihan_tools'));

addpath(genpath('Methods'));
addpath(genpath('Studies-deploy'));

gpuDevice(2); % ZN: for out-of-memory issue happen on gadgetron07-pc

% SELECT STUDIES
%List with all the names of the studies
idPath=1;
studies;
idFile={1};
idRef = ones(1,length(idFile));

% Only consider files selected to reconstruct
[pathIn, fileIn, refIn, refBIn, B0In, B1In,...
 fileUnique, refUnique, refBUnique, B0Unique, B1Unique, ...
 isPTUnique, noiseFileUnique, supportReadoutUnique, resRecUnique, RDesiredUnique, facFOVThUnique] ...
 = extractStudies(pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath, idFile, idRef,[],[], isPT, noiseFile, supportReadout, resRec, RDesired, facFOVTh);

assert( length(idFile) == length(fileIn{1}),'File handling not correct.');

% DATA CONVERSION PARAMETERS
invDataExplicit=0;
estSensExplicit=0;

% DATA PREPARATION PARAMETERS
if ~exist('supportReadoutRecon','var'); supportReadoutRecon = fillCell(refIn{1},[]);end
if ~exist('resRecRecon','var'); resRecRecon = fillCell(refIn{1},[]);end
if ~exist('underS','var'); underS = fillCell(refIn{1},[]);end

%% PROCESS DATA 
p=1;
%BUILD THE ACQUISITION DATA - start with this if an acquisition file is also used as a reference, you don't read it ouy twice (dat2Se checks if rec structure exist whereas dat2Rec reads data out by default)
% [~, ~, facFOVThOrig] = setPilotToneConverion(); % ZN: disabled for unable
flag_useCAIPI = 0; % ZN: if using CAIPI acceleration

% to read one file
for f=1:length(fileUnique{p})
    fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
    %Handling noise file
    if ~isempty(noiseFileUnique{p}{f})
        noiseFileTemp = strcat(pathIn{1},filesep,noiseFileUnique{p}{f});
    elseif ~strcmp(fileIn{p}{f},refIn{1}{1})
        noiseFileTemp=strcat(pathIn{1},filesep,refIn{1}{1});
    else
        noiseFileTemp=[];
    end    
    %Convert
    if ~existsFileVar(fileName)|| invDataExplicit
        if ~flag_useCAIPI % ZN: for ordinary accelration pattern, like GRAPPA, mSENSE etc
            dat2Rec(fileName,supportReadoutUnique{p}{f},[],1,[],isPTUnique{p}{f},[],noiseFileTemp,resRecUnique{p}{f},[],RDesiredUnique{p}{f});       
        else % ZN: for CAIPI acceleration pattern
            removeZeros = 0
            dat2Rec_CAIPI(fileName,supportReadoutUnique{p}{f},[],1,[],isPTUnique{p}{f},[],noiseFileTemp,resRecUnique{p}{f},removeZeros,RDesiredUnique{p}{f});       
        end
    end                
end

%BUILD THE REFERENCE DATA - SENSITIVITY MAPS
useACS_flag = 1; % ZN: 1- use ACS line for coil sensitivity map estimation; 0 - use external ref scan for coil sensitivity map estimation
if useACS_flag == 0 % ZN: use external ref
    for f=1:length(refUnique{p})
        fileRef=strcat(pathIn{p},refUnique{p}{f});
        fileRefB=strcat(pathIn{p},refBUnique{p}{f});
        if existsFileVar(fileRef,{'recS'})~=2  || estSensExplicit
            dat2Se(fileRef,fileRefB,[],1,1,[]);
        end
    end
else % ZN: use ACSline
    fileAcq=strcat(pathIn{p},filesep,fileIn{p}{f});%Is allowed to be a cell array
    if ( ~exist(strcat(fileName,'_ACS.mat'), 'file') || ~ismember('recS',who('-file', strcat(fileName,'_ACS.mat')) ) ) || estSensExplicit && ~strcmp(fileName,'_ACS')
        ACS2rec(char(fileAcq),1,1);
    end
end

%% RUN RECONSTRUCTIONS
f = 1
% for f=1:length(idFile)
    %%PREPARE RECONSTRUCTION STRUCTURE
    %%% PREPARE RECONSTRUCTION STRUCTURE
    fileRef=strcat(pathIn{p},filesep,refIn{p}{f});
    fileB1=strcat(pathIn{p},filesep,B1In{p}{f});
    fileB0=strcat(pathIn{p},filesep,B0In{p}{f});
    fileAcq=strcat(pathIn{p},filesep,fileIn{p}{f});%Is allowed to be a cell array - strcat works with this

    %%% MERGE DATA FROM ALL STRUCTURES
    writeSWCCFlag=0;
%     rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f}, underS{f});
    % modified by ZN % to enable ACS line for coil sensitivity map
    % estimation
    if useACS_flag == 0 % ZN: use external ref for coil sensitivity map estimation
        rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f},underS{f}, useACS_flag);
    else % ZN: use ACS line for coil sensitivity map estimation
        rec=prepareRec(fileAcq, [], fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f},underS{f}, useACS_flag);
    end
    
    %%% ADD SPATIAL INFROMATION ABOUT COILS
    rec.Par.preProcessing.coilGeom.RAS = coilCentroids(rec.S,rec.Par.Mine.APhiRec);
    [~,rec.Par.preProcessing.coilGeom.idxIS] = sort(rec.Par.preProcessing.coilGeom.RAS(3,:));
    ccm = eye(size(rec.S,4));
    ccm = ccm(flip(rec.Par.preProcessing.coilGeom.idxIS),:);
    rec.Par.preProcessing.coilGeom.ccmIS = ccm;ccm=[];

    %%% COMPRESS COILS BEFORE CALLING MOCO
%     percCoil = [17];
%     if ~isempty(percCoil)
%         rec.Par.preProcessing.compressCoils.NChaOrig = size(rec.y,4);
%         rec.Par.preProcessing.compressCoils.percCoil = percCoil;
%         [rec.S,rec.y] = compressCoils(rec.S,percCoil,rec.y);
%         rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction
%         fprintf('Coil compression: Number of coils compressed from %d to %d elements (%d%%).\n',rec.Par.preProcessing.compressCoils.NChaOrig,size(rec.y,4),percCoil*100);
%     end
    
    rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction

    %%% GENERAL RECONSTRUCTION PARAMETERS
    rec.Alg.UseSoftMasking=0;%2 if not masking
    rec.Alg.parXT.writeInter=0;

    rec.Dyn.GPUbS = [6 7];%[2 4]
    rec.Dyn.MaxMem = [6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing

    %% MOTION CORRECTION
    if contains(rec.Names.Name, 'T2', 'IgnoreCase', true)
        contast_scheme = 't2w'; % or 't2w' rec.Names.Name
    elseif contains(rec.Names.Name, 'PD', 'IgnoreCase', true)
        contast_scheme = 'pdw';
    elseif contains(rec.Names.Name, 'MPRAGE', 'IgnoreCase', true)
        contast_scheme = 'mprage';
    else
        disp('There is not recognizable contrast in the seq name!');
        pause;
    end
    
    rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
    %Disabling B0 estimation
    rec.Alg.parXB.useSH = 0;
    rec.Alg.parXB.useTaylor = 0;
    rec.Alg.parXB.useSusc = 0;

    %Motion estimation parameters
    rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
    switch contast_scheme
        case 'pdw'
            rec.Alg.parXT.groupSweeps=12;% Factor to group sweeps
        otherwise
            rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
    end
    
    
    % ==============
    % ZN: Used in batch processing
    % ==============
    loc_file = 'location.txt';
    FOV_range = compute_FOV_range(rec,loc_file);
    rec.Alg.parXT.redFOVflex = 1; % ZN: get flexible redFOV
    

    rec.Alg.parXT.disableGrouping = 1;
    rec.Alg.parXT.convTransformJoint = 1;
    rec.Alg.parXT.UseGiRingi=0;

    %Optimisation scheme (3 levels)
%     nIt = 7;
%     rec.Alg.resPyr =  [ .5  .5  1];
%     rec.Alg.nExtern = [ nIt   nIt   1];
%     rec.Alg.parXT.estT = [ 1 1 0 ];
%     rec.Alg.parXT.PT.subDiv = [0 0 0];
    

    if isequal(contast_scheme,'t2w') || isequal(contast_scheme,'mprage')
        nIt = 15;
        rec.Alg.resPyr =  [ .5  .5  1];
        rec.Alg.nExtern = [ nIt   nIt   1];
        estimateOnFullResol = 1;

        if estimateOnFullResol 
           rec.Alg.parXT.estT = [ 1 1 1 ];
           rec.Alg.nExtern(end)=rec.Alg.nExtern(1);
        else
            rec.Alg.parXT.estT = [ 1 1 0 ];
        end
        rec.Alg.parXT.PT.subDiv = [1 0 0];
    else
        nIt = 15;
        rec.Alg.resPyr =  [ .5  1];
        rec.Alg.nExtern = [ nIt   1];
        estimateOnFullResol = contains(fileIn{p}{f}{1},'midres');
        if estimateOnFullResol 
           rec.Alg.parXT.estT = [ 1 1 ];
           rec.Alg.nExtern(end)=rec.Alg.nExtern(1);
        else
            rec.Alg.parXT.estT = [ 1 0 ];
        end
        rec.Alg.parXT.PT.subDiv = [0 0]; 
    end
    
    
    %Deal with shot duration and sequence dependencies
    if ~rec.Enc.DISORDER.Flag %Need to specify the shot definition as it will not be detected automatically.
        shotDuration = 2;%In seconds - ad hoc
        numShots = [];
        if ~isempty(shotDuration);rec.Alg.parXT.sampleToGroup = shotDuration/(1e-3*rec.Par.Labels.RepetitionTime);end
        if ~isempty(numShots);rec.Alg.parXT.sampleToGroup = round(length(rec.Assign.z{2})/numShots);end
    end

    if mod(length(rec.Assign.z{2})/round(rec.Par.Labels.TFEfactor),1)~=0%Inconsistent data
        rec.Alg.parXT.sampleToGroup=round(length(rec.Assign.z{2})/(rec.Par.Labels.NShots*size(rec.y,5)));
        rec.Par.Labels.TFEfactor=length(rec.Assign.z{2})/size(rec.y,5);%Fake steady state sequence 
    end 

    %%% Set mask
    % read in the prepared mask (head by nnunet seg)
    segDir = fullfile(rec.Names.pathOu, 'An-Aq', 'seg');
    if ~isfolder(segDir)
        error('segDir not found: %s', segDir);
    end
    allFiles = dir(segDir);
    allFiles = allFiles(~[allFiles.isdir]);
    nameKey = rec.Names.Name;  % required substring
    hit = false(size(allFiles));
    for i = 1:numel(allFiles)
        fn = allFiles(i).name;
        isNii = endsWith(fn, '.nii', 'IgnoreCase', true) || endsWith(fn, '.nii.gz', 'IgnoreCase', true);
        hasKeys = contains(fn, nameKey) && contains(fn, 'msk_head_fullresol');
        hit(i) = isNii && hasKeys;
    end
    cand = allFiles(hit);
    if isempty(cand)
        error('No head fullres mask found in %s with name containing "%s" and "msk_head_fullresol" (.nii/.nii.gz).', ...
            segDir, nameKey);
    end

    % If multiple matches, use the most recent one
    [~, idx] = max([cand.datenum]);
    headMaskPath = fullfile(segDir, cand(idx).name);

    fprintf('[PIPE] head fullres mask: %s\n', headMaskPath);

    % Read with correct geometry
    headInfo = niftiinfo(headMaskPath);
    headMask = niftiread(headInfo);
    headMask = uint8(headMask);   % masks as integer labels (safe)

    % Geometry you may want later:
    % - headInfo.Transform.T : voxel->world affine
    % - headInfo.ImageSize   : volume size
    % - headInfo.PixelDimensions : voxel sizes
    headAffine = headInfo.Transform.T;
    fprintf('[PIPE] headMask size: [%s]\n', num2str(size(headMask)));
    
    % expand the mask a bit and not let it to tight (endore some motion)
    BW = headMask > 0;
    N = 30; % expand mask by N voxels % 10 used to be 
    BWd = imdilate(BW, strel('sphere', N));   % 3D dilation
    BWd = imclose(BWd, strel('sphere', 1));  % light closing
    BWd = imfill(BWd, 'holes');             % fill internal holes (3D)
    headMask_exp = uint8(BWd);         % back to uint8 (0/1)
    
    % set the head mask
    if isequal(size(rec.M),size(headMask_exp))
        rec.M = single(headMask_exp);
    else
        error('The readin mask does not matches the original image size.\n');
    end
    
    % set an elliptical mask to mask out the edges anyway
%     rTh=[3 1.2 1.2];
%     ellip_msk = gather(getEllipsMask(rec.ref_y,rTh)==1);
%     % expand the ACS img generated mask to tolerant some movement
%     mask_expanded = expand_mask(rec.M, 5);
%     % generate the final mask
%     rec.M = mask_expanded.* ellip_msk;

%     switch contast_scheme
%         case 'mprage'
%             rTh=[3 1 1]; 
%             rec.M = getEllipsMask(rec.M,rTh)==1;
%         case 't2w' 
%             rTh=[3 1 1]; 
%             rec.M = getEllipsMask(rec.M,rTh)==1;
%         case 'pdw' 
%             rTh=[3 1 1]; % mprage
%             rec.M = getEllipsMask(rec.M,rTh)==1;
%         otherwise
%             flag_nomask = 1;
%             rec.M = ones(size(rec.M));
%     end
    
%     % do not mask
%     flag_nomask = 1;
%     rec.M = ones(size(rec.M));
%     
    
    % final, show the mask
    plotND([],RSOS(rec.ref_y),defRange(RSOS(rec.ref_y)),[],0,[],rec.Par.Mine.APhiRec,[],rec.M,{2});

    isMeGRE = size(rec.y,8)>1;
    isMP2RAGE = size(rec.y,10)>1;
    if isMeGRE
        dimToLoop = 8; suff = '_echo';
    elseif isMP2RAGE
        dimToLoop = 10; suff = '_inv';
    else
        dimToLoop = 8;suff = '';
    end
    

    %% Run (full FOV)
    rec.Alg.parXT.redFOVflex = 0; % ZN: get flexible redFOV
    rec.Alg.parXT.redFOV = 0; % force to be 0 for the ddMRI
    suff = '_fullFOV_';
       
    recTemp = rec;
    idToRun = 1:size(rec.y,8);
    doEveryEchoIndividually = 1;
    
    for e=1:length(idToRun)
        if isMeGRE || isMP2RAGE
            recTemp.Plan.SuffOu=sprintf('%s%d',suff,idToRun(e));
        else
            recTemp.Plan.SuffOu='';
        end
        recTemp.y = dynInd(rec.y,idToRun(e),dimToLoop);
        if isMeGRE
            recTemp.Par.Labels.TE = rec.Par.Labels.TE(idToRun(e));
        elseif isMP2RAGE
            recTemp.Par.Labels.inversionTime = rec.Par.Labels.inversionTime(idToRun(e));
        end
        if 1 % enable for space
            if size(rec.y,6)>1
                recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
                for i=2:3
                    recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
                end
            else
                idToRun = 1:size(rec.y,dimToLoop);
                doEveryEchoIndividually = 0;
            end
            numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
            recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);  
        end
        if e==1 || doEveryEchoIndividually
           temp = solveXTB_PT(recTemp,[]);
           T = temp.T;
           temp = [];
        else
           solveXTB_PT(recTemp,T);
        end
        close all
    end
    
    %% Run (upper-jaw-only)
    % upper-jaw only
    rec.Alg.parXT.redFOVflex = 1;
    
    % ZN: use the info from location.txt
    if rec.Alg.parXT.redFOVflex
        if contains(rec.Par.Scan.Mine.RO,'H')
            % general teeth
            rec.Alg.parXT.redFOV_upper=1-FOV_range(1);
            rec.Alg.parXT.redFOV_lower=FOV_range(2);
        else
            rec.Alg.parXT.redFOV=0;
        end
    else
%         if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
        rec.Alg.parXT.redFOV = 0; % force to be 0 for the ddMRI
    end
    suff = '_upperjaw_';
       
    recTemp = rec;
    idToRun = 1:size(rec.y,8);
    doEveryEchoIndividually = 1;
    usePriorMotionStates = 1; % ZN: only for dental MRI: by activating this (=1), 
                              % use motion states estimated by fullFOV as the initial states of the rec.T
                              % so that the bulk head movement could be
                              % better captured
    if usePriorMotionStates % ZN: when set as this mode, first check whether the fullFOV result has been computed
        targetDir = [rec.Names.pathOu,'An-Ve/'];
        targetStr = rec.Names.Name;

        files = dir(fullfile(targetDir, '*.mat'));

        matchList = {};

        for k = 1:numel(files)
            fname = files(k).name;
            if contains(fname, '_Tr_MotCorr') && contains(fname, targetStr)
                matchList{end+1} = fullfile(targetDir, fname); 
            end
        end

        if isempty(matchList)
            fprintf('No matching fullFOV results found - initial motion states set all zeros!\n');
            usePriorMotionStates = 0;
        else
            fprintf('Found %d matching fullFOV results files:\n', numel(matchList));
            disp(matchList(:));
            try
                rec.init_T = load(matchList{1});
                rec.init_T.useflag = 1;
                fprintf('fullFOV results loaded and about to be used as the inital motion states: %s \n',matchList{1});
            catch ME
                usePriorMotionStates = 0;
                fprintf(2, 'Error in loading the fullFOV results: %s \n',matchList{1});
            end
        end
    end
    
    for e=1:length(idToRun)
        if isMeGRE || isMP2RAGE
            recTemp.Plan.SuffOu=sprintf('%s%d',suff,idToRun(e));
        else
            recTemp.Plan.SuffOu=sprintf(suff);
        end
        recTemp.y = dynInd(rec.y,idToRun(e),dimToLoop);
        if isMeGRE
            recTemp.Par.Labels.TE = rec.Par.Labels.TE(idToRun(e));
        elseif isMP2RAGE
            recTemp.Par.Labels.inversionTime = rec.Par.Labels.inversionTime(idToRun(e));
        end
        if 1 % enable for space
            if size(rec.y,6)>1
                recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
                for i=2:3
                    recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
                end
            else
                idToRun = 1:size(rec.y,dimToLoop);
                doEveryEchoIndividually = 0;
            end
            numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
            recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);  
        end
        if e==1 || doEveryEchoIndividually
           temp = solveXTB_PT(recTemp,[]);
           T = temp.T;
           temp = [];
        else
           solveXTB_PT(recTemp,T);
        end
        close all
    end

    %% Run (lower-jaw-only)
    % lower-jaw only
    rec.Alg.parXT.redFOVflex = 1;
    if rec.Alg.parXT.redFOVflex
        if contains(rec.Par.Scan.Mine.RO,'H')
            % general teetch
            rec.Alg.parXT.redFOV_upper=1-FOV_range(2);
            rec.Alg.parXT.redFOV_lower=FOV_range(3);
        else
            rec.Alg.parXT.redFOV=0;
        end
    else
%         if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
        rec.Alg.parXT.redFOV = 0; % force to be 0 for the ddMRI
    end
    suff = '_lowerjaw_';
       
    recTemp = rec;
    idToRun = 1:size(rec.y,8);
    doEveryEchoIndividually = 1;
    
    for e=1:length(idToRun)
        if isMeGRE || isMP2RAGE
            recTemp.Plan.SuffOu=sprintf('%s%d',suff,idToRun(e));
        else
            recTemp.Plan.SuffOu=sprintf(suff);
        end
        recTemp.y = dynInd(rec.y,idToRun(e),dimToLoop);
        if isMeGRE
            recTemp.Par.Labels.TE = rec.Par.Labels.TE(idToRun(e));
        elseif isMP2RAGE
            recTemp.Par.Labels.inversionTime = rec.Par.Labels.inversionTime(idToRun(e));
        end
        if 1 % enable for space
            if size(rec.y,6)>1
                recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
                for i=2:3
                    recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
                end
            else
                idToRun = 1:size(rec.y,dimToLoop);
                doEveryEchoIndividually = 0;
            end
            numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
            recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);  
        end
        if e==1 || doEveryEchoIndividually
           temp = solveXTB_PT(recTemp,[]);
           T = temp.T;
           temp = [];
        else
           solveXTB_PT(recTemp,T);
        end
        close all
    end

