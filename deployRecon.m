

clear all;
close all;

addpath(genpath('brackenier-tools'));
addpath(genpath('dhcp-repo-release07'));

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
    rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
    %Disabling B0 estimation
    rec.Alg.parXB.useSH = 0;
    rec.Alg.parXB.useTaylor = 0;
    rec.Alg.parXB.useSusc = 0;

    %Motion estimation parameters
    rec.Alg.AlignedRec=2;%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
    rec.Alg.parXT.groupSweeps=1;% Factor to group sweeps
    rec.Alg.parXT.redFOVflex = 1; % ZN: get flexible redFOV
    if rec.Alg.parXT.redFOVflex
        if contains(rec.Par.Scan.Mine.RO,'H')
%             % upper teeth
%             rec.Alg.parXT.redFOV_upper=0.4;
%             rec.Alg.parXT.redFOV_lower=0.4;
%             % lower teeth
%             rec.Alg.parXT.redFOV_upper=0.6;
%             rec.Alg.parXT.redFOV_lower=0.2;
            % general teetch
            rec.Alg.parXT.redFOV_upper=0.32;
            rec.Alg.parXT.redFOV_lower=0.27;
        else
            rec.Alg.parXT.redFOV=0;
        end
    else
        if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
%         rec.Alg.parXT.redFOV = 0; % force to be 0 for the ddMRI
    end

    rec.Alg.parXT.disableGrouping = 1;
    rec.Alg.parXT.convTransformJoint = 1;
    rec.Alg.parXT.UseGiRingi=0;

    %Optimisation scheme (3 levels)
%     nIt = 7;
%     rec.Alg.resPyr =  [ .5  .5  1];
%     rec.Alg.nExtern = [ nIt   nIt   1];
%     rec.Alg.parXT.estT = [ 1 1 0 ];
%     rec.Alg.parXT.PT.subDiv = [0 0 0];
    
%     nIt = 7;
%     rec.Alg.resPyr =  [ .5  .5  1];
%     rec.Alg.nExtern = [ nIt   nIt   1];
%     estimateOnFullResol = contains(fileIn{p}{f}{1},'midres');
% 
%     if estimateOnFullResol 
%        rec.Alg.parXT.estT = [ 1 1 1 ];
%        rec.Alg.nExtern(end)=rec.Alg.nExtern(1);
%     else
%         rec.Alg.parXT.estT = [ 1 1 0 ];
%     end
%     rec.Alg.parXT.PT.subDiv = [1 0 0];
%     
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

    %Set mask
%     rTh=[.9 .8]; % mprage
    rTh=[1 .9]; % space
    rec.M = getEllipsMask(rec.M,rTh)==1;
    
    % do not mask
%     flag_nomask = 1;
%     rec.M = ones(size(rec.M));
%     
%     % set ZNmask
%     rec.Alg.UseSoftMasking=0;% 3 if use ZN mask
%     S_tmp = RSOS(rec.ref_y);
%     S_tmp = permute(S_tmp,[3 2 1]); % masking in the sagittal direction
%     M_zn = simple_mask(S_tmp,0.01,800,10); % thres = 0.3; isolated_region_size = 800; enlarge_factor = 15
%     % Use the processed binary mask as the rec.M
%     M_zn = permute(M_zn,[3 2 1]); % masking in the sagittal direction
%     rec.M = M_zn;
%     S_tmp = []; M_zn = [];
    
    % mask for dental
%     rec.Alg.UseSoftMasking=0;% 3 if use ZN mask
%     rec.M(1:30,:,:) = 0;
%     rec.M(230:256,:,:) = 0; 
% %     rec.M(1:60,:,:) = 0;
% %     rec.M(190:256,:,:) = 0; 
%     S_tmp = []; M_zn = [];
    
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
    
    %% Run (Original)
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
        if 0 % enable for space
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
    

%      %% Run (multi-NSA or single-NSA CAIPI)
%     rec.Alg.disabledisplay = 1; % 1 - not display the image during the recon % ZN
% %     rec.Enc.UnderSampling.R = [2 2]; % ZN: manually set the acceleration factor temporally 
%     recTemp = rec;
%     recTemp.Plan.Suff='_MotCorr';recTemp.Plan.SuffOu='_outl';
%     if size(rec.y,6)>1
%         recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
%         for i=2:3
%             recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
%         end
%     else
%         recTemp = rec;
%         idToRun = 1:size(rec.y,dimToLoop);
%         doEveryEchoIndividually = 0;
%     end
%     numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
% %     numShots = size(recTemp.y,5)*rec.Enc.DISORDERInfo.NShots;
%     recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);          
%     temp = solveXTB_PT_FlexSamp(recTemp);
% 
%     % Uncorrected
%     recTemp.Plan.Suff='_UnCorr';recTemp.Plan.SuffOu='_outl';
%     solveXTB_PT_FlexSamp(recTemp, zerosL(temp.T));
% 
%     %% MOTION + B0 CORRECTION (BASIS FUNCTIONS)
%     rec.Plan.Suff='_MotdB0Corr';
% 
%     rec.Alg.parXB.useSH = 1;%0 if disabled, 1 if in head frame, 2 if in scanner frame
%     rec.Alg.parXB.SHorder = 2; % disable by setting < 0
%     rec.Plan.SuffOu= sprintf('_SH%d', rec.Alg.parXB.SHorder);
%     rec.Alg.parXB.deTaylor = [2 inf];%(0) de-activated / (1) move cr to D / (2) make D estimate with SH and discard cr / (3) if estimated 1 SH term and no motion dependence
%     rec.Alg.parXB.Optim.basisDelay=1;  
% 
%     rec.Alg.parXT.estB = [ 1 1 0 ];
% 
%     %Run
%     recTemp = rec;
%     for e=1:length(idToRun)
%         if isMeGRE || isMP2RAGE
%             recTemp.Plan.SuffOu=sprintf('%s%d',suff,idToRun(e));
%         else
%             recTemp.Plan.SuffOu='';
%         end
%         recTemp.y = dynInd(rec.y,idToRun(e),dimToLoop);
%         if isMeGRE
%             recTemp.Par.Labels.TE = rec.Par.Labels.TE(idToRun(e));
%         elseif isMP2RAGE
%             recTemp.Par.Labels.inversionTime = rec.Par.Labels.inversionTime(idToRun(e));
%         end
%         if 1 % enable for space
%             numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
%             recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);  
%         end
%         if e==1 || doEveryEchoIndividually
%            temp = solveXTB_PT(recTemp,[]);
%            T = temp.T;
%            c = temp.Db.c;
%            temp = [];
%         else
%            solveXTB_PT(recTemp,T,[],c);
%         end
%         close all
%     end

% end

