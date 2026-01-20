
clear all;
addpath(genpath('/home/ybr19/Projects/Simulations/Reconstruction/Methods/synthesizeKSpace'))
addpath(genpath('/home/ybr19/Projects/B0-Shimming/Functions'));
addpath(genpath('/home/ybr19/Projects/Reconstruction'));
addpath(genpath('/home/ybr19/Projects/PilotTone'));
addpath(genpath('/home/ybr19/Software/Utilities'));
addpath(genpath('/home/ybr19/Software/DISORDER/DefinitiveImplementationRelease07'));

addpath(genpath('/home/ybr19/Projects/MP2RAGE/Methods'))

for studyId = [7]
    
    idFile = [];

    clear supportReadoutRecon
    clear resRecRecon

    %% SELECT STUDIES
    %List with all the names of the studies
    idPath=1;
    if studyId == 1
        studies_20240201_SWI;
        if isempty(idFile);idFile={3};end
    elseif studyId == 2
        studies_20240201_SPACE_CAIPI;
        if isempty(idFile);idFile={1};end
    elseif studyId==3
        studies_20240207_SPACE_CAIPI;
        if isempty(idFile);idFile={2};end
    elseif studyId==4
        studies_20240222_CAIPI_phantom;
        if isempty(idFile);idFile={1,2};end
    elseif studyId==5
        studies_20240222_CAIPI_invivo;
        if isempty(idFile);idFile={2};end
    elseif studyId==6
        studies_20240228_CAIPI_invivo;
        if isempty(idFile);idFile={[1,2]};end
    elseif studyId==7
        studies_20240318_FLAIR;
        if isempty(idFile);idFile={1,2,3};end
    end
    idRef = 1*ones(1,length(idFile));

    % Only consider files selected to reconstruct
    [pathIn, fileIn, refIn, refBIn, B0In, B1In,...
     fileUnique, refUnique, refBUnique, B0Unique, B1Unique, ...
     isPTUnique, noiseFileUnique, supportReadoutUnique, resRecUnique, RDesiredUnique, facFOVThUnique] ...
     = extractStudies(pathIn, fileIn, refIn, refBIn, B0In, B1In, idPath, idFile, idRef,[],[], isPT, noiseFile, supportReadout, resRec, RDesired, facFOVTh);

    assert( length(idFile) == length(fileIn{1}),'File handling not correct.');

    %% DATA CONVERSION PARAMETERS
    invDataExplicit=0;
    estSensExplicit=0;
    estB1Explicit=0;
    estB0Explicit=0;

    %% DATA PREPARATION PARAMETERS
    if ~exist('supportReadoutRecon','var'); supportReadoutRecon = fillCell(refIn{1},[]);end
    if ~exist('resRecRecon','var'); resRecRecon = fillCell(refIn{1},[]);end
    if ~exist('underS','var'); underS = fillCell(refIn{1},[]);end
    
    %% PROCESS DATA 
    p=1;
    
    %BUILD THE REFERENCE DATA - SENSITIVITY MAPS
    for f=1:length(refUnique{p})
        fileRef=strcat(pathIn{p},filesep,refUnique{p}{f});
        fileRefB=strcat(pathIn{p},filesep,refBUnique{p}{f});
        if ~existsFileVar(fileRef)|| invDataExplicit
            dat2Rec_CAIPI(fileRef,[],[],1,[],[],[],[],[],[],[],0);
        end
        if existsFileVar(fileRef,{'recS'})~=2  || estSensExplicit
            dat2Se(fileRef,fileRefB,[],1,1,[]);
        end
    end
    
    %BUILD THE ACQUISITION DATA - start with this if an acquisition file is also used as a reference, you don't read it ouy twice (dat2Se checks if rec structure exist whereas dat2Rec reads data out by default)
    [~, ~, facFOVThOrig] = setPilotToneConverion();
    for f=1:length(fileUnique{p})
        fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
        %Handling noise file
        if ~isempty(noiseFileUnique{p}{f})
            noiseFileTemp = strcat(pathIn{1},filesep,noiseFileUnique{p}{f});
        elseif ~strcmp(fileUnique{p}{f},refIn{1}{1})
            noiseFileTemp=strcat(pathIn{1},filesep,refIn{1}{1});
        else
            noiseFileTemp=[];
        end    
        %Handling PT extraction parameters
        setPilotToneConverion([], [], facFOVThUnique{p}{f});
        %Convert
        if ~existsFileVar(fileName)|| invDataExplicit
            removeZeros = 0
            dat2Rec_CAIPI(fileName,supportReadoutUnique{p}{f},[],1,[],isPTUnique{p}{f},[],noiseFileTemp,resRecUnique{p}{f},removeZeros,RDesiredUnique{p}{f});       
        end                
    end
    setPilotToneConverion([], [], facFOVThOrig);


    %% RUN RECONSTRUCTIONS
    for f=1:length(idFile)
        %% PREPARE RECONSTRUCTION STRUCTURE
        %%% PREPARE RECONSTRUCTION STRUCTURE
        fileRef=strcat(pathIn{p},filesep,refIn{p}{f});
        fileB1=strcat(pathIn{p},filesep,B1In{p}{f});
        fileB0=strcat(pathIn{p},filesep,B0In{p}{f});
        fileAcq=strcat(pathIn{p},filesep,fileIn{p}{f});%Is allowed to be a cell array - strcat works with this

        %%% MERGE DATA FROM ALL STRUCTURES
        writeSWCCFlag=0;
        rec=prepareRec(fileAcq, fileRef, fileB1, fileB0, writeSWCCFlag, [], resRecRecon{f}, supportReadoutRecon{f}, underS{f});

        %%% ADD SPATIAL INFROMATION ABOUT COILS
        rec.Par.preProcessing.coilGeom.RAS = coilCentroids(rec.S,rec.Par.Mine.APhiRec);
        [~,rec.Par.preProcessing.coilGeom.idxIS] = sort(rec.Par.preProcessing.coilGeom.RAS(3,:));
        ccm = eye(size(rec.S,4));
        ccm = ccm(flip(rec.Par.preProcessing.coilGeom.idxIS),:);
        rec.Par.preProcessing.coilGeom.ccmIS = ccm;ccm=[];

        %%% COMPRESS COILS BEFORE CALLING MOCO
        percCoil = [17];
        if ~isempty(percCoil)
            rec.Par.preProcessing.compressCoils.NChaOrig = size(rec.y,4);
            rec.Par.preProcessing.compressCoils.percCoil = percCoil;
            [rec.S,rec.y] = compressCoils(rec.S,percCoil,rec.y);
            rec.Alg.parXT.perc = size(rec.y,4)*[1 1 1];%Don't allow for further coil reduction
            fprintf('Coil compression: Number of coils compressed from %d to %d elements (%d%%).\n',rec.Par.preProcessing.compressCoils.NChaOrig,size(rec.y,4),percCoil*100);
        end
        
        %%% GENERAL RECONSTRUCTION PARAMETERS
        rec.Alg.WriteSnapshots=1;       
        rec.Alg.UseSoftMasking=0;%2 if not masking

        rec.Alg.parXT.writeInter=0

        rec.Dyn.GPUbS = [6 7];%[2 4]
        rec.Dyn.MaxMem = [6e6 10e6 1e6];%[6e6 2e6 1e6];%Maximum memory allowed in the gpu: first component, preprocessing, second component, actual CG-SENSE, third component certain elements of preprocessing
        rec.Dyn.BlockGPU = 0;

        %% MOTION CORRECTION
        rec.Plan.Suff='_MotCorr';rec.Plan.SuffOu='';
        %Disabling B0 estimation
        rec.Alg.parXB.useSH = 0;
        rec.Alg.parXB.useTaylor = 0;
        rec.Alg.parXB.useSusc = 0;

        %Motion estimation parameters
        rec.Alg.AlignedRec=2%Type of reconstruction, 1-> no outlier rejection / 2-> outlier rejection / 3-> intra-shot corrections based on outlier detection / 4-> full intra-shot corrections
        rec.Alg.parXT.groupSweeps=2% Factor to group sweeps
        if contains(rec.Par.Scan.Mine.RO,'H');rec.Alg.parXT.redFOV=1/3;else;rec.Alg.parXT.redFOV=0;end
        rec.Alg.parXT.disableGrouping = 1;
        rec.Alg.parXT.convTransformJoint = 1;
        rec.Alg.parXT.traLimXT=[0.03 0.01];%[0.05 0.02]; Respectively translation and rotation limits for XT estimation
        rec.Alg.parXT.meanT=0;
        rec.Alg.parXT.UseGiRingi=0;
        rec.Alg.parXT.fractionOrder=0;%.5;
        rec.Alg.parXT.apod = [.1 0 0 ];

        %Optimisation scheme
        nIt = 7;
        rec.Alg.resPyr =  [ .5  .5  1];
        rec.Alg.nExtern = [ nIt   nIt   1];
        estimateOnFullResol = contains(fileIn{p}{f}{1},'midres');

        if estimateOnFullResol 
           rec.Alg.parXT.estT = [ 1 1 1 ];
           rec.Alg.nExtern(end)=rec.Alg.nExtern(1);
        else
            rec.Alg.parXT.estT = [ 1 1 0 ];
        end
        rec.Alg.parXT.PT.subDiv = [1 0 0];

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
        m = rec.M;
        %rec.M = refineMaskExt(m,[],[],10);
        rTh=1;
        rec.M = getEllipsMask(rec.M,rTh)==1;
        plotND([],RSOS(rec.y,4),defRange(RSOS(rec.y)),[],0,[],rec.Par.Mine.APhiRec,[],rec.M,{2});

        rec.Enc.UnderSampling.Flag = 0
        rec.Alg.AlignedRec=1
        rec.Alg.parXT.maximumDynamics=5;
        
        %% Run
        recTemp = rec;
        recTemp.Plan.Suff='_MotCorr';recTemp.Plan.SuffOu='';
        if size(rec.y,6)>1
            recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
            for i=2:3
                recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
            end
        end
        numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
        recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);          
        temp = solveXTB_PT_MPRAGE(recTemp);
        
        %% Uncorrected
        recTemp.Plan.Suff='_UnCorr';recTemp.Plan.SuffOu='';
        solveXTB_PT_MPRAGE(recTemp, zerosL(temp.T));
        
        rec.Alg.AlignedRec=2
        %% Run
        recTemp = rec;
        recTemp.Plan.Suff='_MotCorr';recTemp.Plan.SuffOu='_outl';
        if size(rec.y,6)>1
            recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
            for i=2:3
                recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
            end
        end
        numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
        recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots);          
        temp = solveXTB_PT_MPRAGE(recTemp);
        
        %% Uncorrected
        recTemp.Plan.Suff='_UnCorr';recTemp.Plan.SuffOu='_outl';
        solveXTB_PT_MPRAGE(recTemp, zerosL(temp.T));
        
        %% extra 
%         recTemp = rec;
%         recTemp.Plan.SuffOu=sprintf('_ave%d',3);
%         recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
%         %recTemp.y = dynInd(rec.y,1,6);
%         for i=2:3
%             recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
%         end
%         numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
%         recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots); 
%         temp = solveXTB_PT_MPRAGE(recTemp);
%         
%         %% alpply
%         dd = makeBins(1:size(temp.T,5),4);
%         T = (dynInd(temp.T,dd==2,5)+dynInd(temp.T,dd==4,5))/2;
%         recTemp = rec;
%         recTemp.Plan.SuffOu=sprintf('_ave%d',3);
%         %recTemp.y = cat(5,dynInd(rec.y,1,6),dynInd(rec.y,2,6)); 
%         recTemp.y = dynInd(rec.y,1,6);
%         %for i=2:3
%         %    recTemp.Assign.z{i} = cat(2,rec.Assign.z{i},rec.Assign.z{i});
%         %end
%         recTemp.y = dynInd(recTemp.y,2,5);
%         for i=2:3
%             recTemp.Assign.z{i} = dynInd(recTemp.Assign.z{i},(length(recTemp.Assign.z{i})/2+1):length(recTemp.Assign.z{i}),2);
%         end
%         numShots = size(recTemp.y,5)*rec.Enc.DISORDER.NShots;
%         recTemp.Alg.parXT.sampleToGroup = round(length(recTemp.Assign.z{2})/numShots); 
%         
%         recTemp.Plan.Suff='_Apply';
%         solveXTB_PT_MPRAGE(recTemp,T);
%         
%         recTemp.Plan.Suff='_UnCorr';
%         solveXTB_PT_MPRAGE(recTemp,zerosL(T));
%         
%         recTemp.Plan.Suff='_MotCorr';
%         solveXTB_PT_MPRAGE(recTemp);                
%         
%         %recTemp.Plan.Suff='_UnCorr';
%         %solveXTB_PT_MPRAGE(recTemp,zerosL(temp.T));

    end
end

