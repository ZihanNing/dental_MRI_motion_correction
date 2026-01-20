
function rec = prepareRec(fileAcq, fileRef, fileB1, fileB0, writeNIIFlag, pathOu, resRec, supportReadout, underS, useACS_flag)

%PREPAREREC prepares a reconstruction structure that can be called by a reconstruction algorithm. It combines acquistition data with the corresponding reference data (Sensitivities/B0/B1). 
%   REC = PREPAREREC(FILEACQ,FILEREF,{FILEB1},{FILEB0},{WRITENIIFLAG},{PATHOU},{RESREC})
%   * FILEACQ is the name of the raw data (.dat) of the acquisition.
%   * FILEREF is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * {FILEB1} is the name of the raw data (.dat) of the reference acquisition used for B1 mapping
%   * {FILEB0} is the name of the raw data (.dat) of the reference acquisition used for B0 mapping
%   * {WRITENIIFLAG} is a flag to write a simple Sensitivity Weighted Coil Combined (SWCC) recon to a NIFTI files. Defaults to 1.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the files are stored.
%   * {RESREC} is the resolution at which to downsample the acquisition.
%   ** REC is the reconstruction structure with all the data arranged.
%   * {useACS_flag} use ACS line (1) or external ref (0) for coil
%   sensitivity map estimation % ZN
%   
%   Yannick Brackenier 2023-07-19

% if nargin<2 || isempty(fileRef);error('prepareRec:: Sensitivities must be provided.'); end
if (nargin<2 || isempty(fileRef)) && useACS_flag == 0; error('prepareRec:: Sensitivities must be provided.'); end % ZN: when there's no external ref provided, ACS_flag need to be set as 1 for using ACS line for coil senstivity map estimation
if ~(nargin<3) && ~isempty(fileB0);[~,fileb0]=fileparts(fileB0);else; fileb0=''; end
if ~(nargin<4) && ~isempty(fileB1);[~,fileb1]=fileparts(fileB1);else; fileb1=''; end

fprintf('<strong>Creating reconstruction structure:</strong>\n');
if ~iscell(fileAcq)%Single acquisition
    [pathOuTemp,file]=fileparts(fileAcq);
    fprintf('  Acquisition:     %s\n',file);
else%Multiple acquisitions
    for i=1:length(fileAcq)
        [pathOuTemp,file]=fileparts(fileAcq{i});
        fprintf('  Acquisition:   %s\n',file);
    end
end

assert(~iscell(fileRef),'prepareRec:: The reference file fileRef cannot contain multiple scans. It must be a singe filename.')
if useACS_flag == 0
    [~,fileref]=fileparts(fileRef);
else % ZN: use ACS line for coil sensitivity map estimation
    fileref = strcat(char(fileAcq),'_ACS'); % ZN: load coil sensitivity map estimated by ACS line
end
fprintf('  Sensitivities:   %s\n',fileref);
if ~isempty(fileb0);fprintf('  B0 map:   %s\n',fileb0);end
if ~isempty(fileb1);fprintf('  B1 map:   %s\n',fileb1);end

if (nargin<2 || isempty(fileRef)) && (useACS_flag == 0); error('matchRec:: Reference data needs to be provided.');end % ZN added - to enable ACS line
if nargin<3 || isempty(fileB1); includeB1 =0;else includeB1 =1;end
if nargin<4 || isempty(fileB0); includeB0 =0;else includeB0 =1;end
% if nargin<2 || isempty(fileRef); error('matchRec:: Reference data needs to be provided.');end
if nargin<5 || isempty(writeNIIFlag); writeNIIFlag=1;end
if nargin<6 || isempty(pathOu); pathOu=[];end
if nargin<7 || isempty(resRec); resRec=[];end
if nargin<8 || isempty(supportReadout); supportReadout=[];end
if nargin<9 || isempty(underS); underS=[];end
if nargin<10; if ~isempty(fileRef);useACS_flag = 0;else; error('No fileRef provided and ACS flag is absent'); end; end % ZN: when useACS_flag is not provided, fileRef must not be empty

if isempty(pathOu); pathOu = pathOuTemp;end

%%% ASSIGN DATA  
if ~iscell(fileAcq)%Single acquisition
    %Load/convert data
    if existsFileVar(fileAcq,{'rec'})~=2; dat2Rec(fileAcq);end
    rec = load(strcat(fileAcq,'.mat'));rec=rec.rec;
    resOrig = rec.Enc.AcqVoxelSize;
    NOrig = size(rec.y,1);
    %Resample if needed
    if ~isempty(resRec) || ~isempty(supportReadout) || ~isempty(underS); rec=resampleRec(rec,resRec,supportReadout,underS);end
else%Combine multiple scans
    for n=1:length(fileAcq)%Run over acquistitions
        %Load/convert data
        if existsFileVar(fileAcq{n},{'rec'})~=2; dat2Rec(fileAcq{n});end
        rec=load(strcat(fileAcq{n},'.mat'));rec=rec.rec;
        resOrig = rec.Enc.AcqVoxelSize;
        NOrig = size(rec.y,1);
        %Resample if needed
        if ~isempty(resRec) || ~isempty(supportReadout) || ~isempty(underS); rec=resampleRec(rec,resRec,supportReadout,underS);end
        %Assign data
        if n==1
           y=rec.y;
           z{2}=rec.Assign.z{2};
           z{3}=rec.Assign.z{3};
           Name=rec.Names.Name;
           if isfield(rec,'PT')
               pSliceImage=rec.PT.pSliceImage;
               pTimeTest=rec.PT.pTimeTest;
           end
        else
           y=cat(5,y,rec.y);
           z{2}=cat(2,z{2},rec.Assign.z{2});
           z{3}=cat(2,z{3},rec.Assign.z{3});
           Name=strcat(Name,rec.Names.Name);
           if isfield(rec,'PT')
               pSliceImage=cat(5,pSliceImage,rec.PT.pSliceImage);
               pTimeTest=cat(2,pTimeTest,rec.PT.pTimeTest);
           end
        end
    end
    rec.y=y;y=[];
    rec.Assign.z{2}=z{2};z{2}=[];
    rec.Assign.z{3}=z{3};z{3}=[];
    rec.Names.Name=Name;
    if isfield(rec,'PT')
        rec.PT.pSliceImage=pSliceImage;
        rec.PT.pTimeTest=pTimeTest;
    end
end
rec.Names.pathOu = pathOu;

%%% LOG THE PRE-PROCESSING STEPS
if ~isempty(resRec)
    rec.Par.preProcessing.resRecRecon.resOrig = resOrig;
    rec.Par.preProcessing.resRecRecon.resNew = resRec;
end
if ~isempty(supportReadout)
    rec.Par.preProcessing.supportReadoutRecon.suppFOV = supportReadout;
    rec.Par.preProcessing.supportReadoutRecon.NOrig = NOrig;
    rec.Par.preProcessing.supportReadoutRecon.NNew = size(rec.y,1);
end
if iscell(fileAcq); rec.Par.preProcessing.fileAcq = fileAcq;end
rec.Par.preProcessing.fileref = fileref;
rec.Par.preProcessing.fileb0 = fileb0;
rec.Par.preProcessing.fileb1 = fileb1;

%%% LOAD COIL SENSITIVITIES
% if existsFileVar(fileRef,{'recS'})~=2; dat2Se(fileRef);end %Check if sensitivities have been estimated
% load(strcat(fileRef,'.mat'), 'recS');
% modified by ZN to enable ACS line
if useACS_flag == 0 % ZN: use external ref for coil sensitivity map estimation
    if existsFileVar(fileRef,{'recS'})~=2; dat2Se(fileRef);end%Check if data has been parsed
    load(strcat(fileRef,'.mat'), 'recS'); 
else % ZN: use ACS line for coil sensitivity map estimation
    fileRef = strcat(char(fileAcq),'_ACS'); % ZN: load coil sensitivity map estimated by ACS line
    if ~existsFileVar(fileRef); error('There is no coil sensitivity map estimated by ACS line'); end
    load(strcat(fileRef,'.mat'), 'recS');if ~exist('recS','var'); error('There is no coil sensitivitity map estimated by ACS line');end
end

%%% COIL ALLOCATION
%Find corresponding coils 
% hdrRef = recS.Varia.hdr;
% hdrAcq = rec.Varia.hdr;
% [idRef,idAcq] = assignCoils(hdrRef, hdrAcq);
% %Extract the relevant coil elements
% recS.S = dynInd(recS.S, idRef,4);
% rec.y = dynInd(rec.y, idAcq,4);
%Compute and apply de-correlation matrix
if isfield(recS,'N') && ~isempty(recS.N)
    recS.N = dynInd(recS.N , idRef,4);
    %rec.N = dynInd(rec.N, idAcq,4);
    [~,~,ccm] = standardizeCoils(ones([1 1 1 size(recS.S,4)],'single'), recS.N);
    recS.S = aplGPU(ccm.',recS.S ,4);
    rec.y = aplGPU(ccm.',rec.y ,4);
end

%%% ASSIGN SENSITIVITIES AND MASK
rec.S=recS.S;
%rec.S=recS.SBart;fprintf('Taking Bart espirit')

%%% SET MASKS
% refimg = permute(abs(recS.x),[2 3 1]); % ZN: assume the HF for RO, here turn imag to axial view for masking
% if contains(recS.Names.Name,'MPRAGE') || contains(recS.Names.Name,'PDw')
%     mask_varagin.p_low = 15; mask_varagin.p_high = 65; mask_varagin.thr_tuning = 0.8;
% else
%     mask_varagin.p_low = 15; mask_varagin.p_high = 55; mask_varagin.thr_tuning = 0.8;
% end
% rec.M = single(dental_background_mask(refimg,mask_varagin));
% rec.M = permute(rec.M,[3 1 2]); % ZN: back to original view
rec.M=recS.x;

NS=multDimSize(recS.S,1:3);
NX=rec.Enc.FOVSize(1:3);

%%% FILTER SENSITIVITIES AND MASK
filterRefData = 0
if filterRefData
    gibbsRing=.5;
    %Create filter
    sp = .5;%was .5
    HS=buildFilter(2*NS,'tukeyIso',sp,0,gibbsRing,1);
    sp = 1;
    HM=buildFilter(2*NS,'tukeyIso',sp,0,gibbsRing,1);
    %Filter
    rec.S=filtering(rec.S,HS,1);
    rec.M=abs(filtering(abs(rec.M),HM,1));        
end

%%% MAP FOV OF ACQUISITION DATA AND REFERENCE DATA
MTS = recS.Par.Mine.APhiRec;MTy = rec.Par.Mine.APhiRec;
% rec.S = mapVolume (rec.S, ones(NX(1:3)), MTS, MTy,[],[],'spline');
% rec.M = mapVolume (rec.M, ones(NX(1:3)), MTS, MTy,[],[],'spline');
% rec.ref_y = mapVolume(RSOS(recS.x),ones(NX(1:3)),MTS,MTy,[],[],'spline'); % by ZN, for masking
% ZN: debug - try geometry map with 'linear'
rec.S = mapVolume (rec.S, ones(NX(1:3)), MTS, MTy,[],[],'linear');
rec.M = mapVolume (rec.M, ones(NX(1:3)), MTS, MTy,[],[],'linear');
rec.ref_y = mapVolume(RSOS(recS.x),ones(NX(1:3)),MTS,MTy,[],[],'linear'); % by ZN, for masking

%%% GENERATE MASK
% refimg = permute(abs(rec.ref_y),[2 3 1]);  % in axial view
% refNX = size(refimg);
% refimg = imgaussfilt3(refimg,1.5);
% refimg = imresize3(refimg,[64 64 64],'linear');
% if contains(rec.Names.Name,'MPRAGE')
%     mask = levelsetseg(refimg, 'preset','MPRAGE');
%     mask = dental_mask_polish(mask); % fill holes
% elseif contains(rec.Names.Name,'PDwSPACE')
%     mask_varagin.p_low = 15; mask_varagin.p_high = 65; mask_varagin.thr_tuning = 0.8;
%     mask = dental_mask_edge_aware_bg_strict(refimg);
% else
%     mask_varagin.p_low = 15; mask_varagin.p_high = 55; mask_varagin.thr_tuning = 0.8;
%     mask = dental_mask_edge_aware_bg_strict(refimg);
% %     mask = levelsetseg(refimg, 'preset','T2wSPACE');
% %     mask = dental_mask_polish(mask); % fill holes
% end
% rec.M = imresize3(mask,refNX,'nearest');
% rec.M = permute(rec.M,[3 1 2]); % ZN: back to original view

%%% CHECK GEOMTRIES
isAliased = isfield(rec.Enc.UnderSampling,'R') && any(rec.Enc.UnderSampling.R~=1);
if ~rec.Enc.UnderSampling.Flag && ~isAliased
%     xTest1 = mapVolume (recS.x, ones(NX(1:3)), MTS, MTy,[],[],'spline');
    xTest1 = mapVolume (recS.x, ones(NX(1:3)), MTS, MTy,[],[],'linear'); % debugging ZN
    %xTest2 = SWCC(multDimMea(rec.y,5:16),rec.S);
    xTest2 = RSOS(dynInd(resSub(rec.y,5:16),1,5),4);
    xTest2 = matchHist(xTest2,xTest1);
    labText = {sprintf('REF scan: %s',replace(recS.Names.Name,'_',' '));sprintf('ACQ scan (rsos): %s',replace(rec.Names.Name,'_',' '))};
    M1 = refineMaskExt(xTest1,0:.2:1,[],0,MT2MS(MTy) ); 
    %M2 = refineMaskExt(rec.M, 0:.2:1,[],10,MT2MS(MTy) ); 
    try
        plotND([], abs(cat(4,xTest1,xTest2)), defRange(xTest1),[],0,[],MTy,labText, M1,{2},1000,'Geometry check between REF and ACQ');
        saveFig(fullfile( rec.Names.pathOu ,'An-Ve_Sn','GeomTesting',rec.Names.Name));
    catch
        xTest1=[];xTest2=[];
    end
end
recS=[];
    
%%% ASSIGN B0 MAP
if ~isempty(fileb0)
    if existsFileVar(fileB0,{'recB0'})~=2; dat2B0(fileRef,fileB0);end
    recB0 = load(strcat(fileB0,'.mat'));recB0 =recB0.recB0; 
    MTB0 = recB0.Par.Mine.APhiRec;
    rec.B0 = mapVolume (dynInd(recB0.B0,1,4), ones(NX(1:3)) , MTB0, MTy);recB0=[]; %Only extract one B0 map   
end

%%% ASSIGN B1 MAP
if ~isempty(fileb1)
    if existsFileVar(fileB1,{'recB0'})~=2; dat2B1(fileRef,fileB1);end
    recB1 = load(strcat(fileB1,'.mat'));recB1 =recB1.recB1; 
    MTB1 = recB1.Par.Mine.APhiRec;
    rec.B = mapVolume (dynInd(recB1.B,1,4), ones(NX(1:3)) , MTB1, MTy);recB1=[];%Only extract one B1 map   
end

%%% RECONSTRUCT IMAGE AND WRITE TO INSPECT DEFAULT PARALLEL IMAGING RECON
if writeNIIFlag && ~rec.Enc.UnderSampling.Flag
    recW = rec;
    x = SWCC( multDimMea(recW.y,5), recW.S);%Sensitivity Weighted Coil Combination where repeats are averaged
    recW.Par.Mine.Modal=7;%Writing data to An-Ve = Volumetric encoding (3D)
    recW.Plan.Suff = '_SWCC'; recW.Plan.SuffOu='';recW.Names.pathOu = pathOu;
    recW.Dyn.Typ2Wri=zeros(1,50);
    recW.x=x;recW.Plan.Types{10}='x';recW.Plan.TypeNames{10}='Aq';recW.Dyn.Typ2Wri(10)=1;recW.Dyn.Typ2Rec=vertcat([],10);
    recW.Alg.OnlyJSON = 0;recW.Fail=0;
    recW.Par.Mine.Proce=1;recW.Par.Mine.Nat=1;recW.Alg.OverDec=ones(1,3);recW.Par.Mine.Signs=[];
    writeData(recW);
    fprintf('Sensitivity Weighted Coil Combination reconstruction saved to NIFTI file.\n')
end

