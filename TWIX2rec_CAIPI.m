
function rec = TWIX2rec_CAIPI(TW, suppFOV, isPT, shimDir, removeOversampling, noiseFileName, fileName, resRec, useGPU, removeZeroSamples, targetR, idTE, idINV)

%TWIX2REC Inverts SIEMENS TWIX data to a (Philips) reconstruction structure.
%   * TW is a TWIX object.
%   * {SUPPFOV} is the support in the readout direction (defaults to empty).
%   * {ISPT} is a flag to detect PT signal in the over-sampled k-space (defaults to 0).
%   * {SHIMDIR} is directory with additional shim calibration information of the scanner.
%   * {REMOVEOVERSAMPLING} is a flag whether to remove the default oversampling (defaults to 1).
%   * {NOISEFILENAME} is the name of the file where noise samples are stored.
%   * {FILENAME} is the name of the file so that we can use it for logging/saving snapshots during parsing.
%   * {RESREC} is the resolution at which to save the reconstruction structure (defaults to native resolution).
%   * {USEGPU} is a flag to enable gpu functionality (defaults to 1).
%   * {REMOVEZEROSAMPLES} is a flag to remove the un-acquired samples (defaults to 1 - only for regular undersampling).
%   * {TARGETR} is the desired undersampling ratio of the acquisition.
%   * {IDTE} is the echo to extract.
%   * {IDINV} is the inverision contrast to extract.
%   ** REC is the converted reconstruction structure.
%
%   Yannick Brackenier 2023-12-04

if nargin<2 || isempty(suppFOV);suppFOV=[];end
if nargin<3 || isempty(isPT);isPT=0;end
if nargin<4 || isempty(shimDir);shimDir='';end%'/home/ybr19/Projects/B0Shimming/'
if nargin<5 || isempty(removeOversampling);removeOversampling=1;end
if nargin<6 || isempty(noiseFileName);noiseFileName='';end
if nargin<7 || isempty(fileName);fileName=[];end
if nargin<8 || isempty(resRec);resRec=[];end
if nargin<9 || isempty(useGPU);useGPU=1;end
if nargin<10 || isempty(removeZeroSamples);removeZeroSamples=1;end
if nargin<11 || isempty(targetR);targetR=[];end
if nargin<12 || isempty(idTE);idTE=[];end
if nargin<13 || isempty(idINV);idINV=[];end

assert(~isempty(TW),'TWIX2rec:: TWIX object cannot be empty.');
removeZeroSamples = 0

% >>>> ZN: recon to recon_resol (increase resolution) or recon to a larger
% FOVmm (more dense kspace sampling) >>>>
flag_reconresol = 1;  % 1-recon to the recon_resol (zero-padding); 0-recon to acq resol
if flag_reconresol == 1
    warning('Reconstruct the image to its recon_resol!');
else
    warning('Reconstruct the image to its acquired resolution!');
end
flag_reconflexFOVmm = 1;  % 1-recon to a larger FOVmm (fake PE oversampling, zero-padding); 0-recon to the acquired FOVmm
if flag_reconflexFOVmm
    flex_reconFOVmm = [250 250 208];
    warning('Reconstruct the image to a larger FOVmm [%s,%s,%s] as planning by adding zero-padded PE oversampling!',num2str(flex_reconFOVmm(1)),num2str(flex_reconFOVmm(2)),num2str(flex_reconFOVmm(3)));
else
    warning('Reconstruct the image to the FOVmm as it acquired!');
end
% <<<<< ZN: recon to recon resol or recon to a larger FOVmm <<<<<<<<<<<<<<

%% SET EXPORTING INFORMATION
if ~isempty(fileName)
    exportFlag=1;
    [exportFolder,acqName]=fileparts(fileName);
    exportFolder = fullfile(exportFolder,'Parsing_Sn'); if ~exist(exportFolder,'dir');mkdir(exportFolder);end
else
    exportFlag=0;acqName='';
end

%%% ======= TEMPORARY 
coilToPlot = [1:3];
useGPU = 0;
%%% ======= TEMPORARY 

TW.image.flagIgnoreSeg=1;%See mapVBVD tutorial

if any(TW.image.dataSize(5:end)>1); exportFlag=0;end%Plotting won't work yet for NY(5:end)>1
Date = TW.hdr.Phoenix.tReferenceImage0; Date = strsplit(Date,'.');Date = strcat( Date{end}(7:8),'/',Date{end}(5:6),'/',Date{end}(1:4));

gpu = useGPU && (gpuDeviceCount>0 && ~blockGPU);
if gpu; gpuDevice(gpuDeviceCount);end%Take last one

typ2Conv={'y'};

%% DEFINE ENCODING PARAMETERS FOR DATA ACQUISITION
%Initialise reconstruction structure
rec=[];
rec.Names.Date=Date;

%Parameters regarding TWIX object
ND=length(TW.image.dataSize); %Number of dimensions that can be called in the TWIX object (see public method subsref)

%Number of receiver channels
NCha=TW.image.NCha;
%NCha=min(NCha,3);warning('Number of channels fixed to lower count.');

%Samples in both phase encoding dimensions
Lin = TW.image.Lin;
Par = TW.image.Par;
if any(ismember(fieldnames(TW.image),'IsRawDataCorrect'));idDataCorr = TW.image.IsRawDataCorrect;else;idDataCorr=[];end
isShotBased = strcmp(TW.hdr.Meas.tScanningSequence,'GR\IR');% TW.hdr.MeasYaps.sFastImaging.lTurboFactor~=1;
isMP2RAGE = isShotBased && TW.image.dataSize(10)>1;
isMPRAGE = isShotBased && TW.image.dataSize(10)==1;
isMeGRE = ~isShotBased && TW.image.dataSize(8)>1;
isSPACE = contains(TW.hdr.Meas.tScanningSequence,'SE');
if isShotBased
    %Shot ordering information
    if TW.hdr.MeasYaps.sKSpace.unReordering==1; TFE = TW.image.dataSize(4); %Linear
    else; TFE = TW.image.dataSize(3);%Linear-Rotated
    end
    numTI=TW.image.dataSize(10);
    %Make sure full inversions are extracted
    idFullInversion = 1: (floor(length(Lin)/(numTI*TFE))*(numTI*TFE));
    Lin = dynInd(Lin,idFullInversion,2); Par = dynInd(Par,idFullInversion,2); 
    if ~isempty(idDataCorr);idDataCorr = dynInd(idDataCorr,idFullInversion,2);end
    idFullInversion=[];
    %Remove duplicate sampling information
    if numTI>1
        id = makeBins(Lin, length(Lin)/TFE );
        id = mod(id,numTI)==1;
        Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); 
        if ~isempty(idDataCorr);idDataCorr=dynInd(idDataCorr,id,2);end
        id=[];
    end
elseif isMeGRE
    %Remove duplicate sampling information
    numTE = TW.image.dataSize(8);
    id = 1:numTE:length(Lin);
    Lin = dynInd(Lin,id,2);Par = dynInd(Par,id,2); 
    if ~isempty(idDataCorr);idDataCorr=dynInd(idDataCorr,id,2);end
    id=[];
end

if TW.image.dataSize(6)>1
    Lin = Lin(TW.image.Ave==1);
    Par = Par(TW.image.Ave==1);
end
    
%Oversampling in both PE dimensions
ROOverSamp = TW.hdr.Dicom.flReadoutOSFactor;
PEOverSamp = [0 0];
if ~isempty(TW.hdr.Meas.flPhaseOS); PEOverSamp(1) = TW.hdr.Meas.flPhaseOS;end
if ~isempty(TW.hdr.Meas.flSliceOS); PEOverSamp(2) = TW.hdr.Meas.flSliceOS;end

%FOV array size 
NImageCols = TW.hdr.Config.NImageCols;%TW.hdr.Config.BaseResolution
% NImageLins = TW.hdr.Config.NImageLins;%Attention, this is "recon FOV", not what you acquired!
% NImagePars = TW.hdr.Config.NImagePar;%TW.hdr.Meas.lImagesPerSlab%Attention, this is "recon FOV", not what you acquired!

%Measured K-Space array size 
if isSPACE
    NKSpaceCols = TW.image.NCol;
else
    NKSpaceCols = TW.hdr.Config.NColMeas;%TW.hdr.Config.RawCol /     %NKSpaceCols=TW.image.NCol
end
% NKSpaceLins = TW.hdr.Config.NLinMeas;%TW.hdr.Config.RawLin  
% NKSpacePars = TW.hdr.Config.NParMeas;%TW.hdr.Config.RawPar  
assert(ROOverSamp*NImageCols==NKSpaceCols,'TWIX2rec:: Parameters involving over-sampling dealing not consistent.')

%Grid K-Space array size 
NKSpaceColsGrid = TW.hdr.Meas.iRoFTLength;
% NKSpaceLinsGrid = TW.hdr.Meas.iPEFTLength;
% NKSpaceParsGrid = TW.hdr.Meas.i3DFTLength;
if ~isSPACE
    if ~flag_reconresol
        NKSpaceLinsGrid = TW.hdr.Config.NImageLins*(1+PEOverSamp(1))*TW.hdr.MeasYaps.sKSpace.dPhaseResolution;
        NKSpaceParsGrid = TW.hdr.Config.NImagePar*(1+PEOverSamp(2))*TW.hdr.MeasYaps.sKSpace.dSliceResolution;
        %NKSpaceParsGrid=TW.hdr.Meas.lImagesPerSlab*(1+PEOverSamp(2))*TW.hdr.MeasYaps.sKSpace.dSliceResolution
    else
        % >>>> ZN: recon image to the recon_resol & larger FOVmm >>>>
        NKSpaceLinsGrid = TW.hdr.Config.NImageLins*(1+PEOverSamp_recon(1)); % do not consider dPhaseResolution or dSliceResolution; consider fake PE oversampling; zero-padding later
        NKSpaceParsGrid=TW.hdr.Meas.lImagesPerSlab*(1+PEOverSamp_recon(2)); %SPACE
        % <<<< ZN: recon image to the recon_resol & larger FOVmm <<<<
    end
else
    sliceInfo = TW.hdr.MeasYaps.sSliceArray.asSlice{1};
    rec.Enc.FOVmm = [sliceInfo.dReadoutFOV , sliceInfo.dPhaseFOV , sliceInfo.dThickness];
    if ~flag_reconresol % recon to acq_resol and acq_FOVmm
        NKSpaceLinsGrid = rec.Enc.FOVmm(2)/(rec.Enc.FOVmm(2)/NImageCols/TW.hdr.MeasYaps.sKSpace.dPhaseResolution)*(1+PEOverSamp(1));
        NKSpaceParsGrid = TW.hdr.Meas.lImagesPerSlab*(1+PEOverSamp(2))*TW.hdr.MeasYaps.sKSpace.dSliceResolution %SPACE        
    else
        % >>>> ZN: recon image to the recon_resol & larger FOVmm >>>>
        NKSpaceLinsGrid = rec.Enc.FOVmm(2)/(rec.Enc.FOVmm(2)/NImageCols)*(1+PEOverSamp(1));
        NKSpaceParsGrid = TW.hdr.Meas.lImagesPerSlab*(1+PEOverSamp(2)); %SPACE
        % <<<< ZN: recon image to the recon_resol & larger FOVmm <<<<
    end
end

if mod(NKSpaceLinsGrid,1)~=0; warning('KSpace array size for PE_1 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceLinsGrid,round(NKSpaceLinsGrid)); NKSpaceLinsGrid = round(NKSpaceLinsGrid);end
if mod(NKSpaceParsGrid,1)~=0; warning('KSpace array size for PE_2 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceParsGrid,round(NKSpaceParsGrid)); NKSpaceParsGrid = round(NKSpaceParsGrid);end

% >>>> ZN: recon image to a larger FOVmm or recon_rescol >>>>
% the NKSpaceLinGrid here is already been modified according to flag_reconresol
NImageLins = NKSpaceLinsGrid;
NImagePars = NKSpaceParsGrid;
% >>>> ZN: recon image to a larger FOVmm or recon_rescol >>>>

%Padding to be performed in k-space
kSpacePad = zeros([2 3]);
%Padding at origin of array
kSpacePad(1,1) = centerIdx(NKSpaceColsGrid) - TW.image.centerCol(1);
kSpacePad(1,2) = centerIdx(NKSpaceLinsGrid) - TW.image.centerLin(1);
kSpacePad(1,3) = centerIdx(NKSpaceParsGrid) - TW.image.centerPar(1);
%Padding after array
kSpacePad(2,1) = NKSpaceColsGrid - TW.image.sqzSize(1) - kSpacePad(1,1);
kSpacePad(2,2) = NKSpaceLinsGrid - TW.image.sqzSize(3) - kSpacePad(1,2);
kSpacePad(2,3) = NKSpaceParsGrid - TW.image.sqzSize(4) - kSpacePad(1,3);
rec.Par.preProcessing.kSpacePad = kSpacePad;
assert(all(kSpacePad(:)>=0),'TWIX2Rec:: Conversion pipeline can not yet deal with negative padding.');
applyPad = kSpacePad>0;
if any(applyPad(:,1));warning('TWIX2Rec:: Needed to pad k-space RO dimension with [ %d ] pre- and [ %d ] post-elements since size(TW.image) was not consistent with physically acquired Columns.',kSpacePad(1,1),kSpacePad(2,1));end
if multDimSum(applyPad(:,2:3))>0;warning('TWIX2Rec:: Needed to pad k-space PE dimensions with [ %.0f %.0f ] pre- and [ %.0f %.0f ] post-elements since size(TW.image) was not consistent with physically acquired Lines/Partitions.',kSpacePad(1,2),kSpacePad(1,3),kSpacePad(2,2),kSpacePad(2,3));end
assymEcho = multDimSum(kSpacePad(:,1))>0;if assymEcho;warning('TWIX2Rec:: Assymetric echo acquired.');end

%If padding in from of array, need to change the indices to extract
if applyPad(1,2); Lin = Lin+kSpacePad(1,2); end
if applyPad(1,3); Par = Par+kSpacePad(1,3); end

        
%Acceleration
if TW.hdr.MeasYaps.sPat.ucPATMode>1; rec.Enc.UnderSampling.Flag = 1;else; rec.Enc.UnderSampling.Flag = 0;end %ucPATMode: No accel (1) GRAPPA (2) or mSENSE (3)
if (rec.Enc.UnderSampling.Flag || ~isempty(targetR)) && removeZeroSamples
    if rec.Enc.UnderSampling.Flag;fprintf('Accelerated scan detected.\n');end
    %Method of acceleration
    if TW.hdr.MeasYaps.sPat.ucPATMode==2; rec.Enc.UnderSampling.accelMethod = 'GRAPPA'; else ; rec.Enc.UnderSampling.accelMethod = 'mSENSE'; end
    %Type of method
    if isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==2
        rec.Enc.UnderSampling.accelType = 'Integrated'; 
    elseif isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==4
        rec.Enc.UnderSampling.accelType = 'Non-Integrated'; 
    else
        rec.Enc.UnderSampling.accelType = []; 
    end

    %Undersampling factors and corresponding k-space indices to extract
    rec.Enc.UnderSampling.R = cat(2, TW.hdr.MeasYaps.sPat.lAccelFactPE, TW.hdr.MeasYaps.sPat.lAccelFact3D );
    if isempty(targetR)
        rec.Enc.UnderSampling.RTarget = rec.Enc.UnderSampling.R./(1+PEOverSamp);
    else
        rec.Enc.UnderSampling.RTarget = targetR;
    end
    
    idxAcqPE = cell(1,2);
    offsetAcq{1} = rec.Enc.UnderSampling.R(1)*floor((min(Lin)-1)/rec.Enc.UnderSampling.R(1));
    offsetAcq{2} = rec.Enc.UnderSampling.R(2)*floor((min(Par)-1)/rec.Enc.UnderSampling.R(2));
    idxAcqPE{1} = min(Lin)-offsetAcq{1}:rec.Enc.UnderSampling.R(1):NKSpaceLinsGrid;
    idxAcqPE{2} = min(Par)-offsetAcq{2}:rec.Enc.UnderSampling.R(2):NKSpaceParsGrid;
    
    %Save in rec to access later
    rec.Enc.UnderSampling.limPE1 = [idxAcqPE{1}(1) idxAcqPE{1}(end)];
    rec.Enc.UnderSampling.limPE2 = [idxAcqPE{2}(1) idxAcqPE{2}(end)];
    %tt = zeros([NKSpaceLinsGrid, NKSpaceParsGrid]); tt = dynInd(tt,idxAcqPE,1:2,1); figure;imshow(tt,[]);
end

%Assign and convert ACS data
% if rec.Enc.UnderSampling.Flag %&& strcmp(rec.Enc.UnderSampling.accelType,'Non-Integrated')
%      %rec.ACS = dynInd(TW.refscan,':',ND);
%      %typ2Conv{end+1}='ACS';
% end

% revised by ZN, to enable ACS line for coil sensitivity map estimation for
% CAIPI acceleration pattern
if rec.Enc.UnderSampling.Flag % ZN: to enable ACS
    rec.ACS = dynInd(TW.refscan,':',ND); % might just support for non-integrated case only - need to check (ZN)
    typ2Conv{end+1}='ACS';
%     if strcmp(rec.Enc.UnderSampling.accelType,'Non-Integrated') % ZN: for non-integrated case
%         rec.ACS = dynInd(TW.refscan,':',ND);
%         typ2Conv{end+1}='ACS';
%     elseif strcmp(rec.Enc.UnderSampling.accelType,'Integrated') % ZN: for integrated case
%         NACS = length(TW.refscan.dataSize);
%         rec.ACS = squeeze(dynInd(TW.refscan,':',NACS));
%         if length(size(rec.ACS)) > 4 % ZN: multi-echo or multi-inv
%             rec.ACS = rec.ACS(:,:,:,:,1); % ZN: use the first echo or inv
%         end
%         typ2Conv{end+1}='ACS';
%     end
end


%% GAIN FACTORS
fileNameDat = TW.image.filename; %.dat file must still co-exist with TW object to extract image data
[~,rec.Par.preProcessing.fftScale, rec.Par.preProcessing.rawCorrectionScale] = CoilScalingFactors(fileNameDat,TW);

%% NOISE READOUT

if 0%existsFileVar(noiseFileName) || existsFileVar(strcat( noiseFileName,'.dat')) %Take noise samples from other data
    %Get the data that should contain the noise files
    if existsFileVar( noiseFileName,'rec')==2 % Converted .mat file exists
        recN = load(noiseFileName,'rec'); 
        recN = recN.rec;
    elseif existsFileVar(strcat( noiseFileName,'.dat')) % Convert .dat file on the fly
        fprintf('Converting noise file on the fly: %s\n',noiseFileName)
        evalc('recN = dat2Rec(noiseFileName,[],[],1,[],[],[],[],[],[],[],0);')
    end 
    %Assign it to the current rec structure if it exists
    if isfield(recN,'N')
        fprintf('Using noise samples from separate acquisition: %s\n',noiseFileName)
        rec.N = recN.N;
        rec.Par.preProcessing.noiseFileName = noiseFileName;
    else
        warning('TWIX2rec:: No noise samples found in %s.mat provided/converted. Noise de-correlation disabled!',noiseFileName);%Noise field does not exist
    end
elseif isfield(TW,'noise') %Use noise for current file
    if ~isPT
        rec.N = TW.noise.unsorted;%rec.N=dynInd(TW.noise,':',ND); 
        typ2Conv{end+1}='N';%Still need to pre-process the noise signal
    else
        warning('TWIX2rec:: Noise samples ignored since PT signal is present. No external noise file was provided.');
    end  
end

if isfield(rec,'N') && gpu;rec.N=gpuArray(rec.N);end

%% K-SPACE READOUT
perm=1:ND;
perm(2:4)=[3:4 2]; %Re-arrange so that 1st = Readout (RO), 2nd = First PE direction (PE1 or Lines), 3rd = Second PE direction (PE2 or Partitions) and 4th = Channel
    
%%% Read data
fprintf('Reading out channels: %d/%d\n',NCha,TW.image.NCha);
rec.y=dynInd(TW.image,{1:NCha ':'},[2 ND]); %RO-Channel-PE1-PE2
if ~isempty(idTE);rec.y = dynInd(rec.y,idTE,8);end
if ~isempty(idINV);rec.y = dynInd(rec.y,idINV,10);end

%%% Pad
if multDimSum(kSpacePad (1,:))>0; rec.y = padArrayND(rec.y, [kSpacePad(1,1) 0 kSpacePad(1,2) kSpacePad(1,3)],[],0,'pre'); end
if multDimSum(kSpacePad (2,:))>0; rec.y = padArrayND(rec.y, [kSpacePad(2,1) 0 kSpacePad(2,2) kSpacePad(2,3)],[],0,'post'); end

%%% Apply raw data correction - TSE
if any(idDataCorr)
    fprintf('Correcting raw data for gain changes.\n');

    NPE = multDimSize(rec.y,3:4);
    idTemp = sub2ind(NPE, Lin, Par);
    rec.y = resSub(rec.y, 3:4);

    idDataCorrT =  resPop(idDataCorr(:),1,[],3) .* resPop(rec.Par.preProcessing.rawCorrectionScale(:),1,[],2); 
    idDataCorrT(idDataCorrT==0)=1;
    
    NCorr = multDimSize(dynInd(rec.y,1,1),1:4);
    corrFac = ones(NCorr,'like',rec.y);
    corrFac = dynInd(corrFac,idTemp,3, idDataCorrT);

    rec.y = rec.y .* corrFac;
    rec.y = resSub(rec.y,3,NPE);
    NPE = [];idTemp=[];corrFac=[];
end

%%% Remove undersampling
if removeZeroSamples && rec.Enc.UnderSampling.Flag; rec.y = dynInd(rec.y, idxAcqPE, 3:4);end

%%% Permute
for n=1:length(typ2Conv)
    t2c=typ2Conv{n}; 
    rec.(t2c)=permute(rec.(t2c),perm);
end

%%% Plot individual coil data in k-space
if exportFlag
    NY=size(rec.y);
    for s=1:NCha
        APhiRecPlot = getPCS2RAS(TW.hdr.Dicom.tPatientPosition) * dynInd( diag(ones([1 4])) , {1:3,1:3},1:2,convertNIIGeom( ones([1,3]), TW.image.slicePos(4:7,1)', 'qForm', 'sForm'));
        [~,APhiRecPlot] = mapNIIGeom([], APhiRecPlot,'permute', [2 1 3]);  
        if isPT && ismember(s,coilToPlot)
            %Plot full k-space
            plotND({angle(dynInd(rec.y,s,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .5*dynInd(abs(dynInd(rec.y,s,4)),centerIdx(NY(1:3)),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_3D_%d',acqName,s));
            %Plot central line of k-space
            line = dynInd(dynInd(rec.y,s,4),centerIdx(NY(2:3)),2:3).';
            visPTSignal (line, [], [],[], 1, [], [], 200, replace(sprintf('\\textbf{k-space:}\n %s - Coil %d',acqName,s),'_',' '), fullfile(exportFolder ,'CoilKSpace') , sprintf('%s_Line_%d',acqName,s),'k-space readout [\#]',{'y','a.u.'});
        end
    end
end
    
%%% To image domain in readout
fprintf('Fourier transform in the readout dimension.\n');
if gpu;rec.y=gpuArray(rec.y);end 
for n=1:length(typ2Conv)
    t2c=typ2Conv{n}; 
    NY=size(rec.(t2c));NY(end+1:ND)=1;  
    %Block size operation due to large array sizes in over-sampled domain
    blSz = NCha;
    for s = 1:blSz:NCha
        vS = s:min(s+blSz-1,NCha);
        for l=1
            rec.(t2c) = dynInd( rec.(t2c), vS, 4,  fftshiftOperator(dynInd(rec.(t2c),vS,4), 1, 0, l));
            rec.(t2c) = dynInd( rec.(t2c), vS, 4,  fftGPU(dynInd(rec.(t2c),vS,4),l));
            rec.(t2c) = dynInd( rec.(t2c), vS, 4,  fftshiftOperator(dynInd(rec.(t2c),vS,4), 2, 0, l));
        end  
    end
end

%%% Plot individual coil data for hybrid data
if exportFlag
    NY=size(rec.y);
    for s=1:NCha
        APhiRecPlot = getPCS2RAS(TW.hdr.Dicom.tPatientPosition) * dynInd( diag(ones([1 4])) , {1:3,1:3},1:2,convertNIIGeom( ones([1,3]), TW.image.slicePos(4:7,1)', 'qForm', 'sForm'));
        [~,APhiRecPlot] = mapNIIGeom([], APhiRecPlot,'permute', [2 1 3]);  
        if isPT && ismember(s,coilToPlot)
            %Plot full hybrid k-space
            %plotND({angle(dynInd(rec.y,s,4)),1,0},abs(dynInd(rec.y,1,4)),[0 .5*dynInd(abs(dynInd(rec.y,s,4)),centerIdx(NY(1:3)),1:3) -pi pi],[],0,{[],2},APhiRecPlot,{'Phase';'Magnitude'},[],[],100,replace(sprintf('\\textbf{Hybrid k-space:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilHybrid') , sprintf('%s_3D_%d',acqName,s));
            %Plot central line of hybrid k-space
            line = dynInd(dynInd(rec.y,s,4),centerIdx(NY(2:3)),2:3).';
            %visPTSignal (line, [], [],[], 1, [], [], 200, replace(sprintf('\\textbf{Hybrid k-space:}\n %s - Coil %d',acqName,s),'_',' '), fullfile(exportFolder ,'CoilHybrid') , sprintf('%s_Line_%d',acqName,s),'Image readout [\#]',{'y','a.u.'});
        end
    end
end

%% EXTRACT PART OF READOUT TO AVOID FFT IN BOTH PE DIRECTIONS

%%% Remove over-sampling
if removeOversampling
    fprintf('Removing oversampling.\n');
    rec.y=resampling(rec.y,NImageCols,2);%Remove over-encoding provided by Siemens (used e.g. in Pilot Tone) 
    if isfield(rec,'N');rec.N=resampling(rec.N,NImageCols,2);end
    if isfield(rec,'ACS');rec.ACS=resampling(rec.ACS,size(rec.ACS,1)/2,2);end
end
NImageColsArray=size(rec.y,1);%The FOV array kept in the readout direction

%%% Reduce FOV in readout dimension
if ~isempty(suppFOV)
    fprintf('Extracting part of FOV with range [%.2f --> %.2f].\n', suppFOV(1), suppFOV(2));
    NY = size(rec.y);
    vr=max(1,round(NY(1)*suppFOV(1))):min( NY(1),round(NY(1)*suppFOV(2))); %Only in the RO direction
    rec.y=dynInd(rec.y,vr,1); %Extract supported FOV
    if isfield(rec,'N');rec.N=dynInd(rec.N,vr,1);end %No need to remove oversampling for noise I think 
    %if isfield(rec,'ACS');rec.ACS=dynInd(rec.ACS,vr,1);end
    rec.Par.preProcessing.supportReadout.suppFOV = suppFOV;
    rec.Par.preProcessing.supportReadout.NOrig = NY(1);
    rec.Par.preProcessing.supportReadout.NNew = length(vr);
else
    vr=[];
end

if gpu;rec.y=gather(rec.y);end

%% TRANSFORM TO IMAGE DOMAIN IN PE DIMENSIONS
%%% Check if gpu can handle all channels
yMem = getSize(rec,1);%In bytes
if gpu
    ss=gpuDevice; 
    safetyFact=.85;%Safety margin for GPU memory
    if yMem > (safetyFact * ss.AvailableMemory); gpu=0; rec=gatherStruct(rec);warning('TWIX2rec:: GPU usage disabled since the required memory exceeds GPU memory.'); end%Disable gpu and take everything out of it
end
if gpu;rec.y=gpuArray(rec.y);end

%%% Fourier transforms in PE dimensions - move all channels to image domain
fprintf('Fourier transform in the phase encode dimensions.\n');
NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    rec.y = fftshiftOperator(rec.y, 1, 0, n);
    rec.y = fftGPU(rec.y,n);
    rec.y = fftshiftOperator(rec.y, 2, 0, n);
    
    if isfield(rec,'ACS') % ZN: to transfer parallel calibration (ACS) to image domain as well
        rec.ACS = fftshiftOperator(rec.ACS, 1, 0, n);
        rec.ACS = fftGPU(rec.ACS,n);
        rec.ACS = fftshiftOperator(rec.ACS, 2, 0, n);
    end
end

%%% Plot coil images to inspect later
if exportFlag
    if ~exist(fullfile( exportFolder ,'CoilImage'),'dir'); mkdir(fullfile(exportFolder ,'CoilImage'));end
    for s=1:coilToPlot
        xPlotRef = abs(dynInd(rec.y,{s,1},[4 8]));
        xPlot = angle(dynInd(rec.y,{s,1},[4 8]));
        plotND({xPlotRef,1,0},xPlot,[-pi pi defRange(xPlotRef) ],[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{Image:} %s - Coil %d',acqName,s),'_',' '),[],[],fullfile(exportFolder ,'CoilImage') , sprintf('3D_%s_%d',acqName,s));
        xPlotRef=[];
    end
    xRSOS = RSOS(multDimSum(rec.y,6:ND));
    plotND([],xRSOS,defRange(xRSOS),[],0,{[],2},APhiRecPlot,[],[],[],100,replace(sprintf('\\textbf{RSOS:} %s',acqName),'_',' '),[],[],fullfile(exportFolder ,'RSOS') , sprintf('%s',acqName));
end


%% DE-CORRELATE NOISE
fprintf('Noise de-correlating disabled.\n');
if 0%isfield(rec,'N')
    fprintf('De-correlating channels from noise samples.\n');
    rec.N = dynInd(rec.N,1:NCha,4);
    [rec.y, rec.Par.preProcessing.deCorrNoise.covMatrix,rec.Par.preProcessing.deCorrNoise.ccm] = standardizeCoils(rec.y,rec.N);
    if isfield(rec,'ACS');[rec.ACS] = standardizeCoils(rec.ACS,rec.N);end
    rec.preProcessing.deCorrNoise.Flag=1;
else
    warning('TWIX2rec:: No noise-decorrelation applied! Might result in sub-optimal performance.');
    rec.Par.preProcessing.deCorrNoise.Flag=0;
end
if exist('covMatrix','var')
    h=figure('color','w'); imshow(abs(rec.Par.preProcessing.covMatrix),[]);title('Noise covariance matrix.');set(h,'color','w','Position',get(0,'ScreenSize'));
    if exportFlag; saveFig(fullfile(exportFolder ,'Noise','Covariance',acqName)); end
end

%% CONVERT TO PRS (PE-RO-SL), CONSISTENT WITH SIEMENS ORIENTATION CONVENTION
sliceInfo = TW.hdr.MeasYaps.sSliceArray.asSlice{1};

%%% DIMENSIONS IN IMAGE DOMAIN: RO-PE-SL
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.FOVmm = [sliceInfo.dReadoutFOV , sliceInfo.dPhaseFOV , sliceInfo.dThickness];
rec.Enc.FOVmm = rec.Enc.FOVmm .* (1 + [0 PEOverSamp]); %Recale the FOV that is returned (not including the over-sampling)
rec.Enc.FOVSize = [NY(1), NImageLins, NImagePars]; %Work on the FOV of the acquired (over-sampled) data

rec.Enc.AcqVoxelSize=rec.Enc.FOVmm ./ rec.Enc.FOVSize;
rec.Enc.AcqSize=[ 2^(ROOverSamp==2)*NY(1) NY(2:3)];

%%% CHANGE TO PE-RO-SL                  
rec.y=gather(rec.y);if isfield(rec,'N');rec.N=gather(rec.N);end
perm=1:ND; perm(1:3)=[2 1 3];
rec.y = permute(rec.y,perm);
if isfield(rec,'ACS'); rec.ACS = permute(rec.ACS,perm); end
if isfield( rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
rec.Enc.AcqVoxelSize = rec.Enc.AcqVoxelSize(perm(1:3));%Need to change spacing as well

%% GEOMETRY COMPUTATION FOR RECONSTRUCTION ARRAY
fprintf('Computing geometry for reconstruction array.\n');
slicePos = TW.image.slicePos(:,1);

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

%ROTATION 
quaternionRaw = slicePos(4:7);
rec.Par.Mine.Arot = eye(4);
rec.Par.Mine.Arot(1:3,1:3) = convertNIIGeom( ones([1,3]), quaternionRaw', 'qForm', 'sForm');%PE-RO-SL to PCS

%TRANSLATION
translationRaw = slicePos(1:3); %in mm for center FOV (I think)
rec.Par.Mine.tablePos = TW.hdr.Config.GlobalTablePosTra; %in mm for table position
if ~isempty(rec.Par.Mine.tablePos);translationRaw(3) = translationRaw(3) + rec.Par.Mine.tablePos;end

rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

%Account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = [ NImageLins, inf ,NImagePars];%Set inf to make sure this element is replaced
N(2) = NImageColsArray;%N(2) will depend if you removed over-sampling
orig = centerIdx(N(1:3));
orig = (orig - [0 0 .5])';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
%The [0 0 .5] is also the difference you find between the gadgetron PACS exported data and Siemens PACS exported data
if ~isempty(vr); orig(2)=orig(2)-(vr(1)-1); end%YB:orig(2) since this is readout/only vr(1)-1 elements removed from array/minus sign since should be addition and below already negative sign

rec.Par.Mine.Atra(1:3,4)= rec.Par.Mine.Atra(1:3,4)...
                        - rec.Par.Mine.Arot(1:3,1:3)*rec.Par.Mine.Asca(1:3,1:3)*orig;
                    
%COMBINED MATRIX
rec.Par.Mine.MTT=eye(4);%YB: not sure what this does --> de-activated
rec.Par.Mine.APhiRec=rec.Par.Mine.MTT*rec.Par.Mine.Atra*rec.Par.Mine.Arot*rec.Par.Mine.Asca;

% % Alternative way (works): 
% [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'translate',orig);
% if ~isempty(vr);[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'dynInd',{':',vr,':'},NYOrig(1:3),NY(1:3));end

%MOVE TO RAS
rec.Par.Mine.patientPosition = TW.hdr.Dicom.tPatientPosition;
rec.Par.Mine.PCS2RAS = getPCS2RAS(rec.Par.Mine.patientPosition);
rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
[~,rec.Par.Mine.APhiRec] = mapNIIGeom([], rec.Par.Mine.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;%Backup as how was read in

%DEDUCE ACQUISITION ORDER
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane] = acquisitionOrder(rec.Par.Mine.APhiRec);
rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
fprintf('Slice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
fprintf('    Readout: %s\n', rec.Par.Scan.Mine.RO );
fprintf('    1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
fprintf('    2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

%% GEOMETRY COMPUTATION FOR AUTOCALIBRATION ARRAY
if isfield(rec,'ACS')
    fprintf('Computing geometry for autocalibration array.\n');
    YSize = [ NImageLins, NImageColsArray ,NImagePars];%Set inf to make sure this element is replaced
    ACSSize = multDimSize(rec.ACS,1:3);
    [rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS] = mapNIIGeom([], rec.Par.Mine.APhiRec,'resampling',[],YSize,ACSSize);%TO CHECK:AD HOC - REVERSE ENGINEERED
end

%% MAKE RO-PE-SL FOR ALIGNED SENSE PIPELINE CONVENTION
perm=1:ND; perm(1:3) = [2 1 3 ];
rec.y = permute(rec.y, perm);
if isfield(rec,'ACS');rec.ACS = permute(rec.ACS, perm);end
if isfield(rec,'PT'); rec.PT.pSliceImage = permute(rec.PT.pSliceImage, perm);end
[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'permute', perm);
if isfield(rec,'ACS');[rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS] = mapNIIGeom(rec.Enc.UnderSampling.ACSVoxelSize, rec.Par.Mine.APhiACS,'permute', perm);end

%%% Store the permutations performed from the PRS space
rec.Par.Mine.permuteHist = [];rec.Par.Mine.permuteHist{1} = perm(1:4);

%%% Make MPS consisitent with RO-PE-SL
MPStemp = rec.Par.Scan.MPS ;
rec.Par.Scan.MPS(1:2) = MPStemp(4:5);
rec.Par.Scan.MPS(4:5) = MPStemp(1:2);

%% SEQUENCE INFORMATION
fprintf('Reading sequence parameters.\n');

rec.Par.Labels.RepetitionTime=TW.hdr.MeasYaps.alTR{1}/1000;%In ms
numTE=size(rec.y,8);
rec.Par.Labels.TE=cat(2,TW.hdr.MeasYaps.alTE{1:numTE})/1000;%In ms
rec.Par.Labels.FlipAngle=cat(2,TW.hdr.MeasYaps.adFlipAngleDegree{1:size(rec.y,10)});%In degrees
if isfield(TW.hdr.MeasYaps.sFastImaging,'ulEnableRFSpoiling')
    if strcmp(TW.hdr.MeasYaps.sFastImaging.ulEnableRFSpoiling,'0x1');rec.Par.Labels.RFSpoiling=1;%Hexadecimal code
    else; rec.Par.Labels.RFSpoiling = TW.hdr.MeasYaps.sFastImaging.ulEnableRFSpoiling;
    end
else
    rec.Par.Labels.RFSpoiling=0;
end
rec.Par.Labels.ScanDuration = TW.hdr.MeasYaps.lScanTimeSec;%In seconds
rec.Par.Labels.dwellTime = 2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;%In seconds %*2 /1e+9 since alDwellTime{1} in nanoseconds and for oversampling
rec.Par.Labels.Bandwidth= 1 / rec.Par.Labels.dwellTime;%In Hz/FOV
rec.Par.Labels.voxBandwidth = rec.Par.Labels.Bandwidth / NImageCols;%In Hz/voxel in the readout direction
rec.Par.Labels.mmBandwidth = rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm

rec.Par.Labels.fieldStrength = TW.hdr.Dicom.flMagneticFieldStrength; %In T

rec.Par.Labels.sequenceType = TW.hdr.Meas.tScanningSequence;
if isShotBased
    %Change sequence name
    if isMPRAGE; rec.Par.Labels.sequenceType = strcat(rec.Par.Labels.sequenceType,'-MPRAGE');elseif isMP2RAGE;rec.Par.Labels.sequenceType = strcat(rec.Par.Labels.sequenceType,'-MP2RAGE');end
    
    %Change TFE factor
    if TW.hdr.MeasYaps.sKSpace.unReordering==1 %Linear --> Each shot covers all partition samples = conventional MP-RAGE
        rec.Par.Labels.shotLayout='Linear/shot-covers-PE2';
        rec.Par.Labels.TFEfactor = size(rec.y,3);
    elseif TW.hdr.MeasYaps.sKSpace.unReordering==256 %Linear rotated --> Each shot covers all line samples
        rec.Par.Labels.shotLayout='Linear-Rotated/shot-covers-PE1';
        rec.Par.Labels.TFEfactor = size(rec.y,2);%TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor  
    end
    rec.Par.Labels.NShots = length(Par)/rec.Par.Labels.TFEfactor; %rec.Par.Labels.TFEfactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor
    rec.Par.Labels.NShotsFS = prod(NY(2:3))/rec.Par.Labels.TFEfactor; %rec.Par.Labels.TFEfactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor
    fprintf('Shot-based sequence with %d shots (%d samples per shot).\n',rec.Par.Labels.NShots, rec.Par.Labels.TFEfactor)
    %assert(mod(rec.Par.Labels.NShots,1)==0,'TWIX2Rec:: Number of samples per shot is not an integer. Probably a corrupted scan.');
    rec.Par.Labels.RepetitionTime_Long=rec.Par.Labels.RepetitionTime;
    rec.Par.Labels.RepetitionTime=rec.Par.Labels.TE*2;%Old TR was the TR per shot --> assumed TR is twice the echo time
    rec.Par.Labels.inversionTime = cat(2,TW.hdr.MeasYaps.alTI{1:length(TW.hdr.MeasYaps.alTI)})/1000;%In ms
    rec.Par.Labels.ScanDuration = rec.Par.Labels.RepetitionTime_Long*rec.Par.Labels.NShots/1000;%In sec
else
    rec.Par.Labels.TFEfactor=length(Par); 
end
rec.Par.Labels.ZReconLength=1;

if isSPACE
    rec.Par.Labels.turboFactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;
    rec.Par.Labels.echoTrainDuration = TW.hdr.MeasYaps.sFastImaging.lEchoTrainDuration;
    rec.Par.Labels.echoSpacing = TW.hdr.MeasYaps.sFastImaging.lEchoTrainDuration/rec.Par.Labels.turboFactor;
    rec.Par.Labels.RepetitionTime_Long = TW.hdr.MeasYaps.alTR{1}/1000 ;%In ms;
    rec.Par.Labels.RepetitionTime = rec.Par.Labels.echoSpacing;%In ms;
    rec.Par.Labels.RFEchoTrainLength = TW.hdr.Meas.RFEchoTrainLength;
end

%% SAMPLING INFORMATION
fprintf('Reading sampling information.\n');
NY=size(rec.y);NY=NY(1:3);

%%% Take out under-sampling intervals
if rec.Enc.UnderSampling.Flag && removeZeroSamples 
    Lines = ceil(Lin/rec.Enc.UnderSampling.R(1));
    Partitions = ceil(Par/rec.Enc.UnderSampling.R(2));
else
    Lines = Lin;
    Partitions = Par;
end
if isShotBased; LinesOrig = Lines; PartitionsOrig=Partitions;end

%%% Flip sampling and offset (because Siemens fft vs. ifft convention)
if mod(NY(2),2)==0;offset =-1;else;offset=0;end
Lines = mod( Lines-1 + offset, NY(2) ) +1;
if mod(NY(3),2)==0;offset =-1;else;offset=0;end
Partitions = mod( Partitions -1 + offset, NY(3) ) +1;

%%% Shift to make compatible with DISORDER recon
kShift = floor((NY)/2) + 1; %Not NY+1 since used like this in solveXT (floor((NY)/2) + 1 == floor((diff(kRange,1,2)+1)/2)+1 ) 
rec.Assign.z{2}= (NY(2) - (Lines-1)) - kShift(2); % 2nd PE 
rec.Assign.z{3}= (NY(3) - (Partitions-1)) - kShift(3); % 3rd PE = slices
%ATTENTION: this mod() operator gives big spikes in the sampling. This is just a processing implication.

rec.Enc.kRange={[-NY(1) NY(1)-1],[-NY(2)/2 NY(2)/2-1],[-NY(3)/2 NY(3)/2-1]};

%%%DISORDER information
if isfield(TW.hdr.MeasYaps,'sWipMemBlock') && isfield(TW.hdr.MeasYaps.sWipMemBlock,'alFree');alFree = TW.hdr.MeasYaps.sWipMemBlock.alFree;;end
if exist('alFree','var') && alFree{1}==2%DISORDER
    alFree = TW.hdr.MeasYaps.sWipMemBlock.alFree;
    rec.Enc.DISORDER.Flag=true;
    if isShotBased %DISORDER parameters not reliable when having a shot based sequence (To be fixed)
        rec.Enc.DISORDER.tileSize=estimateTileSize(cat(2,LinesOrig(:),PartitionsOrig(:)),rec.Par.Labels.NShots);
        %LinesOrig=[];PartitionsOrig=[];       
    else
        rec.Enc.DISORDER.tileSize=cat(2,alFree{2},alFree{3});
    end
    if alFree{4}==2;rec.Enc.DISORDER.tileOrdering = 'Lines in Partitions';else; rec.Enc.DISORDER.tileOrdering = 'Partitions in Lines';end
    if alFree{5}==2;rec.Enc.DISORDER.shotOrdering = 'Pseudo-random';else; rec.Enc.DISORDER.shotOrdering = 'Random';end
    if length(alFree)>5;rec.Enc.DISORDER.seed = alFree{6};end
            
    kMax = pi./rec.Enc.AcqVoxelSize;%rad/mmm
    dK = 2*kMax./rec.Enc.AcqSize;%rad/mm
    rec.Enc.DISORDER.tileSize_radmm=rec.Enc.DISORDER.tileSize .* dK(2:3);%physical units
    rec.Enc.DISORDER.NShots=prod(rec.Enc.DISORDER.tileSize);
    rec.Enc.DISORDER.segmentDuration = size(rec.Assign.z{2},2)/rec.Enc.DISORDER.NShots*rec.Par.Labels.RepetitionTime/1000;%In s - Use rec.Assign.z since NY does not take into account shutter
    if isShotBased && rec.Par.Labels.NShots~=rec.Enc.DISORDER.NShots;warning('TWIX2rec:: DISORDER encoded shot-based sequence should have compatible shot/segment durations.');end
else
    rec.Enc.DISORDER.Flag=false;
end

%%% Shutter
rec.Enc.ellipticalShutter = ~isempty(TW.hdr.Meas.ucEnableEllipticalScanning) && TW.hdr.Meas.ucEnableEllipticalScanning==1;

%%% Plot trajectory
permTraj=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
kTraj=zeros([length(rec.Assign.z{2}) 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{permTraj(n)}(:);end
kTrajSS=reshape(kTraj,[length(rec.Assign.z{2}) 1 2]);
visTrajectory(kTrajSS,0,[],[],206);%h=figure('color','w'); imshow(timeMat , []);title('Time indices of k-space sampling'); set(h,'color','w','Position',get(0,'ScreenSize'));
sgtitle(replace(sprintf('\\textbf{Trajectory:} %s\n',acqName),'_',' '),'Interpreter','latex','Color',[1 1 1]*0,'FontSize',20);
if exportFlag; saveFig(fullfile(exportFolder ,'Trajectory',acqName)); end

%% ENCODING COMPONENTS INFO
rec.Enc.RFInfo=[]; 
if isfield(TW.hdr.Meas,'ucExcitMode');[rec.Enc.RFInfo.RF_EXC, rec.Enc.RFInfo.SS_EXC] = convertTWIX_ucExcitMode(TW.hdr.Meas.ucExcitMode);end
if isfield(TW.hdr.Meas,'ucInversion');[rec.Enc.RFInfo.RF_INV, rec.Enc.RFInfo.SS_INV] = convertTWIX_ucInversion(TW.hdr.Meas.ucInversion);end

%%% Add information in case provided
if exist('alFree','var')
    %%% Gradient info in case blocked
    rec.Enc.GradInfo=[]; 
    rec.Enc.GradInfo.RO = 1;
    rec.Enc.GradInfo.PE = 1;
    rec.Enc.GradInfo.SS = 1;
    if length(alFree)>9; rec.Enc.GradInfo.RO = alFree{10}==2;end
    if length(alFree)>10; rec.Enc.GradInfo.PE = alFree{11}==2;end  
    if length(alFree)>11; rec.Enc.GradInfo.SS = alFree{12}==2;end  

    %%% RF info in case blocked
    if length(alFree)>12 && alFree{13}~=2;rec.Enc.RFInfo.RF_EXC = 0;rec.Enc.RFInfo.SS_EXC='NA';end
end

%%% Reference voltage
TXSpec = TW.hdr.MeasYaps.sTXSPEC;
if isfield(rec.Enc.RFInfo,'RF_INV') && rec.Enc.RFInfo.RF_INV==1
    rec.Enc.RFInfo.refVolt_INV=0;
    for run=1:length(TXSpec.aRFPULSE)
        if isfield(TXSpec.aRFPULSE{run},'tName') && strcmp(TXSpec.aRFPULSE{run}.tName,'SLoopIRns4') && isfield(TXSpec.aRFPULSE{run},'flAmplitude')
            rec.Enc.RFInfo.refVolt_INV = TXSpec.aRFPULSE{run}.flAmplitude;
            break
        end
    end
end
if isfield(rec.Enc.RFInfo,'RF_EXC') && rec.Enc.RFInfo.RF_EXC==1
    rec.Enc.RFInfo.refVolt_EXC = 0;
    for run=1:length(TW.hdr.MeasYaps.sTXSPEC.aRFPULSE)
        if isfield(TXSpec.aRFPULSE{run},'tName') && strcmp(TXSpec.aRFPULSE{run}.tName,'SRFExcit') && isfield(TXSpec.aRFPULSE{run},'flAmplitude')
            rec.Enc.RFInfo.refVolt_EXC = TXSpec.aRFPULSE{run}.flAmplitude;
            break
        end
    end
end

%TW.hdr.Dicom.flTransRefAmpl

%% SHIM SETTINGS
fprintf('Reading shim information.\n');
%%% Extract shim currents
rec.Par.Labels.Shim.shimCurrents_au = zeros(13,1);
if isfield(TW.hdr.Dicom, 'lFrequency'); rec.Par.Labels.Shim.shimCurrents_au(1) = TW.hdr.Dicom.lFrequency;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetX'); rec.Par.Labels.Shim.shimCurrents_au(2) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetX;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetY'); rec.Par.Labels.Shim.shimCurrents_au(3) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetY;end
if isfield(TW.hdr.Phoenix.sGRADSPEC, 'asGPAData') && isfield(TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}, 'lOffsetZ'); rec.Par.Labels.Shim.shimCurrents_au(4) = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1}.lOffsetZ;end
if isfield(TW.hdr.Meas, 'alShimCurrent'); rec.Par.Labels.Shim.shimCurrents_au(5:13) = TW.hdr.Meas.alShimCurrent(1:9).';end
                                  
rec.Par.Labels.Shim.shimOrder = 3; %Given the 13 parameters in the 7T Terra console - for 3T this will be fine as the non-existent coils with have zero shim current                          

%%% Load correctionFactors and crossTerms
if ~isempty(shimDir)
    correctionFactorsFileName = strcat(shimDir,'Calibration', filesep, 'Results',filesep, 'correctionFactors.mat');
    if existsFileVar(correctionFactorsFileName)==1
        ss = load(correctionFactorsFileName); rec.Par.Labels.Shim.correctFact_mA2au = ss.correctFact_mA2au;%Multiplicative correction
        rec.Par.Labels.Shim.shimCurrents_mA = rec.Par.Labels.Shim.shimCurrents_au ./ rec.Par.Labels.Shim.correctFact_mA2au;
    else
       fprintf('No Calibration data (correctionFactrs.mat) found in the shim directory.\n')  
    end

    crossTermFileName = strcat(shimDir,'Calibration', filesep, 'Results',filesep, 'crossTerms.mat');
    if existsFileVar(crossTermFileName)==1 
        ss = load(crossTermFileName); rec.Par.Labels.Shim.crossTerms = ss.crossTerms;
    else
       fprintf('No Calibration data (crossTerms.mat) found in the shim directory.\n')  
    end
else
   fprintf('No shim directory provided and hence no Calibration data read out.\n') 
end

%% PREPARE TEMPORAL PT SIGNAL (OPTIONAL)
if isPT
    tt = rec.PT.pSliceImage;
        
    %[kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);
    %filterStruct =[];filterStruct.medFiltKernelWidthidx=400/ rec.Par.Labels.RepetitionTime; 
    %tt = filterPTSlice (tt, filterStruct, NY, kIndex);

    for n=2:3
        tt=fftshiftOperator(tt,2,1,n);
        tt=fftGPU(tt,n)/NY(n);
        tt=fftshiftGPU(tt,n);%fftshift since in solveXT there is an iffthift on timeIndex
    end

    [kIndex, idx, hitMat, timeMat] = PESamplesFromRecPTTest(rec);

    tt = resSub(tt,5:ND);
    tt = resSub(tt,2:3);
    rec.PT.pTimeTest=[];
    for i=1:size(tt,4)%Multiple repeats
       rec.PT.pTimeTest = cat(2, rec.PT.pTimeTest, permute(dynInd(tt, {idx,i}, [2 4]),[3 2 4 1]));
    end    
    
    coilPlot=1:min(5,size(rec.PT.pTimeTest,1));
    plotType = 1;% (1) magn/phase (2) real/image (3) all
    %Plot original PT signal
    visPTSignal(dynInd(rec.PT.pTimeTest,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 200, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',acqName)) ; end
    visPTSignal(normm(rec.PT.pTimeTest(:,:,rec.PT.idxMB),[],1), [], [],[], plotType, [], [], 201, replace(sprintf('\\textbf{Norm original PT signal:}\n %s',acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',strcat(acqName,'_Norm'))) ; end
    %visPTSignalImage (rec.PT.pTimeTest(:,1:100,rec.PT.idxMB), [], [],[], plotType, [], [], 202, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    %if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','Trace',strcat(acqName,'_Image'))) ; end

    %Pre-processed signal
    tt = rec.PT.pTimeTest;
    tt = bsxfun(@rdivide, tt , sqrt(normm(dynInd(rec.PT.pTimeTest,{':',rec.PT.idxMB},[1 3]),[],1)) );
    if rec.PT.preProcessing.isRelativePhase; tt = bsxfun(@minus, tt, multDimMea((tt),2) );end%This only works well if you have referenced the phase to a channel
    visPTSignal (dynInd(tt,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 203, replace(sprintf('\\textbf{Pre-processed PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','TracePreProcessed',acqName)) ; end
    %visPTSignalImage (tt(:,1:100,rec.PT.idxMB), [], [],[], plotType, [], [], 205, replace(sprintf('\\textbf{Pre-processed PT signal for first %d coils:}\n %s',length(coilPlot),acqName),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    
    %if exportFlag; saveFig(fullfile(exportFolder ,'PilotTone','TracePreProcessed',strcat(acqName,'_Image'))) ; end
end

%% ADD ADDITIONAL PT SIGNAL INFORMATION (OPTIONAL)
% if isPT
%     centreFOV = rec.Par.Mine.APhiRec * [ceil( (rec.Enc.FOVSize(1:3)-1 -.5 +1)/2) 1]';%-1 since NIFTI convention and not MATLAB 
%     centreFOV = rec.Par.Mine.APhiRec * [(rec.Enc.FOVSize(1:3)/2 -[ 0 0 .5]) 1]';%-1 since NIFTI convention and not MATLAB
%     RODir = rec.Par.Mine.APhiRec * [1;0;0;1] - rec.Par.Mine.APhiRec * [0;0;0;1];
% 
%     offsetFOV = multDimSum(RODir(1:3).*centreFOV(1:3))./sqrt(normm(RODir(1:3)));%Projection on RO direction
% 
%     rec.PT.transmitFreqHz = setPTFreq(...
%        rec.Par.Labels.Shim.shimCurrents_au(1) ,...
%        NImageCols,...
%        rec.Par.Labels.voxBandwidth ,...
%        rec.PT.factorFOV, ...
%        offsetFOV,...
%        rec.Enc.AcqVoxelSize(1))/1e6;
% end

%% MODIFY UNDERSAMPLING
if removeZeroSamples && (rec.Enc.UnderSampling.Flag || ~isempty(targetR)) && ~isequal(rec.Enc.UnderSampling.RTarget,rec.Enc.UnderSampling.R) % ZN: normally, for the CAIPI acceleration, we do not remove undersampling, so do not run this session
    %Store original parameters if ever needed
    rec.Enc.UnderSampling.Acq.FOVSize = rec.Enc.FOVSize;
    rec.Enc.UnderSampling.Acq.AcqSize = rec.Enc.AcqSize;
    rec.Enc.UnderSampling.Acq.R = rec.Enc.UnderSampling.R;
    rec.Enc.UnderSampling.Acq.FOVmm = rec.Enc.FOVmm;
    rec.Enc.UnderSampling.Acq.PEOverSamp = PEOverSamp;
    %Assign new values
    rec.Enc.FOVSize(2:3) = round(rec.Enc.UnderSampling.RTarget .* rec.Enc.AcqSize(2:3));
    rec.Enc.FOVmm = rec.Enc.FOVSize .* rec.Enc.AcqVoxelSize;
    rec.Enc.UnderSampling.R = rec.Enc.FOVSize(2:3)./rec.Enc.AcqSize(2:3);
    %Change orientation information accordingly
    for l=2:3
        [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'ipermute',rec.Par.Mine.permuteHist{1});
        if rec.Enc.UnderSampling.Acq.FOVSize(l)>rec.Enc.FOVSize(l)%Extract
            vr=(centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(l))  - centerIdx(rec.Enc.FOVSize(l)))+[1:rec.Enc.FOVSize(l)];
            dynIndParam = {':',':',':'};dynIndParam(l) = {vr};
            [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec,'dynInd',dynIndParam,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
            [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'dynInd',dynIndParam,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
        else%Pad
            padDim = [0 0 0];
            padDim(l) = (centerIdx(rec.Enc.FOVSize(l)) - centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(l))) ;%Only need padding before
            [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'padArrayND',padDim,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
            [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'padArrayND',padDim,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
        end
        [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'permute',rec.Par.Mine.permuteHist{1});
    end
    fprintf('Changing undersampling from R=%.1fx%.1f to R=%.2f-%.2f by changing the FOV size.\n', rec.Enc.UnderSampling.Acq.R(1),rec.Enc.UnderSampling.Acq.R(2),rec.Enc.UnderSampling.R(1),rec.Enc.UnderSampling.R(2));
else
    if flag_reconflexFOVmm % ZN: for the special case: we'd like to reconstruct the image to a larger FOVmm, then even with removeZeroSamples as 0, we still need to recalculate the FOVmm
        %Store original parameters if ever needed
        rec.Enc.UnderSampling.Acq.FOVSize = rec.Enc.FOVSize;
        rec.Enc.UnderSampling.Acq.AcqSize = rec.Enc.AcqSize;
        rec.Enc.UnderSampling.Acq.FOVmm = rec.Enc.FOVmm;
        rec.Enc.UnderSampling.Acq.PEOverSamp = PEOverSamp;
        %Assign new values
        rec.Enc.FOVmm = flex_reconFOVmm; % ZN: change to FOVmm to a larger FOVmm; while rec.Enc.AcqVoxelSize should remain unchanged
        rec.Enc.FOVSize = round(rec.Enc.FOVmm./rec.Enc.AcqVoxelSize); % the AcqVoxelSize is not changed
        %Change orientation information accordingly
        for l=2:3
        [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'ipermute',rec.Par.Mine.permuteHist{1});
        if rec.Enc.UnderSampling.Acq.FOVSize(l)>rec.Enc.FOVSize(l)%Extract
            vr=(centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(l))  - centerIdx(rec.Enc.FOVSize(l)))+[1:rec.Enc.FOVSize(l)];
            dynIndParam = {':',':',':'};dynIndParam(l) = {vr};
            [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize, rec.Par.Mine.APhiRec,'dynInd',dynIndParam,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
            [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'dynInd',dynIndParam,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
        else%Pad
            padDim = [0 0 0];
            padDim(l) = (centerIdx(rec.Enc.FOVSize(l)) - centerIdx(rec.Enc.UnderSampling.Acq.FOVSize(l))) ;%Only need padding before
            [rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'padArrayND',padDim,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
            [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'padArrayND',padDim,rec.Enc.UnderSampling.Acq.FOVSize,rec.Enc.FOVSize);
        end
        [~,rec.Par.Mine.APhiRecOrig]=mapNIIGeom([],rec.Par.Mine.APhiRecOrig,'permute',rec.Par.Mine.permuteHist{1});
    end
    fprintf('Changing reconstructed FOVmm from [%d,%d,%d] to [%d,%d,%d].\n', rec.Enc.UnderSampling.Acq.FOVmm(1),rec.Enc.UnderSampling.Acq.FOVmm(2),rec.Enc.UnderSampling.Acq.FOVmm(3),rec.Enc.FOVmm(1),rec.Enc.FOVmm(2),rec.Enc.FOVmm(3));
    end
end

%% RESAMPLE STRUCTURE TO A SPECIFIED RESOLUTION
if ~isempty(resRec)
    if length(resRec)==1; resRec = resRec * ones([1 3]); else; resRec(end+1:3) = rec.Enc.AcqVoxelSize((length(resRec)+1):3);end
    fprintf('Resampling reconstruction to new resolutions %.2fmm (RO) / %.2fmm (PE1) %.2fmm (PE2).\n', resRec)
    rec.Par.preProcessing.resRec.resOrig = rec.Enc.AcqVoxelSize;
    rec = resampleRec(rec,resRec); %Resamples the arrays but also changes the sampling pattern accordingly
    rec.Par.preProcessing.resRec.resNew = rec.Enc.AcqVoxelSize;
end

%% VARIA
rec.Varia.hdr = TW.hdr;
if isfield(rec,'noise')
    rec.Varia.noise.hdr = TW.noise;
    rec.Varia.noise.N = TW.noise.unsorted;
    [~, rec.Varia.noise.covMatrix,rec.Varia.noise.ccm] = standardizeCoils(ones([1 1 1 NCha],'single'), rec.Varia.noise.N);
end

%% CLOSE
rec=gatherStruct(rec);

%% REMOVE MULTIPLE CONTRASTS
removeContrasts = 0;
if removeContrasts && size(rec.y,8)>1
    warning('TWIX2rec:: Multiple echoes detected. Removed additional echoes to reduce memmory.');
    rec.y = dynInd(rec.y,1,8);
    rec.Par.Labels.TE = rec.Par.Labels.TE(1);
end
if removeContrasts && size(rec.y,10)>1
    warning('TWIX2rec:: Multiple inversion contrasts detected. Removed additional inversion contrasts to reduce memmory.');
    rec.y = dynInd(rec.y,1,10);
    rec.Par.Labels.inversionTime = rec.Par.Labels.inversionTime(1);
end

%% RESOURCES
%FFT array sizes
%https://www.magnetom.net/t/number-of-lines-partitions-vs-fftlength/3012
%https://www.magnetom.net/t/anisotropic-in-plane-resolution-in-ima-header/797

%Bandwidth
%https://www.magnetom.net/t/how-to-get-the-spectroscopy-bandwidth-value-from-meas-as/320/4
%https://www.magnetom.net/t/dwelltime-in-header-dat-and-rda-files-incorrect-for-spectroscopy/3342

%Geometry
%asSlice=TW.hdr.Phoenix.sSliceArray.asSlice{1}; %See invert7T_old to see how to use this
%sNormal=asSlice.sNormal; sPosition=asSlice.sPosition;

%Eddy currents
%rec.Varia.ECC = TW.hdr.Phoenix.sGRADSPEC.asGPAData{1,1};

