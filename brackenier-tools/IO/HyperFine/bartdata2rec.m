
function rec = bartdata2rec(fileName, suppFOV, isPT, shimDir, removeOversampling, noiseFileName, fileName2, resRec, useGPU, removeZeroSamples, targetR, idTE, idINV)

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
if nargin<7 || isempty(fileName2);fileName2=[];end
if nargin<8 || isempty(resRec);resRec=[];end
if nargin<9 || isempty(useGPU);useGPU=1;end
if nargin<10 || isempty(removeZeroSamples);removeZeroSamples=1;end
if nargin<11 || isempty(targetR);targetR=[];end
if nargin<12 || isempty(idTE);idTE=[];end
if nargin<13 || isempty(idINV);idINV=[];end

%% SET EXPORTING INFORMATION
if ~isempty(fileName)
    exportFlag=1;
    [exportFolder,acqName]=fileparts(fileName);
    exportFolder = fullfile(exportFolder,'Parsing_Sn'); if ~exist(exportFolder,'dir');mkdir(exportFolder);end
else
    exportFlag=0;acqName='';
end

useGPU = 0;
Date = '00/00/0000';

gpu = useGPU && (gpuDeviceCount>0 && ~blockGPU);
if gpu; gpuDevice(gpuDeviceCount);end%Take last one

typ2Conv={'y'};

%% Load data
ss = loadHF_BartData(fileName);
yAll = ss.y;
kAll = ss.k;
flipOrder = [2 1 3];
kAll = dynInd(kAll,flipOrder,1);%Because of geom
ss.N = ss.N (flipOrder);
ss.info.fov = ss.info.fov (flipOrder);

ss.MT(3,4) = -195;
ss.MT(3,3) = -ss.MT(3,3);

rrOrig = BARTKSpaceRange(ss.N);
rrData = zerosL(rrOrig);

isVD = contains(fileName, 'vd') || contains(fileName, 'HF');

 %% Extract data
if ~isempty(idTE);yAll = dynInd(yAll,idTE,6);end

k = dynInd(kAll,1,6);%First anyways since duplicates
k = k(:,:,:,:);
for i=1:3
    [~, ~, rrData(:,i)] = getRange(k(i,:,:));
    rrData(:,i) = [-1 1]* max( abs(rrData(:,i)) );
end

%% Remove samples that fall off equidistant grid
if ~isVD
    [idx, grid] = findOnGrid(k,rrData);
    k = dynInd(k, idx==1,3);
    yAll = dynInd(yAll, idx==1,3);
    NCart = [length(grid{1}) length(grid{2}) length(grid{3}) ];
else
    NCart = ss.N;
end

rrCart = cartesianKSpaceRange(NCart);

%% Correct traj
correctTraj=1;
if correctTraj
    for i=1:3; k = dynInd(k, i,1, rescaleND(dynInd(k,i,1),rrOrig(:,i),rrData(:,i))); end
end

%% Remove double sampled points
removeDuplicates = 1;
if removeDuplicates
    [idx] = removeDuplicateSamp(squeeze(k(:,1,:)));

    k = dynInd(k,idx,3);
    yAll = dynInd(yAll, idx,3);
end

%% Map onto Cartesian grid
if ~isVD
    [yCart, kCart, z] = convertToCartesianGrid(yAll,k, NCart,rrOrig);
    Lin = z(:,1);
    Par = z(:,2);
    
    kCart = squeeze(kCart(2:3,1,:));
    kCart = kCart - centerIdx(multDimSize(yCart,2:3).');
    m = sqrt(normm(kCart,[],1));
    th = .95;
    idx = abs(m)>th*max(abs(m));
    noise = yAll(1,:,idx,:,1,1);
    noise = resPop(noise,2:3,[],1:2);
else
    fprintf('Gridding using BART\n');
    yCart = bart(sprintf('nufft -a -x %d:%d:%d',NCart),k,yAll);
    for dim=1:3; yCart = fftc(yCart,dim);end
    
    Lin = mod([1:prod(NCart(2:3))]-1,NCart(2))+1;
    Par = makeBins(prod(NCart(2:3)),NCart(3));
    
    m = sqrt( (Lin-centerIdx(NCart(2))).^2  + (Par-centerIdx(NCart(3))).^2 ); 
end

%% DEFINE ENCODING PARAMETERS FOR DATA ACQUISITION
%Initialise reconstruction structure
rec=[];
rec.Names.Date=Date;

%Parameters regarding TWIX object
ND=16;

%Number of receiver channels
NCha=size(yAll,4);
%NCha=min(NCha,3);warning('Number of channels fixed to lower count.');

%Oversampling in both PE dimensions
ROOverSamp = 1;%TW.hdr.Dicom.flReadoutOSFactor;
PEOverSamp = [0 0];
%if ~isempty(TW.hdr.Meas.flPhaseOS); PEOverSamp(1) = TW.hdr.Meas.flPhaseOS;end
%if ~isempty(TW.hdr.Meas.flSliceOS); PEOverSamp(2) = TW.hdr.Meas.flSliceOS;end

%FOV array size 
NImageCols = NCart(1);
% NImageLins = TW.hdr.Config.NImageLins;%Attention, this is "recon FOV", not what you acquired!
% NImagePars = TW.hdr.Config.NImagePar;%TW.hdr.Meas.lImagesPerSlab%Attention, this is "recon FOV", not what you acquired!

%Measured K-Space array size 
NKSpaceCols = NCart(1);
% NKSpaceLins = TW.hdr.Config.NLinMeas;%TW.hdr.Config.RawLin  
% NKSpacePars = TW.hdr.Config.NParMeas;%TW.hdr.Config.RawPar  
assert(ROOverSamp*NImageCols==NKSpaceCols,'TWIX2rec:: Parameters involving over-sampling dealing not consistent.')

%Grid K-Space array size 
NKSpaceColsGrid = NCart(1);
% NKSpaceLinsGrid = TW.hdr.Meas.iPEFTLength;
% NKSpaceParsGrid = TW.hdr.Meas.i3DFTLength;
relRes = [1 1 ];
NKSpaceLinsGrid = NCart(2)*(1+PEOverSamp(1))*relRes(1);
NKSpaceParsGrid = NCart(3)*(1+PEOverSamp(2))*relRes(2);

if mod(NKSpaceLinsGrid,1)~=0; warning('KSpace array size for PE_1 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceLinsGrid,round(NKSpaceLinsGrid)); NKSpaceLinsGrid = round(NKSpaceLinsGrid);end
if mod(NKSpaceParsGrid,1)~=0; warning('KSpace array size for PE_2 based on over-sampling and relative resolution not an integer. Rounded array size used ( %.3f --> %d ).',NKSpaceParsGrid,round(NKSpaceParsGrid)); NKSpaceParsGrid = round(NKSpaceParsGrid);end

NImageLins = ss.N(2);%NKSpaceLinsGrid;%We don't care about the recon FOV. We reconstruct at the resol/FOV we acquire.
NImagePars = ss.N(3);%NKSpaceParsGrid;

%Padding to be performed in k-space
kSpacePad = zeros([2 3]);
%Padding at origin of array
kSpacePad(1,1) = 0;%centerIdx(NKSpaceColsGrid) - TW.image.centerCol(1);
kSpacePad(1,2) = 0;%centerIdx(NKSpaceLinsGrid) - TW.image.centerLin(1);
kSpacePad(1,3) = 0;%centerIdx(NKSpaceParsGrid) - TW.image.centerPar(1);
%Padding after array
kSpacePad(2,1) = 0;%NKSpaceColsGrid - TW.image.sqzSize(1) - kSpacePad(1,1);
kSpacePad(2,2) = 0;%NKSpaceLinsGrid - TW.image.sqzSize(3) - kSpacePad(1,2);
kSpacePad(2,3) = 0;%NKSpaceParsGrid - TW.image.sqzSize(4) - kSpacePad(1,3);
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
rec.Enc.UnderSampling.Flag = 0%%%if TW.hdr.MeasYaps.sPat.ucPATMode>1; rec.Enc.UnderSampling.Flag = 1;else; rec.Enc.UnderSampling.Flag = 0;end %ucPATMode: No accel (1) GRAPPA (2) or mSENSE (3)
if (rec.Enc.UnderSampling.Flag || ~isempty(targetR)) && removeZeroSamples
    if rec.Enc.UnderSampling.Flag;fprintf('Accelerated scan detected.\n');end
%     %Method of acceleration
%     if TW.hdr.MeasYaps.sPat.ucPATMode==2; rec.Enc.UnderSampling.accelMethod = 'GRAPPA'; else ; rec.Enc.UnderSampling.accelMethod = 'mSENSE'; end
%     %Type of method
%     if isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==2
%         rec.Enc.UnderSampling.accelType = 'Integrated'; 
%     elseif isfield(TW.hdr.MeasYaps.sPat,'ucRefScanMode') && TW.hdr.MeasYaps.sPat.ucRefScanMode==4
%         rec.Enc.UnderSampling.accelType = 'Non-Integrated'; 
%     else
%         rec.Enc.UnderSampling.accelType = []; 
%     end

    %Undersampling factors and corresponding k-space indices to extract
    %rec.Enc.UnderSampling.R = cat(2, TW.hdr.MeasYaps.sPat.lAccelFactPE, TW.hdr.MeasYaps.sPat.lAccelFact3D );
    rec.Enc.UnderSampling.R = [1 1]
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
if rec.Enc.UnderSampling.Flag && strcmp(rec.Enc.UnderSampling.accelType,'Non-Integrated')
     %rec.ACS = dynInd(TW.refscan,':',ND);
     %typ2Conv{end+1}='ACS';
end
    
%% NOISE READOUT
if existsFileVar(noiseFileName) || existsFileVar(strcat( noiseFileName,'.h5')) %Take noise samples from other data
    %Get the data that should contain the noise files
    if existsFileVar( noiseFileName,'rec')==2 % Converted .mat file exists
        recN = load(noiseFileName,'rec'); 
        recN = recN.rec;
    elseif existsFileVar(strcat( noiseFileName,'.h5')) % Convert .dat file on the fly
        fprintf('Converting noise file on the fly: %s\n',noiseFileName)
        evalc('recN = rrdf2Rec(noiseFileName,[],[],1,[],[],[],[],[],[],[],0);')
    end 
    %Assign it to the current rec structure if it exists
    if isfield(recN,'N')
        fprintf('Using noise samples from separate acquisition: %s\n',noiseFileName)
        rec.N = recN.N;
        rec.Par.preProcessing.noiseFileName = noiseFileName;
    else
        error('TWIX2rec:: No noise samples found in %s.mat provided/converted. Noise de-correlation disabled!',noiseFileName);%Noise field does not exist
    end
else %Use noise for current file
    rec.N = noise;%rec.N=dynInd(TW.noise,':',ND); 
    typ2Conv{end+1}='N';%Still need to pre-process the noise signal
    
end

if isfield(rec,'N') && gpu;rec.N=gpuArray(rec.N);end

%% K-SPACE READOUT
perm=1:ND;
    
%%% Read data
fprintf('Reading out channels: %d/%d\n',NCha,NCha);
%rec.y=dynInd(TW.image,{1:NCha ':'},[2 ND]); %RO-Channel-PE1-PE2
rec.y = yCart;

%%% Pad
if multDimSum(kSpacePad (1,:))>0; rec.y = padArrayND(rec.y, [kSpacePad(1,1) 0 kSpacePad(1,2) kSpacePad(1,3)],[],0,'pre'); end
if multDimSum(kSpacePad (2,:))>0; rec.y = padArrayND(rec.y, [kSpacePad(2,1) 0 kSpacePad(2,2) kSpacePad(2,3)],[],0,'post'); end

%%% Remove undersampling
if removeZeroSamples && rec.Enc.UnderSampling.Flag
    %Remove Lin/Par
    id = ismember(Lin , idxAcqPE{1}) & ismember(Par ,  idxAcqPE{2});
    Lin=Lin(id);
    Par=Par(id);
    id=[];

    %Remove zero-samples from k-space
    rec.y = dynInd(rec.y, idxAcqPE, 3:4);
end

%%% Permute
for n=1:length(typ2Conv)
    t2c=typ2Conv{n}; 
    rec.(t2c)=permute(rec.(t2c),perm);
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
%%% Fourier transforms in PE dimensions - move all channels to image domain
fprintf('Fourier transform in the phase encode dimensions.\n');
NY=size(rec.y);NY(end+1:ND)=1;
for n=2:3
    rec.y = fftshiftOperator(rec.y, 1, 0, n);
    rec.y = fftGPU(rec.y,n);
    rec.y = fftshiftOperator(rec.y, 2, 0, n);
    
    if isPT
        rec.PT.pSliceImage = fftshiftOperator(rec.PT.pSliceImage, 1, 0, n);
        rec.PT.pSliceImage = fftGPU(rec.PT.pSliceImage,n);
        rec.PT.pSliceImage = fftshiftOperator(rec.PT.pSliceImage, 2, 0, n);
    end
    if isfield(rec,'ACS')
        rec.ACS = fftshiftOperator(rec.ACS, 1, 0, n);
        rec.ACS = fftGPU(rec.ACS,n);
        rec.ACS = fftshiftOperator(rec.ACS, 2, 0, n);
    end
end

%% DE-CORRELATE NOISE
if isfield(rec,'N')
    fprintf('De-correlating channels from noise samples.\n');
    rec.N = dynInd(rec.N,1:NCha,4);
    [rec.y, rec.Par.preProcessing.deCorrNoise.covMatrix,rec.Par.preProcessing.deCorrNoise.ccm] = standardizeCoils(rec.y,rec.N);
    if isfield(rec,'ACS');[rec.ACS] = standardizeCoils(rec.ACS,rec.N);end
    rec.preProcessing.deCorrNoise.Flag=1;
else
    warning('bartdata2rec:: No noise-decorrelation applied! Might result in sub-optimal performance.');
    rec.Par.preProcessing.deCorrNoise.Flag=0;
end
if exist('covMatrix','var')
    h=figure('color','w'); imshow(abs(rec.Par.preProcessing.covMatrix),[]);title('Noise covariance matrix.');set(h,'color','w','Position',get(0,'ScreenSize'));
    if exportFlag; saveFig(fullfile(exportFolder ,'Noise','Covariance',acqName)); end
end

%% CONVERT TO PRS (PE-RO-SL), CONSISTENT WITH SIEMENS ORIENTATION CONVENTION
%%% DIMENSIONS IN IMAGE DOMAIN: RO-PE-SL
NY=size(rec.y);NY(end+1:ND)=1;
rec.Enc.FOVmm = [ss.info.fov(1)*1000 , ss.info.fov(2)*1000 , ss.info.fov(3)*1000];
rec.Enc.FOVmm = rec.Enc.FOVmm .* (1 + [0 PEOverSamp]); %Recale the FOV that is returned (not including the over-sampling)
rec.Enc.FOVSize = [NY(1), NImageLins, NImagePars]; %Work on the FOV of the acquired (over-sampled) data

rec.Enc.AcqVoxelSize=rec.Enc.FOVmm ./ rec.Enc.FOVSize;
rec.Enc.AcqSize=[ 2^(ROOverSamp==2)*NY(1) NY(2:3)];

%% GEOMETRY COMPUTATION FOR RECONSTRUCTION ARRAY
fprintf('Computing geometry for reconstruction array.\n');

%SCALING
rec.Par.Mine.Asca=diag([rec.Enc.AcqVoxelSize 1]);

%ROTATION 
rec.Par.Mine.Arot = diag([-1 -1 1 1]);

%TRANSLATION
translationRaw = [0 0 0 ];%slicePos(1:3); %in mm for center FOV (I think)
rec.Par.Mine.tablePos = 0;%TW.hdr.Config.GlobalTablePosTra; %in mm for table position
if ~isempty(rec.Par.Mine.tablePos);translationRaw(3) = translationRaw(3) + rec.Par.Mine.tablePos;end

rec.Par.Mine.Atra=eye(4); rec.Par.Mine.Atra(1:3,4) = translationRaw';

%Account for fact that tranRaw is referred to the center of FOV, not the first element in the array
N = [ inf, NImageLins,NImagePars];%Set inf to make sure this element is replaced
N(1) = NImageColsArray;%N(2) will depend if you removed over-sampling
orig = centerIdx(N(1:3));
orig = (orig - [0 0 0])';%.5 from fact that centreFOV not in logical units, but in physical (so need to go back to centre of first voxel). Not sure why not for RO/PE --CHECK
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
rec.Par.Mine.patientPosition = 'HFS';%TW.hdr.Dicom.tPatientPosition;
%rec.Par.Mine.PCS2RAS = getPCS2RAS(rec.Par.Mine.patientPosition);
%rec.Par.Mine.APhiRec = rec.Par.Mine.PCS2RAS  * rec.Par.Mine.APhiRec;% From PE-RO-SL to RAS
%[~,rec.Par.Mine.APhiRec] = mapNIIGeom([], rec.Par.Mine.APhiRec,'translate',[-1 -1 -1]);%TO CHECK:AD HOC - REVERSE ENGINEERED
rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;%Backup as how was read in

%DEDUCE ACQUISITION ORDER
[~,MT_acqOrder] = mapNIIGeom([], rec.Par.Mine.APhiRec, 'permute', [2 1 3], NY);%acquisitionOrder takes MT as input that is defined in the PRS to RAS
[rec.Par.Scan.MPS, rec.Par.Scan.Mine.slicePlane ] = acquisitionOrder(MT_acqOrder);clear MT_acqOrder;
rec.Par.Scan.Mine.RO = rec.Par.Scan.MPS(4:5);
rec.Par.Scan.Mine.PE1 = rec.Par.Scan.MPS(1:2);
rec.Par.Scan.Mine.PE2 = rec.Par.Scan.MPS(7:8); 
fprintf('Slice orientation: %s\n', rec.Par.Scan.Mine.slicePlane);
fprintf('    Readout: %s\n', rec.Par.Scan.Mine.RO );
fprintf('    1st Phase encode direction: %s\n', rec.Par.Scan.Mine.PE1 );
fprintf('    2nd Phase encode (slice) direction: %s\n\n', rec.Par.Scan.Mine.PE2 );

rec.Par.Labels.FoldOverDir = rec.Par.Scan.MPS(1:2);%First PE direction
rec.Par.Labels.FatShiftDir = rec.Par.Scan.MPS(5);%Positive RO (I think because it is the direction where fat has a negative shift)

%% Modify MS,MT if other FOV hyperfine

NOld = rec.Enc.FOVSize;
NNew = NY(1:3);
rec.Enc.FOVSize=NNew;
rec.Enc.AcqSize(1) = 2*rec.Enc.AcqSize(1);

rec.Enc.FOVmm= rec.Enc.AcqVoxelSize .* rec.Enc.FOVSize;

[rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec] = mapNIIGeom(rec.Enc.AcqVoxelSize , rec.Par.Mine.APhiRec,'padArrayND',(NNew-NOld)/2,NOld(1:3),NNew(1:3));

%% MAKE RO-PE-SL FOR ALIGNED SENSE PIPELINE CONVENTION
%%% Store the permutations performed from the PRS space
rec.Par.Mine.permuteHist = [];%rec.Par.Mine.permuteHist{1} = perm(1:4);

%% SEQUENCE INFORMATION
fprintf('Reading sequence parameters.\n');

rec.Par.Labels.RepetitionTime=ss.info.TR*1000;%%In ms
rec.Par.Labels.TE=ss.info.TE*1000;%%In ms
rec.Par.Labels.FlipAngle=0;%cat(2,TW.hdr.MeasYaps.adFlipAngleDegree{1:size(rec.y,10)});%In degrees
rec.Par.Labels.RFSpoiling=NaN;

rec.Par.Labels.ScanDuration = 0;%TW.hdr.MeasYaps.lScanTimeSec;%In seconds
rec.Par.Labels.dwellTime = 0;%2 * TW.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e+9;%In seconds %*2 /1e+9 since alDwellTime{1} in nanoseconds and for oversampling
rec.Par.Labels.Bandwidth= 0;%1 / rec.Par.Labels.dwellTime;%In Hz/FOV
rec.Par.Labels.voxBandwidth = 0;%rec.Par.Labels.Bandwidth / NImageCols;%In Hz/voxel in the readout direction
rec.Par.Labels.mmBandwidth = 0;%rec.Par.Labels.voxBandwidth / (rec.Enc.AcqVoxelSize(1)); %In Hz/mm

rec.Par.Labels.fieldStrength = 0.064; %In T

rec.Par.Labels.sequenceType = 'SE';
isShotBased=1;
if isShotBased    
    %Change TFE factor
    rec.Par.Labels.TFEfactor = ss.info.TR/ss.info.esp;
   
    rec.Par.Labels.NShots = length(Par)/rec.Par.Labels.TFEfactor; %rec.Par.Labels.TFEfactor = TW.hdr.MeasYaps.sFastImaging.lTurboFactor;%Number of samples per shot / TW.hdr.Config.TurboFactor
    fprintf('Shot-based sequence with %d shots (%d samples per shot).\n',rec.Par.Labels.NShots, rec.Par.Labels.TFEfactor)
    %assert(mod(rec.Par.Labels.NShots,1)==0,'TWIX2Rec:: Number of samples per shot is not an integer. Probably a corrupted scan.');
    rec.Par.Labels.RepetitionTime_Long=ss.info.TR*1000;
    rec.Par.Labels.RepetitionTime=ss.info.esp*1000;%Old TR was the TR per shot --> assumed TR is twice the echo time
    rec.Par.Labels.inversionTime = ss.info.TI*1000;%In ms
    rec.Par.Labels.ScanDuration = rec.Par.Labels.RepetitionTime_Long*rec.Par.Labels.NShots/1000;%In sec
else
    rec.Par.Labels.TFEfactor=length(Par); 
end
rec.Par.Labels.ZReconLength=1;

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
disScan = contains(fileName,'dis');

if disScan%DISORDER
    rec.Enc.DISORDER.Flag=true;

    ff = 1;
    kTrajs = cat(2,LinesOrig,PartitionsOrig);
    kTrajs = kTrajs - mean(kTrajs,1);
    kTrajs = resampling(kTrajs, ff*size(kTrajs,1),2);
    kTrajs = abs(fct(kTrajs)); %YB: By doing this frequency analysis, you don't actually need the sampling patter procided to the scanner!
   
    Ntt=size(kTrajs);
    [~,iAs]=max(dynInd(kTrajs,1:floor(Ntt(1)/2),1),[],1);    
    iM=min(iAs,[],2);
    iM = iM /ff;
    NSamples=(2*length(Lin)/(iM-1));
    NShots = length(Lin)/NSamples;
    rec.Par.Labels.NShots = NShots;
    
    rec.Enc.DISORDER.tileSize=estimateTileSize(cat(2,LinesOrig(:),PartitionsOrig(:)));

    %if alFree{4}==2;rec.Enc.DISORDER.tileOrdering = 'Lines in Partitions';else; rec.Enc.DISORDER.tileOrdering = 'Partitions in Lines';end
    %if alFree{5}==2;rec.Enc.DISORDER.shotOrdering = 'Pseudo-random';else; rec.Enc.DISORDER.shotOrdering = 'Random';end
    %if length(alFree)>5;rec.Enc.DISORDER.seed = alFree{6};end
            
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
rec.Enc.ellipticalShutter = NaN;%~isempty(TW.hdr.Meas.ucEnableEllipticalScanning) && TW.hdr.Meas.ucEnableEllipticalScanning==1;

%%% Plot trajectory
permTraj=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
kTraj=zeros([length(rec.Assign.z{2}) 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{permTraj(n)}(:);end
kTrajSS=reshape(kTraj,[length(rec.Assign.z{2}) 1 2]);
visTrajectory(kTrajSS,0,[],[],206);%h=figure('color','w'); imshow(timeMat , []);title('Time indices of k-space sampling'); set(h,'color','w','Position',get(0,'ScreenSize'));
sgtitle(replace(sprintf('\\textbf{Trajectory:} %s\n',acqName),'_',' '),'Interpreter','latex','Color',[1 1 1]*0,'FontSize',20);
if exportFlag; saveFig(fullfile(exportFolder ,'Trajectory',acqName)); end

%% MODIFY UNDERSAMPLING
if removeZeroSamples && (rec.Enc.UnderSampling.Flag || ~isempty(targetR)) && ~isequal(rec.Enc.UnderSampling.RTarget,rec.Enc.UnderSampling.R) 
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
end

%% RESAMPLE STRUCTURE TO A SPECIFIED RESOLUTION
if ~isempty(resRec)
    if length(resRec)==1; resRec = resRec * ones([1 3]); else; resRec(end+1:3) = rec.Enc.AcqVoxelSize((length(resRec)+1):3);end
    fprintf('Resampling reconstruction to new resolutions %.2fmm (RO) / %.2fmm (PE1) %.2fmm (PE2).\n', resRec)
    rec.Par.preProcessing.resRec.resOrig = rec.Enc.AcqVoxelSize;
    rec = resampleRec(rec,resRec); %Resamples the arrays but also changes the sampling pattern accordingly
    rec.Par.preProcessing.resRec.resNew = rec.Enc.AcqVoxelSize;
end

%% CLOSE
rec=gatherStruct(rec);

