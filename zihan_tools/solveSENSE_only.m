function rec=solveSENSE_only(rec, Tenf, Tinit, cEnf, setOuterToZero)

if nargin<2 || isempty(Tenf);Tenf=[];end%If -1, just set zero motion
if nargin<3 || isempty(Tinit);Tinit=[];end
if nargin<4 || isempty(cEnf);cEnf=[];end
if nargin<5 || isempty(setOuterToZero);setOuterToZero=0;end

tsta=tic;
close all

%SET DEFAULT PARAMETERS
rec=disorderAlgorithm(rec);
rec=B0Algorithm(rec);
rec=pilotToneAlgorithm(rec);

%FILENAME
folderSnapshots=strcat(rec.Names.pathOu, filesep, 'An-Ve_Sn');%YB: Anatomical - Volumetric
[~,fileName]=fileparts(generateNIIFileName(rec));

%LOGGING INFO
if rec.Dyn.Log
    logDir = strcat(rec.Names.pathOu,filesep,'An-Ve_Log', filesep);
    logName = strcat(logDir, fileName,rec.Plan.Suff,rec.Plan.SuffOu, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~exist(logDir,'dir'); mkdir(logDir);end
    diary(logName) %eval( sprintf('diary %s ', logName))
end

if rec.Dyn.Debug>=2
    fprintf('\n=================================================================================\n');
    fprintf('Study name: %s \n', rec.Names.pathOu);
    fprintf('Running file: %s \n', fileName);
    fprintf('Running suff: %s \n', replace(rec.Plan.Suff,'_',''));
    fprintf('Running suffOu: %s \n', replace(rec.Plan.SuffOu,'_',' '));
%     fprintf('Date of acquisition: %s\n', rec.Names.Date);
    c=clock; fprintf('Date of reconstruction: %.2d/%.2d/%.4d\n', c(3:-1:1))
    fprintf('=================================================================================\n');
end

%%% REPORT PRE-PROCESSING
if rec.Dyn.Debug>=2
    fprintf('\nPRE-PROCESSED STEPS:\n');
    fprintf('Sensitivity file: %s\n', rec.Par.preProcessing.fileref);
    if isfield(rec.Par.preProcessing,'compressCoils');fprintf('Coil compression: %d --> %d\n', rec.Par.preProcessing.compressCoils.NChaOrig, size(rec.y,4));end
    if isfield(rec.Par.preProcessing,'resRec');fprintf('Down-sampling during DAT-conversion: %.2f/%.2f/%.2f --> %.2f/%.2f/%.2f\n', rec.Par.preProcessing.resRec.resOrig, rec.Par.preProcessing.resRec.resNew);end
    if isfield(rec.Par.preProcessing,'resRecRecon');fprintf('Down-sampling during REC-preparation: %.2f/%.2f/%.2f --> %.2f/%.2f/%.2f\n', rec.Par.preProcessing.resRecRecon.resOrig, rec.Par.preProcessing.resRecRecon.resNew);end
    if isfield(rec.Par.preProcessing,'supportReadout');fprintf('Readout extraction during DAT-conversion: %d --> %d (suppFOV = %.2f --> %.2f)\n', rec.Par.preProcessing.supportReadout.NOrig, rec.Par.preProcessing.supportReadout.NNew,rec.Par.preProcessing.supportReadout.suppFOV);end
    if isfield(rec.Par.preProcessing,'supportReadoutRecon');fprintf('Readout extraction during REC-preparation: %d --> %d (suppFOV = %.2f --> %.2f)\n', rec.Par.preProcessing.supportReadoutRecon.NOrig, rec.Par.preProcessing.supportReadoutRecon.NNew,rec.Par.preProcessing.supportReadoutRecon.suppFOV);end
end

%%% REPORT ENCODING
if rec.Dyn.Debug>=2
    fprintf('\nACQUISITION PARAMETERS:\n');
    fprintf('   AcqVoxelSize = %.2f/%.2f/%.2f mm\n',rec.Enc.AcqVoxelSize);
    fprintf('   FOVmm = %.2f/%.2f/%.2f mm\n',rec.Enc.FOVmm);
    fprintf('   AcqSize = %.0f/%.0f/%.0f\n',rec.Enc.AcqSize);
    fprintf('   FOVSize = %.0f/%.0f/%.0f\n',rec.Enc.FOVSize);
    fprintf('   Acquisition directions: %s\n',rec.Par.Scan.MPS);
    if rec.Enc.UnderSampling.Flag
        fprintf('   R = %.2f/%.2f\n',rec.Enc.UnderSampling.R);
    else
        fprintf('   R = %.2f/%.2f\n',[1 1]);
    end
    if rec.Enc.DISORDER.Flag; fprintf('   Sampling: DISORDER \n');else;fprintf('   Sampling: Linear \n');end
end
%%% FILTER IF SLAB DETECTION
[filterSize, isSlab, slabDim] = filterForSlab(rec.Enc.FOVmm);
isSlab = 0;
if isSlab 
    if rec.Enc.UnderSampling.Flag && isfield(rec.Enc.UnderSampling,'Acq')&& rec.Enc.UnderSampling.Acq.R(slabDim-1)~=1
        %Don't do anything
    else
        %Only filter if it's not an accelerated scan
        NY = size(rec.y);
        rec.y=bsxfun(@times,rec.y,ifftshift(buildFilter(NY(1:3),'tukey',filterSize,[],.4))); %Apodize y in image domain
        warning('Slab filtering activated. Make sure this is only done when estimating motion using this slab and not applying motion parameters.')
    end
end

%SHORTCUTS FOR COMMONLY USED STRUCTURE FIELDS AND INITIALIZERS
voxSiz=rec.Enc.AcqVoxelSize;%Acquired voxel size
parXT=rec.Alg.parXT;parXB=rec.Alg.parXB;
gpu=rec.Dyn.GPU;gpuIn=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);if gpuIn;gpuFIn=2;else gpuFIn=0;end 
gpu=single(gpuDeviceCount && ~rec.Dyn.BlockGPU);%YB: set gpu to single(gpuDeviceCount && ~rec.Dyn.BlockGPU)

if size(rec.y,5)>=parXT.maximumDynamics;fprintf('Problem too big to fit in memory (%d dynamics for a limit of %d)\n',size(rec.y,5),parXT.maximumDynamics);rec.Fail=1;return;end
on=cell(1,5);for n=1:5;on{n}=ones(1,n);end
typ2Rec=rec.Dyn.Typ2Rec;

%RECONSTRUCTION PLAN
[resPyr,L,estT,resIso,estB]=pyramidPlan(voxSiz,parXT.resolMax,parXT.NWend,parXT.accel);% YB: resPy is used for resAni, which in its turn is used in downSampingOperators.m
resPyrMax = [ .25 .5  1  ];
resPyr=rec.Alg.resPyr;%resPyr = resPyrMax(end+1-min(length(rec.Alg.resPyr),length(resPyr)):end) ;
resIso=sqrt(2)*((prod(voxSiz).^(1/3))./resPyr);%resIso = [resIso(end-(length(resPyr)-1):end)]
L = length(resIso);
estT = parXT.estT;%estT = ones(size(resIso));
if isfield(parXT,'estB');estB = parXT.estB;else;estB=zerosL(estT);end

rec = makePTParamsCompatible(rec, resPyr, resIso, estT, estB);
parXT=rec.Alg.parXT;parXB=rec.Alg.parXB;%Repeat because you recompute some fields in makePTParamsCompatible.m

%%% CHANGE LEVELS IF CALIBRATION MATRIX PROVIDED PTHandling
if (any(parXT.PT.usePT) && parXT.PT.Calibration.externalFit.Flag==1 && ~parXT.PT.Calibration.externalFit.forInit) || ~isempty(Tenf)
    L=2;
    estT=[0 0];
    estB=[1 1];
    resPyr=resPyr(end-1:end);
    resIso=resIso(end-1:end);
    rec.Alg.nExtern = [1 1];
    
    if any(parXT.PT.usePT)
        parXT.PT.subDiv=parXT.PT.subDiv(end-1:end);
        parXT.PT.usePT=[1 1];

        rec.Alg.parXT.PT.profileClustering.Flag = rec.Alg.parXT.PT.profileClustering.Flag(end-1:end);
            parXT.PT.profileClustering.Flag = parXT.PT.profileClustering.Flag(end-1:end);
        rec.Alg.parXT.PT.signalUsage.eigTh = rec.Alg.parXT.PT.signalUsage.eigTh(end-1:end);
            parXT.PT.signalUsage.eigTh = parXT.PT.signalUsage.eigTh(end-1:end);
        %rec.PT.NChaPTRec = rec.PT.NChaPTRec(end);%depends on eigTh so will autmoatically have right size
    else
        parXT.PT.usePT=[0 0];
        parXT.PT.subDiv=[0 0];
    end
    rec = makePTParamsCompatible(rec, resPyr, resIso, estT, estB);
end

if rec.Dyn.Debug>=2
    fprintf('\nRECONSTRUCTION PLAN:\n')
    fprintf('   Levels:             %s \n',sprintf(' %d ',1:L))
    fprintf('   Effective levels:   %s \n',sprintf(' %d ',(1:L) - sum(resPyr~=1)) )
    fprintf('   Resolution:         %s \n',sprintf(' %g ',resPyr))
    fprintf('   Motion estimation:  %s \n',sprintf(' %d ',estT))
    fprintf('   B0 estimation:      %s \n',sprintf(' %d ',estB))
    fprintf('   PT usage flag:      %s \n',sprintf(' %d ',parXT.PT.usePT))
    fprintf('   PT subdividing:     %s \n',sprintf(' %d ',parXT.PT.subDiv))
    fprintf('   PT clustering:      %s \n\n',sprintf(' %d ',parXT.PT.profileClustering.Flag))
end

%ROI COMPUTATION AND EXTRACTION IN THE READOUT DIRECTION
fprintf('warning:: Disable ROI extraction!!!!\n')
rec.Enc.ROI=computeROI(onesL(rec.M),[],[0 0 0],[1 0 0]); %YB: At this point non-logical values
if rec.Dyn.Debug>=2
    fprintf('ROI:\n%s',sprintf('   %d %d %d %d %d %d\n',rec.Enc.ROI'));
end
if ~isfield(rec.Par.Mine,'APhiRecOrig'); rec.Par.Mine.APhiRecOrig = rec.Par.Mine.APhiRec;end
if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = []; end
[~,MT] = mapNIIGeom([],rec.Par.Mine.APhiRec,'dynInd', {rec.Enc.ROI(1,1):rec.Enc.ROI(1,2),rec.Enc.ROI(2,1):rec.Enc.ROI(2,2),rec.Enc.ROI(3,1):rec.Enc.ROI(3,2)}, size(rec.M),[rec.Enc.ROI(1,4),rec.Enc.ROI(2,4),rec.Enc.ROI(3,4)]);

for n=typ2Rec';datTyp=rec.Plan.Types{n};
    if ~ismember(n,[5 12]);rec.(datTyp)=extractROI(rec.(datTyp),rec.Enc.ROI,1,1);end%YB: 1 at the end means 'forward'. At the end also inverse operation. Not extracting ROI for rec.N and rec.x
end

%PERMUTE PHASE ENCODES TO THE FIRST DIMENSIONS
perm=1:rec.Plan.NDims;perm(1:3)=[3 2 1]; %YB: At this point becomes:readout is in 3rd dimension; Becomes LR-AP-HF (see invert7T.m and main_recon.m)
for n=typ2Rec'        
    if ~ismember(n,[5 12]);datTyp=rec.Plan.Types{n};%YB: 5 and 12 are noise N and reconstruction x. rec.x does not exist yet, so don't permute.
        rec.(datTyp)=permute(rec.(datTyp),perm);
    end
end
if multDimSum(parXT.PT.usePT)>0 && isfield(rec,'PT') && ~isempty(rec.PT.pSliceImage)%PTHandling
    rec.PT.pSlice = permute(rec.PT.pSliceImage,perm);
    if gpuIn; rec.PT.pSlice = gpuArray(rec.PT.pSlice);end
end
voxSiz=voxSiz(perm(1:3));parXT.apod=parXT.apod(perm(1:3));
if ~isfield(rec.Par.Mine,'permuteHist'); rec.Par.Mine.permuteHist = [];end
rec.Par.Mine.permuteHist{end+1} = perm(1:4);

MS = voxSiz; 
[~,MT] = mapNIIGeom([],MT,'permute',rec.Par.Mine.permuteHist{end});%Don't change MS since already permuted in voxSiz
isValidGeom(MS,MT);

%COIL ARRAY COMPRESSION AND RECONSTRUCTED DATA GENERATION
NX=size(rec.M);
[S,y,eivaS]=compressCoils(rec.S,parXT.perc,rec.y);
NS=size(S);
if gpuIn;y=gpuArray(y);end
if rec.Dyn.Debug>=2 && ~isempty(parXT.perc)
    fprintf('\nCOIL COMPRESSION:\n')
    if parXT.perc(1)<1
        fprintf('   Number of compressed coil elements at%s%%: %d (%s )\n',sprintf(' %0.2f',parXT.perc*100),NS(4),sprintf(' %d',eivaS));
    else
        fprintf('   Number of compressed coil elements: %d (%s )\n',NS(4),sprintf(' %d',eivaS));
    end
    fprintf('   Number of coils used for intermediate image estimation: %d\n', eivaS(1));
    fprintf('   Number of coils used for motion estimation: %d\n', eivaS(2));
    fprintf('   Number of coils used for final image estimation: %d\n', eivaS(3));
end

%APODIZE + REARRANGE DATA + CORRECT FOR INTERLEAVED REPEATS %YB: apodisation means windowing
NY=size(y);NY(end+1:rec.Plan.NDims)=1;
y=bsxfun(@times,y,ifftshift(buildFilter(NY(1:3),'tukey',[10 10 1],gpuIn,parXT.apod))); %Apodize % YB: y in image domain
[y,NY]=resSub(y,5:rec.Plan.NDims);NY(end+1:5)=1;%Different averages along the 5th dimension
y=gather(y);

%ACQUISITION STRUCTURE
if rec.Dyn.Debug>=2;fprintf('\nSAMPLING TRAJECTORY:\n');end
NEchos=max(rec.Par.Labels.TFEfactor,rec.Par.Labels.ZReconLength);
NProfs=numel(rec.Assign.z{2});
if NProfs<=1 || NProfs~=numel(rec.Assign.z{3});fprintf('SEPARABLE Y-Z TRAJECTORIES. ALIGNED RECONSTRUCTION IS NOT PERFORMED\n');rec.Fail=1;return;end
if rec.Dyn.Debug>=2
    if sum(cdfFilt(abs(diff(rec.Assign.z{2}(:))),'med'))<sum(cdfFilt(abs(diff(rec.Assign.z{3}(:))),'med')) %YB: z{2}=PE1 and z{3}=PE2
        fprintf('   Slow direction is: %s\n',rec.Par.Scan.MPS(4:5));
    else; fprintf('   Slow direction is: %s\n',rec.Par.Scan.MPS(7:8));
    end
end

kTraj=zeros([NProfs 2],'single');
for n=1:2;kTraj(:,n)=rec.Assign.z{perm(n)}(:);end%YB: z{2}=PE1 and z{3}=PE2
if NEchos==1;NEchos=NProfs;end
NShots=NProfs/NEchos;
if mod(NShots,1)~=0;fprintf('NUMBER OF SHOTS: %.2f, PROBABLY AN INTERRUPTED SCAN\n',NShots);rec.Fail=1;return;end
NRepeats=NY(5);
if NShots<NRepeats
    NShots=NRepeats;
    NProfs=NShots*NEchos;
    kTraj=repmat(kTraj,[NRepeats 1]);
    NEchos=NProfs/NShots;
end
if rec.Dyn.Debug>=2
    fprintf('   Number of echoes: %d\n',NEchos);
    fprintf('   Number of shots: %d\n',NShots);
end
kTrajSS=reshape(kTraj,[NEchos NShots 2]);% YB: SS for snapshot I think
%PLOT TRAJECTORIES
% if rec.Alg.WriteSnapshots;visTrajectory(kTrajSS,2,strcat(folderSnapshots,filesep,'Trajectories'),strcat(fileName));end%YB:used to have ,rec.Plan.SuffOu
kRange=single(zeros(2,2));
for n=1:2;kRange(n,:)=rec.Enc.kRange{perm(n)};end
kShift=(floor((diff(kRange,1,2)+1)/2)+1)';
kSize=(diff(kRange,1,2)+1)';
kIndex=bsxfun(@plus,kTraj,kShift);
rec.Par.Mine.DISORDER.kTraj=gather(kTraj);
rec.Par.Mine.DISORDER.kRange=gather(kRange);
rec.Par.Mine.DISORDER.kShift=gather(kShift);
rec.Par.Mine.DISORDER.kSize=gather(kSize);
rec.Par.Mine.DISORDER.kIndex=gather(kIndex);

%STEADY-STATE CORRECTIONS
isSteadyState=0;
if mod(NShots,NRepeats)==0
    NShotsUnique=NShots/NRepeats;  
    isSteadyState=(NShotsUnique==1);%No shots have been identified
end
if rec.Dyn.Debug>=2;fprintf('   Steady-state: %d\n',isSteadyState);end
if isSteadyState
    kTrajs=abs(fct(kTraj)); %YB: By doing this frequency analysis, you don't actually need the sampling patter procided to the scanner!
    N=size(kTrajs);
    [~,iAs]=max(dynInd(kTrajs,1:floor(N(1)/2),1),[],1);    
    iM=min(iAs,[],2);
    NSamples=NProfs/(iM-1);
    if isfield(rec.Alg.parXT,'sampleToGroup') && ~isempty(rec.Alg.parXT.sampleToGroup);NSamples=rec.Alg.parXT.sampleToGroup;end%YB: added to force a certain number of samples to group
    if isfield(rec.Enc,'DISORDER') && ~rec.Enc.DISORDER.Flag;assert(~isempty(rec.Alg.parXT.sampleToGroup),'Shot definition should be defined by user in case of sequential sampling. Please provide value in rec.Alg.parXT.sampleToGroup.');end
    NSamplesPerAcquisitionSweep=NSamples; %YB added to comapre sample grouping to original acquisition
    NSweepsAcquisition=round(NProfs/(NSamplesPerAcquisitionSweep));
    if rec.Dyn.Debug>=2
        fprintf('   Number of acquisition sweeps: %d\n',NSweepsAcquisition);
        fprintf('   Temporal resolution per acquisition sweep: %.2f sec\n',NSamplesPerAcquisitionSweep*rec.Par.Labels.RepetitionTime(1)/1000);
    end
    NSweeps=round(NProfs/(NSamples*parXT.groupSweeps));% YB: Groupsweeps only for the first estimation, as afterwards the grouping will be based on the motoincompression
    NSamples=NProfs/NSweeps; %YB: Not necessary a round number as profiles per segment might slightly deviate    
    if rec.Dyn.Debug>=2
        fprintf('   Number of acquisition sweeps to group: %d\n',parXT.groupSweeps);
        fprintf('   Number of sweeps after grouping: %d\n',NSweeps);
        fprintf('   Number of samples per sweep: %.1f\n',NSamples);   
        fprintf('   Temporal resolution per sweep: %.2f sec\n',NSamples*rec.Par.Labels.RepetitionTime(1)/1000);  
    end
else
    NSamples=NEchos;
    NSweeps=NShots;
    NSamplesPerAcquisitionSweep=NSamples; %YB added to comapre sample grouping to original acquisition
    NSweepsAcquisition=round(NProfs/(NSamplesPerAcquisitionSweep));  
end
NStates=NSweeps;

%SWEEP SUBDIVISION AND TIME INDEX
sweepSample=ceil(NSweeps*(((1:NProfs)-0.5)/NProfs));% YB: I think the 0.5 is just a small number (compared to NProfs) so that ceil returns the right numbers
stateSample=sweepSample;
if isfield(rec.Enc, 'stateSample') && ~isempty(rec.Enc.stateSample)
    sweepSample = rec.Enc.stateSample;
    stateSample = sweepSample;
    NSweeps = max(stateSample);
    NStates=NSweeps;
    NSamples=NProfs/NSweeps; 
    NSamplesPerAcquisitionSweep=NSamples; %YB added to comapre sample grouping to original acquisition
    NSweepsAcquisition=round(NProfs/(NSamplesPerAcquisitionSweep)); 
    fprintf('                      !!!!!!!!!!!  WARNING:: Changed statesample. !!!!!!!!!!!!!!!\n')
    if rec.Dyn.Debug>=2
        %fprintf('                      Number of acquisition sweeps to group: %d\n',parXT.groupSweeps);
        fprintf('                      Number of sweeps after grouping: %d\n',NSweeps);
        fprintf('                      Number of samples per sweep: %.1f\n',NSamples);   
        fprintf('                      Temporal resolution per sweep: %.2f sec\n',NSamples*rec.Par.Labels.RepetitionTime(1)/1000);  
     	fprintf('                      !!!!!!!!!!!  WARNING:: Changed statesample. !!!!!!!!!!!!!!!\n')
    end
end

timeIndex=zeros([NY(1:2) NY(5)]);
if gpu;timeIndex=gpuArray(timeIndex);end
for n=1:size(kIndex,1)
    for s=1:NY(5) % YB: over repeats, so if filled for one, break as they all contain same timing
        if timeIndex(kIndex(n,1),kIndex(n,2),s)==0
            timeIndex(kIndex(n,1),kIndex(n,2),s)=n;
            break;
        end
    end
end
for m=1:2;timeIndex=ifftshift(timeIndex,m);end%YB: This goes to downSamplingOperators.m which expects a sampling operator A in the non-shifted fft domain
if isSteadyState
    timeSample=(0:NProfs-1)*rec.Par.Labels.RepetitionTime(1)/1000;
else
    TPerShot=rec.Par.Labels.ScanDuration/NShots;    
    timeSample=(sweepSample-1)*TPerShot;
    if strcmp(rec.Par.Scan.Technique,'TSE') || strcmp(rec.Par.Scan.Technique,'TIR')
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*(rec.Par.Labels.TE(1)/(1000*NEchos/2));
    else
        timeSample=timeSample+repmat(0:NEchos-1,[1 NShots])*rec.Par.Labels.RepetitionTime(1)/1000;
    end
end
sweepToSample=cell(1,NSweeps);
sampleToSweep=1:NSweeps;
for n=1:NSweeps;sweepToSample{n}=n;end%YB: 

%CREATE RECONSTRUCTION DATA ARRAY
NX=size(rec.M);NX(end+1:3)=1;
rec.x=zeros(NX,'single');%Uncorrected non-regularized
rec.d=zeros([NX L-sum(resPyr~=1)],'single');%Corrected non-regularized  %YB: For all the last (full resolution levels), the reconstructions are stored
if parXT.computeCSRecon;rec.r=zeros([NX L-sum(resPyr~=1)],'single');end %Corrected regularized
x=zeros(NX,'single');if gpuIn;x=gpuArray(x);end
typ2RecI=[12;16]; % 12-> Reconstruction (x) 16-> Volumetric alignment (d)
if isfield(rec.Alg,'bockWrAq') && rec.Alg.bockWrAq==1
    typ2RecI=[16]; 
    warning('Disabled saving Acq NIFTI files since it is the same for all the options for MoCo.');
end
if parXT.computeCSRecon;typ2RecI=[typ2RecI;18];end %18-> Filtered reconstruction (Re)
if parXB.useTaylor;typ2RecF=29;else; typ2RecF=[];end %29-> Linear coefficients of B0 fields with motion (Dl)
for n=typ2RecI'
    if ~any(rec.Dyn.Typ2Rec==n);rec.Dyn.Typ2Rec=vertcat(rec.Dyn.Typ2Rec,n);rec.Dyn.Typ2Wri(n)=1;end %YB: TODO should add typ2RecF here as well
end
typ2Rec=rec.Dyn.Typ2Rec;

%INITIALIZE TRANSFORM, NUMBER OF EXTERN ITERATIONS, TOLERANCES, IDENTIFY STEADY STATE SEQUENCES 
%Motion
rec.T=zeros([on{4} NSweeps 6],'single'); % defualt
if isfield(rec, 'init_T') && isfield(rec.init_T, 'useflag') % ZN: able to use other initial motion states
    if rec.init_T.useflag % ZN: for dental imaging, use the loaded fullFOV states as the initial motion states
        if isequal(size(rec.init_T.T), size(rec.T)) % ZN: double check the size of the motion states
            fprintf('Initialize the motion states with the loaded result! \n');
            rec.T = rec.init_T.T;
        end
    end
end

%B0 basis functions
if parXB.useSH>0
    SHNumCoef = size(SH_basis(2*on{3}, parXB.SHorder),2); 
    E.Db.c = zeros([on{4} NSweeps SHNumCoef],'like',real(x));E.Db.cr=E.Db.c;%E.NMs = NSweeps;
    E.Db.TE = rec.Par.Labels.TE * 1e-3;%Seconds
    if parXB.useSH==2; E.Dbs = E.Db; E = rmfield(E,'Db');end%Define fields in scanner frame
    E = updateBasis(E,parXB,NX,MT);
end

%B1 basis functions
if ~isfield(parXB,'useB1');parXB.useB1=0;end
if parXB.useB1
    E = updateBasisB1(E,parXB, NX,MT);
    E.B1m.c = zeros([on{4} NSweeps size(E.B1m.B,2) ],'like',real(x));E.B1m.cr=E.B1m.c;%E.NMs = NSweeps;
end

%Regularisation
if isfield(parXT,'corrFact'); E.corrFact = parXT.corrFact;end

%Optimisation
subT=ones(1,NStates);%To perform subdivisions of T
nExtern=99;
tolType={'Energy','RelativeResidual2Error'};%For CG--without vs with motion correction
tol=[parXT.tolerSolve 0];%For CG
nIt=[15 1];%For CG
nowithin=0;%To stop estimation of within motion

%ENERGY
En=cell(1,L);Re=cell(1,L);EnX=[];
mSt=cell(1,L);iSt=cell(1,L);
tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('\nTime arranging: %.3f s\n',tend);end

%SOLVE WITHOUT MOTION
tsta=tic;
resAni=single(on{3});
[~,indMinRes]=min(voxSiz);        
resAni(indMinRes)=voxSiz(indMinRes)/resPyr(end);
indNoMinRes=1:3;indNoMinRes(indMinRes)=[];
resAni(indNoMinRes)=max(resAni(indMinRes),voxSiz(indNoMinRes));   
resAni=voxSiz./resAni;
BlSz=max(sqrt(prod(voxSiz(1:2)./resAni(1:2))),1);%Block sizes baseline
if ~isfield(rec.Dyn,'GPUbS'); rec.Dyn.GPUbS = [2 4];end
E.bS=round([rec.Dyn.GPUbS].^log2(BlSz));E.oS=on{2};%Block sizes for our GPU
if rec.Dyn.Debug>=2;fprintf('Block sizes:%s\n',sprintf(' %d',E.bS));end
E.Je=1;%To use joint encoding/decoding
E.Uf=cell(1,3);for n=1:3;E.Uf{n}.NY=NY(n);E.Uf{n}.NX=NX(n);end%Folding
EH.Ub=cell(1,3);for n=1:3;EH.Ub{n}.NY=NY(n);EH.Ub{n}.NX=NX(n);end%Unfolding
[E.UAf,EH.UAb]=buildFoldM(NX(1:3),NY(1:3),gpuIn,1);%Folding/Unfolding matrix (empty is quicker!)
R=[];%Regularizer
ncx=1:eivaS(1);

%MASKING
%rec.M = onesL(rec.M);fprintf('warning:: Mask set to ones!!!!!!')
CX.Ma=rec.M;
if rec.Alg.UseSoftMasking==0 && ~all(ismember(CX.Ma(:),[0 1])) %Need to convert rec.M into discrete mask
   parS.Otsu=0:.2:1;%Binary vector of components for multilevel extraction (it picks those components with 1)
   parS.nDilate=10;%Dilation for masking (in mm)
   parS.conComp=2;
   CX.Ma=refineMask(rec.M,parS,voxSiz);
end
%Report
% if ~isequal(multDimSize(rec.y,1:3), multDimSize(rec.M,1:3))
%     plotND([], rec.M,defRange(rec.M),[],0,[],MT,{'Soft mask from reference scan'},CX.Ma==1,{2},2,'Mask used for reconstruction');
% else
%     plotND([], RSOS(multDimMea(rec.y,5:16)),defRange(RSOS(multDimMea(rec.y,5:16))),[],0,[],MT,{'Uncorrected image'},CX.Ma==1,{2},2,'Mask used for reconstruction');
% end
% if ~exist(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask'),'dir'); mkdir(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask'));end
% if rec.Alg.WriteSnapshots; saveFig(fullfile( rec.Names.pathOu ,'An-Ve_Sn','Mask',rec.Names.Name)) ;end
%Make mask consistent
rec.M = CX.Ma;

if gpuIn;CX.Ma=gpuArray(CX.Ma);end
if rec.Alg.UseSoftMasking==2;CX.Ma(:)=1;end %UseSoftMasking==2 means not masking
discMa=all(ismember(CX.Ma(:),[0 1])); %YB: discrete masking - if every element is rather 0 or 1
if ~discMa %Soft masking
    fprintf('Using soft masking\n');
    M=buildFilter(2*[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX],'tukeyIso',1,gpuIn,1,1);
    CX.Ma=abs(filtering(CX.Ma,M,1));M=[];
    if size(S,6)>1 || isfield(rec,'W')%We do not unfold in the PE direction
        NMa=[E.Uf{1}.NX E.Uf{2}.NX E.Uf{3}.NX];
        for m=1:3
            if ~isempty(E.Uf{m});CX.Ma=fold(CX.Ma,n,E.Uf{m}.NX,E.Uf{m}.NY,E.UAf{m});end
        end         
        CX.Ma=resampling(CX.Ma,NMa,2);
    end
    if isfield(rec,'W');CX.Ma=bsxfun(@times,CX.Ma,rec.W);end
    Ti.la=CX.Ma;
    R.Ti.la=1./(abs(Ti.la).^2+0.01);%YB: This is an example where the regularisation is precomputed as R'R and not specifically computed as R'R in regularize.m  
    CX=[];
else
    fprintf('Using discrete masking\n');
    if all(ismember(CX.Ma(:),[1]));fprintf('Using mask of all ones\n');end
end

%SENSITIVITIES
E.Sf=dynInd(S,ncx,4);E.dS=[1 ncx(end)];
if gpuIn;E.Sf=gpuArray(E.Sf);end

%%% B0 Shim
if ~isfield(parXB, 'modelShim');parXB.modelShim=0;end
if ~parXT.exploreMemory && isfield(rec.Par.Labels,'Shim') && ( isfield(parXB, 'modelShim') && parXB.modelShim==1 )% Also computed when est(B)==0
    E = updateShim(E, rec.Par.Labels.Shim, NX, MT);
    E.Shim.TE = rec.Par.Labels.TE * 1e-3;%in seconds
%     plotND({abs(x),1,0}, E.Shim.B0,[-200 200],[],1,{[],2}, MT, {'Shim B0 field'},[],[],9);
end

%PRECONDITION
if ~discMa;P.Se=(normm(E.Sf,[],4)+R.Ti.la).^(-1);%Precondition
else; P.Se=(normm(E.Sf,[],4)+1e-9).^(-1);
end

yX=mean(dynInd(y,ncx,4),5);
if gpuIn;yX=gpuArray(yX);end
gibbsRing=parXT.UseGiRingi*parXT.GibbsRingi;
if gibbsRing~=0 && gibbsRing<=1;yX=filtering(yX,buildFilter(NY(1:3),'tukeyIso',[],gpuIn,gibbsRing));end
nX=nIt(1);
if parXT.exploreMemory;nX=0;end%To test flow without running the main methods       

if nX~=0
    fprintf('\nSENSE RECONSTRUCTION:\n')
    [xRes,EnX]=CGsolver(yX,E,EH,P,CX,R,x,nX,tolType{1},tol(1));
    rec.x=gather(xRes);
    xUncorr = gather(xRes);
end
if ~isempty(EnX)%Print and store final energy
    EnX=EnX/numel(yX);
    if rec.Dyn.Debug>=2;fprintf('Energy evolution last step:%s\n',sprintf(' %0.6g',sum(EnX,1)));end  
    EnX=[];
end;yX=[];

if gpuIn;y=gpuArray(y);end
Residuals=encode(rec.x,E)-y;
for m=1:2;Residuals=fftGPU(Residuals,m)/sqrt(NY(m));end
Residuals=normm(Residuals,[],3:4);
if rec.Alg.WriteSnapshots
%     % below: mask the mouth region
%     if strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=1:floor((1-parXT.redFOV)*NX(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=floor(1+parXT.redFOV*NX(3)):NX(3);else; warning('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);end
    % below: mask the brain region (for dental MRI test)
    if parXT.redFOVflex
        if strcmp(rec.Par.Labels.FatShiftDir,'H') || strcmp(rec.Par.Labels.FatShiftDir,'F')
            E.nF=floor(1+parXT.redFOV_lower*NX(3)):floor((1-parXT.redFOV_upper)*NX(3));
        else
            warning('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);e
        end
    else
        % below: mask the brain region (for dental MRI test)
        if strcmp(rec.Par.Labels.FatShiftDir,'H');E.nF=1:floor((1-parXT.redFOV)*NX(3));elseif strcmp(rec.Par.Labels.FatShiftDir,'F');E.nF=floor(1+parXT.redFOV*NX(3)):NX(3);else; warning('Readout direction is not FH, performance probably suboptimal\n');E.nF=1:NX(3);end
    end
    
    %visSegment(flip(permute(dynInd(rec.x,E.nF,3),[2 1 3]),2),[],2,1,[],[],'Uncorrected',strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.Suff,'_Aq'));
%     plotND([], dynInd(abs(rec.x),E.nF,3), defRange(rec.x), [], [], [], MT, {'Uncorrected'}, [], [], [], [], [], [],  strcat(folderSnapshots,filesep,'Reconstructions'),strcat(fileName,rec.Plan.Suff,rec.Plan.SuffOu,'_Aq'));
%     visResiduals([],[],[],Residuals,kTraj,2,{strcat(folderSnapshots,filesep,'ResidualsStates'),strcat(folderSnapshots,filesep,'ResidualsSpectrum')},strcat(fileName,rec.Plan.Suff,rec.Plan.SuffOu,sprintf('_l=%d',0)));
end
E.Sf=[];

tend=toc(tsta);if rec.Dyn.Debug>=2;fprintf('Time computing non-corrected reconstruction: %.3f s\n',tend);end

x=resampling(xRes,NX); %YB: Start first level with the uncorrected guess

%%% PILOT TONE
if multDimSum(parXT.PT.usePT)>0 && isfield(rec,'PT') %PTHandling
    fprintf('\nPILOT TONE PRE-PROCESSING:\n')
    %Initialise variables
    rec.PT.ssPTConversion=[]; rec.PT.ssPTConversion.Geom.APhiRecOrig = MT;rec.PT.ssPTConversion.Geom.N=size(x);%For transformation conversion
    if parXT.PT.Calibration.useRASMotion==2;rec.PT.ssPTConversion.Geom.useLog=1;end
    %Shift to be consistent with sampling
    for n=1:2; rec.PT.pSlice=fftshiftOperator(rec.PT.pSlice,2,1,n);end
    %Pre-process: MB handling/Filtering/Normalisation/Phase-Whitening/TODO:phase-referencing
    rec.PT = preProcessPT(rec.PT,parXT.PT,NY,kIndex,[],rec.Par.preProcessing, NEchos);
    %Create calibration matrix
    [rec.PT.AfFull, rec.PT.AbFull] = createCalibrationMat(max(rec.PT.NChaPTRec), parXT.PT.Calibration.calibrationOffset, rec.PT.useRealImag);
    %Plot
    visPTSignal(rec.PT.pSlice,'pSlice', NY, kIndex,2,[],[],182,sprintf('Pre-processed PT signal'));pause(eps);
    if rec.Alg.WriteSnapshots;saveFig(strcat(folderSnapshots,filesep,'PTSignal',filesep,strcat(fileName,rec.Plan.Suff,rec.Plan.SuffOu,'_PreProcessed') ));end
end

if isfield(rec.Par.Labels ,'H0') && rec.Par.Labels.H0 ==3;MTT=[-1 0 0 0;0 -1 0 0;0 0 1 0;0 0 0 1]; else; MTT = eye(4);end
rec.Par.Mine.APhiRec = MTT*rec.Par.Mine.APhiRec;
MT = MTT * MT;

%%% SAVE NIFTI
fprintf('Writing NIFTI files.\n');
pathOuTemp = rec.Names.pathOu;
pathOuNII = strcat(pathOuTemp, '/An-Aq/'); if ~exist( pathOuNII,'dir');mkdir(pathOuNII);end
file = rec.Names.Name;
fileSave = strcat( pathOuNII,file);
xW=[];xW{1} = permute(rec.x,[3 2 1]); 
MSW=[];MSW{1} = rec.Enc.AcqVoxelSize; 
MTW =[]; MTW{1} = rec.Par.Mine.APhiRec;
writeNII(fileSave, {'Aq_womsk'},xW, MSW, MTW);

end