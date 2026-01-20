

function [rec] = solveSensit_ESPIRiT_BART(rec, recB, param, saveRaw, computeCoilCompres, removeB0, suff, sensPadding)

if nargin<2 || isempty(recB);recB=[];end
if nargin<3 || isempty(param);param=[];end
if nargin<4 || isempty(saveRaw);saveRaw=0;end
if nargin<5 || isempty(computeCoilCompres);computeCoilCompres=1;end
if nargin<6 || isempty(removeB0);removeB0=0;end
if nargin<7 || isempty(suff);suff='';end
if nargin<8 || isempty(sensPadding);sensPadding=[];end

%% ADD PATHS
addBARTPathKCL();

%% MAKE PARAMETERS COMPATIBLE
if ~isfield(param,'ACLines') || isempty(param.ACLines); param.ACLines=24;end
if ~isfield(param,'kernelSize') || isempty(param.kernelSize); param.kernelSize=5;end
if ~isfield(param,'nMaps') || isempty(param.nMaps); param.nMaps=1;end
if ~isfield(param,'deb') || isempty(param.deb); param.deb=0;end
if ~isfield(param,'cropValue') || isempty(param.cropValue); param.cropValue=0;end%Crop the sensitivities if the eigenvalue is smaller than {crop_value}.
if ~isfield(param,'NNew') || isempty(param.NNew); param.NNew=[];end

%% ESTIMATE SENSITIVITIES
%%% Prepare array
%NY = size(rec.y);
x = rec.y;
N = multDimSize(x,1:3);
if ~isempty(sensPadding);NEspirit = N+2*sensPadding;else;NEspirit=N;end
NCoils = size(x,4);

if size(x,5)>1 && removeB0
    fprintf('solveSensitGE_ESPIRiT_BART:: Extracting first TE and subtractin B0 induced phase.\n');
    %Estimate B0 from raw coil image
    for i=1:3;x = ifftc(x,i);end
    B0 = multiCoilB0mapping(dynInd(x,1:2,5), rec.Par.Labels.TEs(1:2)*1000,rec.MS);
    %Take first TE and remove B0 phase
    x = dynInd(x,1,5);
    x = x.* conj(exp(1i*2*pi*rec.Par.Labels.TEs(1)*B0));%Remove B0 phase from first echo
    for i=1:3;x = fftc(x,i);end
elseif size(x,5)>1
    fprintf('solveSensitGE_ESPIRiT_BART:: Extracting first TE.\n');
    %Take first TE
    x = dynInd(x,1,5);
end

scalingFac = eps;
blockVC = 0;
if ~isempty(recB) && ~blockVC  
    fprintf('solveSensitGE_ESPIRiT_BART:: Adding body coild with scaled magnitude.\n');
    xB = dynInd(recB.y,1,5);%Take first TE
    x = cat(4,scalingFac*xB,x);
elseif ~blockVC   
    fprintf('solveSensitGE_ESPIRiT_BART:: Adding virtual body coild with scaled magnitude.\n');
    xB = x;
    for i=1:3;xB = ifftc(xB,i);end
    xB = compressCoils(xB,1);
    xB = RSOS(xB,4).*sign(compressCoils(xB,1));
        %for i=1:3;x = ifftc(x,i);end
        %x = x ./ (xB +eps);
        %for i=1:3;x = fftc(x,i);end
        %blockVC = 1;
    for i=1:3;xB = fftc(xB,i);end
    x = cat(4,scalingFac*xB,x);
end

if ~isempty(sensPadding)
    for i=1:3;x = ifftc(x,i);end%To image domain
    x = resampling(x,NEspirit(1:3),2);
    for i=1:3;x = fftc(x,i);end%To k-space domain
end

%%% Force k-space samples outside ACS to be zero
if length(param.ACLines)~=3; param.ACLines = cat(2,param.ACLines,param.ACLines(end)*ones([1 3-length(param.ACLines)]));end
if length(param.kernelSize)~=3; param.kernelSize = cat(2,param.kernelSize,param.kernelSize(end)*ones([1 3-length(param.kernelSize)]));end

x = resampling(x,param.ACLines,2);%x already in shifted Fourier domain
x = resampling(x,NEspirit(1:3),2);

%%% Run BART
if ~isempty(param.cropValue)
    cmdBART=sprintf('ecalib -r %d,%d,%d -k %d,%d,%d -m %d -d %d -c %.1f ', param.ACLines, param.kernelSize,param.nMaps,param.deb,param.cropValue);
else
    cmdBART=sprintf('ecalib -r %d,%d,%d -k %d,%d,%d -m %d -d %d -a ', param.ACLines,param.kernelSize,param.nMaps,param.deb);
end

fprintf('Running ESPIRiT using BART command:\n%s\n',cmdBART)
rec.S = bart(cmdBART, x);    
rec.S = dynInd(rec.S,(1:NCoils)+single(~blockVC),4);
rec.S = resampling(rec.S,N(1:3),2);

%% REMOVE OTHER DATA
rec.y=[];
rec.Ref=[];
rec.Nav=[];

%% COIL COMPRESSION CALCULATIONS
if computeCoilCompres
    if ~blockVC;x = dynInd(x,(1:NCoils)+single(~blockVC),4);end
    cutoff = .45;
    passBand = .1;
    [~,dimZ] = max(abs(rec.MT(3,1:3)),[],2);
    if rec.MT(3,dimZ)<0;dimZ=-dimZ;end
    if dimZ>0
        ROIInd = round( (cutoff+passBand/2)*N(abs(dimZ))):N(abs(dimZ));
        ROSInd = 1: round( (cutoff-passBand/2)*N(abs(dimZ)));
    else
        ROIInd = 1:round( ((1-cutoff)-passBand/2)*N(abs(dimZ)));
        ROSInd = round( ((1-cutoff)+passBand/2)*N(abs(dimZ))):N(abs(dimZ));
    end
    ROI = zeros(N(1:3),'single');
    ROI = dynInd(ROI,ROIInd,abs(dimZ),1);
    ROS = zeros(N(1:3),'single');
    ROS = dynInd(ROS,ROSInd,abs(dimZ),1);
    perc = size(x,4);

    rec.CoilCompression = [];
    rec.CoilCompression.ROVir.ROIInd = ROIInd;
    rec.CoilCompression.ROVir.ROSInd = ROSInd;
    rec.CoilCompression.ROVir.criteria = 'ROS';
    for i=1:3;x = ifftc(x,i);end%To image domain
    x = resampling(x,N(1:3),2);
    [~,~,~,rec.CoilCompression.ccmROVir]=compressCoilsROVir(x,perc,[],ROI,ROS,'ROS');
    plotND([],RSOS(x),[],[],0,[],rec.MT,[],ROS+2*ROI,{4},100);
end

%% DOWNSAMPLE (to save memmory)
if ~isempty(param.NNew)
    NNew = round(param.NNew);
    rec.S = resampling(rec.S, NNew);
    [rec.MS, rec.MT] = mapNIIGeom(rec.MS,rec.MT, 'resampling',[],N,NNew);
    rec.N = NNew;
end

%%% Save
if ~isempty(saveRaw) && saveRaw==1
    save(strcat(rec.Names.pathOu,rec.Names.Name,suff,'_sens_BART.mat'),'rec','-v7.3');
end
%% also estimate compression and store the compression matrix
%% Documentions
% Usage: ecalib [-t f] [-c f] [-k ...] [-r ...] [-m d] [-S] [-W] [-I] [-1] [-P] [-v f] [-a] [-d d] <kspace> <sensitivites> [<ev-maps>]
% 
% Estimate coil sensitivities using ESPIRiT calibration.
% Optionally outputs the eigenvalue maps.
% 
% -t threshold      	This determined the size of the null-space.
% -c crop_value      	Crop the sensitivities if the eigenvalue is smaller than {crop_value}.
% -k ksize      	kernel size
% -r cal_size      	Limits the size of the calibration region.
% -m maps      	Number of maps to compute.
% -S		create maps with smooth transitions (Soft-SENSE).
% -W		soft-weighting of the singular vectors.
% -I		intensity correction
% -1		perform only first part of the calibration
% -P		Do not rotate the phase with respect to the first principal component
% -v variance      	Variance of noise in data.
% -a		Automatically pick thresholds.
% -d level      	Debug level
% -h		help

%Another example of how to use it in matlab: https://mrirecon.github.io/bart/examples.html


%Move to normalised Fourer space. (For ecalib you need to take iffshift in image domain first before dong fft and fftshift) THIS IS VERIFIED
% for i=1:3
%     x = ifftshiftGPU(x,i);
%     x = fftGPU(x,i)/sqrt(N(i));
%     x = fftshiftGPU(x,i);
% end
