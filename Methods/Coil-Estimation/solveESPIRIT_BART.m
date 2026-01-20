
function rec = solveESPIRIT_BART(rec)

%%% ADD BART PATH AND ENC
di=1;
rec.Alg.parE.BARTPath = addBARTPathKCL(di);

%%% COMPUTE PARAMETERS
%%% Define constansts and variables
parE = rec.Alg.parE;
N = size(rec.y);
ND = 3;%Number of dimensions for estimation

%%% Define calibration region
DeltaK=(1./(N(1:ND).*rec.Enc.AcqVoxelSize(1:ND)));%1/FOV - YB: in units of 1/mm
DeltaK=DeltaK(1:ND);

KMaxCalib = 1./(2*rec.Alg.parE.NC); %YB: in units of 1/mm
parE.NC=2*KMaxCalib./DeltaK;%Resolution (array size) of k-space points for calibration
parE.NC=min(parE.NC,N(1:ND)-1);%Never bigger than the image resolution %Now parE.NC is in units of voxels
parE.NC=parE.NC-mod(parE.NC,2)+1;%Nearest odd integer
%YB:Now NC is in units of voxels to include in AC region and not resolution anymore

%%% Define kernel size
parE.K=(1./parE.K)./DeltaK;%Resolution of target coil profiles 
parE.K=max(parE.K,parE.Kmin);
parE.K=min(parE.K,N(1:ND));%Never bigger than the image resolution

%%% CONCATENATE BODY COIL
scalingFac = eps('single');
x = cat(4,scalingFac*rec.x,rec.y);
X = x;
for i=1:3
    X = fftc(X,i);
end

%%% PREPARE BART COMMAND
ACLines = N(1:ND);%size of calibration k-space fopr ESPIRiT
kernelSize = parE.K;% kernel size
nMaps=1;
cropValue=0;
deb=2;

kernelSize = [4 4 4]
fprintf('Using BART ACS area of %d-%d-%d.\n',ACLines);
fprintf('Using BART kernel size of %d-%d-%d.\n',kernelSize);

%%% RUN BART
cmdBART=sprintf('ecalib -r %d,%d,%d -k %d,%d,%d -m %d -d %d -c %.1f', ACLines, kernelSize, nMaps, deb, cropValue);
fprintf('Running BART command: %s\n',cmdBART);
[S, W] = bart(cmdBART, gather(X));

%%% ASSIGN
rec.SBart = dynInd(S, 1, 5);%First eigenmap
rec.SBart = dynInd(rec.SBart, 2:size(rec.SBart,4), 4);%Exclude body coil
rec.WBart = W;

size(rec.SBart)

%%% REMOVE BART FROM ENV
di=0;
addBARTPathKCL(di);

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
