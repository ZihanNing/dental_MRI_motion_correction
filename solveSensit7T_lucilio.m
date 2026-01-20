
function rec=solveSensit7T_lucilio(rec, recB)

%SOLVESENSIT7T Solves for the sensitivities in 7T data by building a virtual body coil.
%   * REC is a reconstruction structure. At this stage it may contain the
%   naming information (rec.Names), the status of the reconstruction
%   (.Fail), the .lab information (rec.Par), the fixed plan information 
%   (rec.Plan), the dynamic plan information (rec.Dyn), the data 
%   information (rec.(rec.Plan.Types)), the information for correction
%   (rec.Corr.(rec.Plan.Types)), the informaton for sorting
%   (rec.Assign(rec.Plan.Types)) and the enoding information (rec.Enc)
%   ** REC is a reconstruction structure with estimated sensitivities rec.S
%   and eigenmaps rec.W
%

if nargin < 2; recB=[];end
gpu=gpuDeviceCount;
fprintf('Lucilio version ESPIRIT will be running.\n');

%RECURSIVE ESPIRIT (from Lucilio's code)
rec.Alg.ThreeDEspirit=1;%To use 3D ESPIRIT
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=2*ones(1,3);else rec.Alg.parE.NC=2*ones(1,2);end%Resolution (mm) of calibration area to compute compression
if rec.Alg.ThreeDEspirit;rec.Alg.parE.NC=1.5*ones(1,3);else rec.Alg.parE.NC=1.5*ones(1,2);end%Resolution (mm) of calibration area to compute compression
if rec.Alg.ThreeDEspirit;rec.Alg.parE.eigSc=[0.98 0.3];else rec.Alg.parE.eigSc=[0.9 0.3];end%Possible cut-off for the eigenmaps in soft SENSE, used the first for mask extraction by thresholding the eigenmaps%Default was [0.85 0.3]
if rec.Alg.ThreeDEspirit;rec.Alg.parE.Ksph=200;else rec.Alg.parE.Ksph=50;end%Number of points for spherical calibration area, unless 0, it overrides other K's%Default was 200
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.125;else rec.Alg.parE.factScreePoint=0.25;end%Factor over threshold computed with scree point
if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=0.5;else rec.Alg.parE.factScreePoint=1;end%Factor over threshold computed with scree point
%if rec.Alg.ThreeDEspirit;rec.Alg.parE.factScreePoint=2;else rec.Alg.parE.factScreePoint=4;end%Factor over threshold computed with scree point
rec.Alg.parE.mirr=[0 0 0];%Whether to mirror along a given dimension%Default was [8 8 8]
rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
%rec.Alg.parE.virCo=0;%Flag to use the virtual coil to normalize the maps%Default was 1, use 5 to get the Siemens contrast 
rec.Alg.parE.factorBody=1;%Factor for virtual body coil in ESPIRIT (it was 1)
%if strcmp(field,'refscan')
    %rph.Alg.parE.mirr=[0 0 2];
    %rph.Alg.parE.factorBody=10;
%end

%SOLVE ESPIRIT (from Lucilip's code)
% rph.x=B;rph.y=S;rph.Enc.AcqVoxelSize=rec.Enc.S.AcqDelta;rph.Alg.parE=rec.Alg.parE;
rec.x=blockCompressCoils(rec.y,1);
rec=solveESPIRIT_lucilio(rec);


%POST
rec = gatherStruct(rec);
%%% PLOT PARAMETERS FOR LOGGING
rec.Alg.parE

