
function [xReg, paramsReg, deformField] = elastixRegistration( x1, x2, MT1, MT2, params)

if nargin<2; error('elastixRegistration:: At least 2 images need to be provided');end
if nargin<3 || isempty(MT1); MT1=eye(4);end
if nargin<4 || isempty(MT2); MT2=eye(4);end
if nargin<5 || isempty(params); params='affine';end
[x1, x2, MT1, MT2] = parUnaFun({x1,x2,MT1,MT2},@gather);%Elastix does not support GPU functionality

%%% SET DIRECTORIES
elastixdir = elastixDir;
dirSave = temporaryResultsDir;
    if ~exist(dirSave,'dir');mkdir (dirSave);end
    if ~exist(strcat(dirSave,'/Input'),'dir');mkdir (strcat(dirSave,'/Input'));end
    if ~exist(strcat(dirSave,'/OutputEuler'),'dir');mkdir (strcat(dirSave,'/OutputEuler'));end
    if ~exist(strcat(dirSave,'/OutputAffine'),'dir');mkdir (strcat(dirSave,'/OutputAffine'));end
    if ~exist(strcat(dirSave,'/OutputBSpline'),'dir');mkdir (strcat(dirSave,'/OutputBSpline'));end
    if ~exist(strcat(dirSave,'/Transformation'),'dir');mkdir (strcat(dirSave,'/Transformation'));end
dirparams = transformationParametersDir;

regSort = 'single';%Hard coded 
if isstruct(params);regType = params.regType;else regType=params;end
assert(all(isreal(x1)) && all(isreal(x2)), 'elastixRegistration cannot be complex'); 

%%% REGRID ARRAYS
regridVolume = ~isequal(MT1, MT2);%To improve conditioning of the optimisation
useMask = size(x1,4)>1 ||  size(x2,4)>1;
if useMask; M1 = dynInd(x1,2,4);M2 = dynInd(x2,2,4);end
if regridVolume
    x2 = mapVolume(x2,x1,MT2, MT1);
    if useMask;M2 = mapVolume(M2,x1,MT2, MT1);end
end

MS1 = sqrt(sum(MT1(1:3,1:3).^2,1));
MS2 = sqrt(sum(MT2(1:3,1:3).^2,1));

%%% SAVE NII IMAGES
%%% First image
niftiIm=make_nii(x1,MS1(1:3));
niftiIm.hdr.hist.srow_x=MT1(1,:);niftiIm.hdr.hist.srow_y=MT1(2,:);niftiIm.hdr.hist.srow_z=MT1(3,:);niftiIm.hdr.hist.sform_code=1;         
save_nii(niftiIm,sprintf('%s/%s.nii',strcat(dirSave,'/Input'),'imageFixed'));

%%% Second image
niftiIm=make_nii(x2,MS2(1:3));
niftiIm.hdr.hist.srow_x=MT2(1,:);niftiIm.hdr.hist.srow_y=MT2(2,:);niftiIm.hdr.hist.srow_z=MT2(3,:);niftiIm.hdr.hist.sform_code=1;         
save_nii(niftiIm,sprintf('%s/%s.nii',strcat(dirSave,'/Input'),'imageMoving'));

if useMask
    %%% First mask
    niftiIm=make_nii(M1,MS1{n}(1:3));
    niftiIm.hdr.hist.srow_x=MT1(1,:);niftiIm.hdr.hist.srow_y=MT1(2,:);niftiIm.hdr.hist.srow_z=MT1(3,:);niftiIm.hdr.hist.sform_code=1;         
    save_nii(niftiIm,sprintf('%s/%s.nii',dirSave,'maskFixed'));

    %%% Second mask
    niftiIm=make_nii(M2,MS2{n}(1:3));
    niftiIm.hdr.hist.srow_x=MT2(1,:);niftiIm.hdr.hist.srow_y=MT2(2,:);niftiIm.hdr.hist.srow_z=MT2(3,:);niftiIm.hdr.hist.sform_code=1;         
    save_nii(niftiIm,sprintf('%s/%s.nii',dirSave,'maskMoving'));
end

%%% RUN ELASTIX REGISTRATION
%Parameters
if ~ischar(params)
    if strcmp(regSort,'single')
        write_parameters_single_euler ( params.HistBins, params.NumRes,params.NumIter,params.SpatSam );
        write_parameters_single_affine ( params.HistBins, params.NumRes,params.NumIter,params.SpatSam );
         if strcmp( regType,'bspline')
             write_parameters_single_bspline( params.HistBins, params.NumRes,params.NumIter,params.SpatSam,params.GridSpac );
         end
    else
        error('elastixRegistration:: Pair-wise registration not implemented yet');
    end
end

nameparameterEuler = strcat( dirparams,  'parameters_',regSort , '_euler.txt');
nameparameterAffine = strcat( dirparams,  'parameters_', regSort , '_affine.txt');
nameparameterBSpline = strcat( dirparams,  'parameters_', regSort , '_bspline.txt');

%Image 
fixed = fullfile(dirSave, 'Input', 'imageFixed.nii');
moving = fullfile(dirSave, 'Input','imageMoving.nii');

mask_fixed   = fullfile(dirSave, 'Input','maskFixed.nii');
mask_moving  = fullfile(dirSave, 'Input','maskMoving.nii');

outDir1 = fullfile(dirSave,'OutputEuler');
outDir2 = fullfile(dirSave,'OutputAffine'); 
outDir3 = fullfile(dirSave,'OutputBSpline');

if useMask 
    cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -fMask %s -mMask %s -out %s -p %s', fixed , moving, mask_fixed ,mask_moving , outDir1,nameparameterEuler));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
    fprintf('Euler registration finished.\n')
    cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -fMask %s -mMask %s -out %s -p %s -t0 %s',fixed, moving ,  mask_fixed, mask_moving, outDir2,nameparameterAffine, fullfile(outDir1,'TransformParameters.0.txt')));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
    fprintf('Affine registration finished.\n')
    if strcmp( regType,'bspline')
        cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -fMask %s -mMask %s -out %s -p %s -t0 %s',fixed, moving, mask_fixed, mask_moving, outDir3,nameparameterBSpline, fullfile(outDir2,'TransformParameters.0.txt')));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
        fprintf('Non-rigid registration finished.\n')  
    end
else
    cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -out %s -p %s', fixed , moving, outDir1,nameparameterEuler));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
    fprintf('Euler registration finished.\n')
    cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -out %s -p %s -t0 %s',fixed, moving,  outDir2,nameparameterAffine, fullfile(outDir1,'TransformParameters.0.txt')));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
    fprintf('Affine registration finished.\n')
    if strcmp( regType,'bspline')
        cmd = strcat(elastixdir,sprintf('elastix -f %s -m %s -out %s -p %s -t0 %s',fixed, moving, outDir3,nameparameterBSpline, fullfile(outDir2,'TransformParameters.0.txt')));[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Elastix');end
        fprintf('Non-rigid registration finished.\n')  
    end
end
if ~strcmp( regType,'bspline')
    registrationResults = fullfile(outDir2,'TransformParameters.0.txt');
else
    registrationResults = fullfile(outDir3,'TransformParameters.0.txt');
end

%%% RUN TRANSFORMIX
GETDEFORMFIELD = nargout>2;
paramsReg = readTXT(registrationResults);
outDir = fullfile(dirSave,'Transformation');
imageIn = moving;
if GETDEFORMFIELD
    cmd = strcat(elastixdir, sprintf('transformix -def all -in %s -out %s -tp %s',imageIn,outDir,registrationResults));
else
    cmd = strcat(elastixdir, sprintf('transformix -in %s -out %s -tp %s',imageIn,outDir,registrationResults));
end
[status,~]=dos(cmd);if status ~= 0; warning( 'elastixRegistration: Error occured during Transformix'); end
fprintf('Tranformation finished.\n')  

%%% LOAD TRANSFORMED IMAGE 
nii=load_untouch_nii(fullfile(outDir,'result.nii') );
x2Reg=nii.img;

xReg = cat(4, x1, x2Reg);
if useMask; xReg = cat(4, xReg, M1);end

if GETDEFORMFIELD
    nii=load_untouch_nii(fullfile(outDir,'deformationField.nii') );
    deformField=nii.img;
    %xReg = cat(4, xReg, deformField);
end

%%% CLEAN UP FIRECTORIES
fclose('all') ;
%rmdir(dirSave);
