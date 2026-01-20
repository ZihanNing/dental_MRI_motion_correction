

function [B, coefNamePaper,coefNameTerra ] = shimBasis(N, MT, shimOrder, includeCorr)

%SHIMBASIS creates the set of basis functions for the active shimming used on the SIEMENS 7T scanner
%   [B,COEFIDX]=SHIMBASIS(X,ORDIM,{M})
%   * N is the array size for which to get the basis.
%   * MT if the orientation information in the NIFTI convention
%   * {SHIMORDER} is the shim order used. Defaults to 3
%   * {INCLUDECORR} Indicates whether to apply the cross term correction. Defaults to 0
%   ** B is the set of basis functions, stored in the 4th dimension (dimStore)
%   ** COEFNAMEPAPER are the shim terms, as reported in "Chang P, Nassirpour S, Henning A. Modeling Real Shim Fields for Very High Degree (and Order) B0 Shimming of the Human Brain at 9.4 T. Magnetic Resonance in Medicine 79:529–540 (2018)"
%   ** COEFNAMETERRA are the shim terms, as named in the Siemened MAGNETOM Terra scanner console
%
%   Yannick Brackenier 2023-07-04

gpu=(gpuDeviceCount>0 && ~blockGPU);

if nargin <1 || isempty(N); error('shimBasis:: Array size must be provided.');end
if nargin <2 || isempty(MT); error('shimBasis:: Array orientation must be provided.');end
if nargin <3 || isempty(shimOrder); shimOrder=3; warning('shimBasis:: Shim order not provided. Default value = 3 (Siemens 7T system).');end
if nargin <4 || isempty(includeCorr); includeCorr=0;end

%%% CONVERT NIFTI HEADER MT TO XYZ COORDINATES INSTEAD OF RAS
T_XYZ2RAS = diag([-1 1 -1 1]);
MT_XYZ = T_XYZ2RAS*MT;%Because around moving axis

%%% GENERATE XYZ GRIDS
rdGrid=generateGrid(N,gpu,N,ones(1,3));
[dGrid{1},dGrid{2},dGrid{3}]=ndgrid(rdGrid{1}(:),rdGrid{2}(:),rdGrid{3}(:));dGrid{4}=dGrid{3};dGrid{4}(:)=1;
Grid_ijk=vertcat(dGrid{1}(:)',dGrid{2}(:)',dGrid{3}(:)',dGrid{4}(:)');dGrid{4}=[];%NIFTI coordinates

Grid_xyz=MT_XYZ*Grid_ijk; Grid_ijk=[];%World coordinates
Grid_xyz = resSub(Grid_xyz.', 1,N);
X = dynInd(Grid_xyz,1,4)*1e-3;%In meters
Y = dynInd(Grid_xyz,2,4)*1e-3;%In meters
Z = dynInd(Grid_xyz,3,4)*1e-3;%In meters
Grid_xyz=[];%save memory

%%% CREATE SHIMMING SOLID HARMONICS
B=[]; b=0;
dimStore=4;
R = sqrt(X.^2 + Y.^2);

    % ZEROTH ORDER
    if shimOrder>=0
        B = cat(dimStore,B,ones(N));b=b+1; coefNamePaper{b}='Z0';coefNameTerra{b}='F/A00';
    end
    % FIRST ORDER
    if shimOrder>=1
        B = cat(dimStore,B,X);b=b+1; coefNamePaper{b}='X';coefNameTerra{b}='X/A11';
        B = cat(dimStore,B,Y);b=b+1; coefNamePaper{b}='Y';coefNameTerra{b}='Y/B11';
        B = cat(dimStore,B,Z);b=b+1; coefNamePaper{b}='Z';coefNameTerra{b}='Z/A10';
    end
    % SECOND ORDER
    if shimOrder>=2
        B = cat(dimStore,B,Z.^2-0.5.*(X.^2+Y.^2));b=b+1; coefNamePaper{b}='Z2';coefNameTerra{b}='Z^2/A20';
        B = cat(dimStore,B,Z.*X);b=b+1; coefNamePaper{b}='ZX';coefNameTerra{b}='ZX/A21';
        B = cat(dimStore,B,Z.*Y);b=b+1; coefNamePaper{b}='ZY';coefNameTerra{b}='ZY/B21';
        B = cat(dimStore,B,X.^2-Y.^2);b=b+1; coefNamePaper{b}='C2';coefNameTerra{b}='X^2-Y^2/A22';
        B = cat(dimStore,B,2*X.*Y);b=b+1; coefNamePaper{b}='S2';coefNameTerra{b}='XY/B22';
    end

    % THIRD ORDER
    if shimOrder>=3
        B = cat(dimStore,B,Z.*(Z.^2-3/2.*(X.^2+Y.^2)));b=b+1; coefNamePaper{b}='Z3';coefNameTerra{b}='Z^3/A30';
        B = cat(dimStore,B,X.*(Z.^2-1/4.*(X.^2+Y.^2)));b=b+1; coefNamePaper{b}='Z2X';coefNameTerra{b}='Z^2X/A31';
        B = cat(dimStore,B,Y.*(Z.^2-1/4.*(X.^2+Y.^2)));b=b+1; coefNamePaper{b}='Z2Y';coefNameTerra{b}='Z^2Y/B31';
        B = cat(dimStore,B,Z.*(X.^2-Y.^2));b=b+1; coefNamePaper{b}='ZC2';coefNameTerra{b}='Z(X^2-Y^2)/A32';
        %Next ones are the other 3rd order solid harmonics not included in the Siemens console
        B = cat(dimStore,B,2*Z.*X.*Y);b=b+1; coefNamePaper{b}='ZS2';coefNameTerra{b}='ZXY/None';
        B = cat(dimStore,B,X.*(X.^2-3*Y.^2));b=b+1; coefNamePaper{b}='C3';coefNameTerra{b}='X3/None';
        B = cat(dimStore,B,Y.*(3*X.^2-Y.^2));b=b+1; coefNamePaper{b}='S3';coefNameTerra{b}='Y3/None';       
    end
    
    % FOURTH ORDER
    if shimOrder>=4
        B = cat(dimStore,B, Z.^4-3*Z.^2.*R.^2 + 3/8.*R.^2 );b=b+1; coefNamePaper{b}='Z4';coefNameTerra{b}='Z4/None';
        B = cat(dimStore,B, Z.*X.*(Z.^2-3/4.*R.^2) );b=b+1; coefNamePaper{b}='Z3X';coefNameTerra{b}='Z3X/None';
        B = cat(dimStore,B, Z.*Y.*(Z.^2-3/4.*R.^2) );b=b+1; coefNamePaper{b}='Z3Y';coefNameTerra{b}='Z3Y/None';
        B = cat(dimStore,B, (X.^2-Y.^2).*(Z.^2-1/6.*R.^2) );b=b+1; coefNamePaper{b}='Z2C2';coefNameTerra{b}='Z2C2/None';
        B = cat(dimStore,B, 2*X.*Y.*(Z.^2-1/6.*R.^2) );b=b+1; coefNamePaper{b}='Z2S2';coefNameTerra{b}='Z2S2/None';
        B = cat(dimStore,B, Z.*X.*(X.^2-3.*Y.^2) );b=b+1; coefNamePaper{b}='ZC3';coefNameTerra{b}='ZC3/None';
        B = cat(dimStore,B, Z.*Y.*(3*X.^2-Y.^2) );b=b+1; coefNamePaper{b}='ZS3';coefNameTerra{b}='ZS3/None';
        B = cat(dimStore,B, (X.^2-Y.^2).^2 - 4.* X.^2.*Y.^2 );b=b+1; coefNamePaper{b}='C4';coefNameTerra{b}='X4/None';
        B = cat(dimStore,B, 4*X.*Y.*(X.^2-Y.^2) );b=b+1; coefNamePaper{b}='S4';coefNameTerra{b}='Y4/None';
    end    
      
%plot_inner_prod(B, 1, [0 10^(-3)]);

if includeCorr
   warning('shimBasis:: Correction not implemented in basis functions.')
   %Load cross terms and add to the basis
   %Can also apply correction to coefficients --> Done at this point
end

%%% CHECK FOR ERRORS
if shimOrder>=4; warning('shimBasis:: 4th order solid harmonics not implemented. Not added to basis.');end
if any(isnan(B(:))) && ~isequal(N,[1 1 1]); error('shimBasis::Basis function contain NaN');end
assert(all(isreal(B(:))),'shimBasis::Basis must be real');

end