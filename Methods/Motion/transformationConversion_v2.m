
function [TOut] = transformationConversion_v2(TIn, MTIn, MTOut, NIn, NOut, centGridIn, centGridOut)

%TRANSFORMATIONCONVERSION_V2 Converts transformation parameters between geometries.
%   [TOUT]=TRANSFORMATIONCONVERSION(TIN, MTIN, MTOUT, NIN, NOUT, CENTGRIDIN, CENTGRIDOUT)
%   * TIN the transformation parameters.
%   * {MTIN} the S-FORM of the input parameters.
%   * {MTOUT} the S-FORM of the output parameters.
%   * {NIN} the array size of the input parameters.
%   * {NOUT} the array size of the output parameters.
%   * {CENTGRIDIN} the center of rotation for the input array.
%   * {CENTGRIDOUT} the center of rotation for the output array.
%   ** TOUT the converted transformation parameters.
%
%   Yannick Brackenier 2023-08-29

if nargin < 2 || isempty(MTIn); MTIn = eye(4); end
if nargin < 3 || isempty(MTOut); MTOut = eye(4); end
if nargin < 4 || isempty(NIn); NIn = zeros([1 3]); end
if nargin < 5 || isempty(NOut); NOut = zeros([1 3]); end
if nargin < 6 || isempty(centGridIn); centGridIn = ceil((NIn(1:3)+1)/2); end
if nargin < 7 || isempty(centGridOut); centGridOut = ceil((NOut(1:3)+1)/2); end

%%% INITIALISE
dimM = numDims(TIn);
dimS = dimM-1;
NStates = size(TIn,dimS);

if isequal(MTIn,MTOut) && isequal(NIn,NOut) && isequal(centGridIn,centGridOut)
    fprintf('transformationConversion_v2:: The geometries are the same, so the input parameters are returned.\n')
    TOut=TIn;
    return
end

%%% RUN OVER MOTION STATES
TOut = zerosL(TIn);
for i = 1:NStates
    
    %%% EXTRACT
    Ti = dynInd(TIn,i,dimS);
    
    %%% TRANSFORMATION TO AFFINE MATRIX
    R = convertT(Ti,'T2R'); %transformation matrix for sincRigidTransform.m
    
    %%% FROM INPUT TO RAS
    R = Rt(centGridIn) * R * Rt(-centGridIn) ; 
    RNew = (MTIn * R)/MTIn; %same as RNew = MTIn * R * inv(MTIn);
      
    %%% FROM RAS TO OUTPUT
    RNew = MTOut\(RNew * MTOut); %same as inv(MTOut) * R * MTOut;
    RNew = Rt(-centGridOut) * RNew * Rt(centGridOut) ; 
    
    %%% AFFINE MATRIX TO TRANSFORMATION
    TiNew = convertT(RNew,'R2T');
    
    %%% STORE
    TOut = dynInd(TOut, i, dimS, TiNew);
end

