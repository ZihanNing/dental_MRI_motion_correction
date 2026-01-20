
function [x, p]=dephaseBasis(B, c, NX, TE, returnB0)

%DEPHASEBASIS   Computes the dephasing for a given set of basis functions and coefficients.
%   X=DEPHASEBASIS(B, C, NX, {TE}, {RETURNB0)
%   * B the set of basis functions ordered as [prod(NX) NCoef].
%   * c are the coefficients for corresponing motion states and basis functions.
%   * NX is the array size.
%   * {TE} is the echo time (in seconds).
%   * {RETURNB0} a flag to return the B0 field in Hz instead of the depahsed image.
%   ** x is a complex image with complex phase added.
%   ** p is the unwrapped phase.
%
%   Yannick Brackenier 2023-07-14

if nargin<4 || isempty(TE); TE=1;end
if nargin<5 || isempty(returnB0); returnB0=0;end

%%% GPU HANDLING
gpu = isa(B,'gpuArray');
if gpu && ~isa(c,'gpuArray'); c = gpuArray(c);end

%%% PREPARE ARRAYS
nT = size(c,5);
nC = size(c,6);
c = permute(c, [6 5 1:4]);%[NCoef NStates]
NX = NX(1:3);%Only spatial dimensions

%%% COMPUTE B0 FIELD
x = zeros(prod(NX) , nT, 'like', B); 
BlSzC = 5; %For coefficients
for i= 1:BlSzC:nC
    vI=i:min(i+BlSzC-1,nC);
    x = x + dynInd(B, vI,2) * dynInd(c,vI,1);
end
x=reshape(x , [NX 1 nT] );%Hz
if returnB0; return; end

%%% COMPUTE PHASE
if nargout>1; p = 2*pi*TE * x;end %Unwrapped phase [radians]
x=exp(+1i *2*pi*TE * x);%Wrapped phase

