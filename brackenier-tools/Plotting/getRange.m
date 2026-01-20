
function [rrMagn, rrPhase, rrReal, rrImag] = getRange(x, excludeZero, getRealImag)

%GETRANGE finds the range of the elements in an array .
%   [RRMAGN,RRPHASE,RRREAL,RRIMAG]=INRANGE(X, {EXCLUDEZERO}, {GETREALIMAG})
%   * X is the array.
%   * {EXCLUDEZERO} a flag to not include zeros when defining the range.
%   * {GETREALIMAG} a flag to return the real and imaginary ranges as th first 2 output arguments.
%   ** RRMAGN is the range of the magnitude of X. This will be RRREAL if GETREALIMAG==1.
%   ** RRPHASE is the range of the phase of X. This will be RRIMAG if GETREALIMAG==1.
%   ** RRREAL is the range of the real component of X.
%   ** RRIMAG is the range of the imaginary component of X.
%
%   Yannick Brackenier 2022-08-30

if nargin<2 || isempty(excludeZero);excludeZero=0;end
if nargin<3 || isempty(getRealImag);getRealImag=0;end

%%% Exlude zeros
if excludeZero; M = x==0;x(M)=[];end

%%% Flatten array and set default in case it is empty
x=x(:);
if isempty(x)%Add elements to get a default range
    x = [0 inf];%[0 inf] for the magnitude.
    x = [x exp(1i*pi) exp(-1i*pi)];%[-pi pi] for the phase.
    x = [x -inf inf];%[-inf inf] for the real part.
    x = [x 1i*inf  -1i*inf];%[-inf inf] for the imaginary part.
end

%%% Compute ranges for all the data types
rrMagn = [ min(abs(x))  max(abs(x))];
rrPhase = [ min(angle(x))  max(angle(x))];

rrReal = [ min(real(x)) max(real(x))];
rrImag = [ min(imag(x))  max(imag(x))];

%%% Adjust output if needed
if getRealImag
    rrMagn = rrReal;
    rrPhase = rrImag;
end