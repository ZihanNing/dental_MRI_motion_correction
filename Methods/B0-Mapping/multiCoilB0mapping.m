

function [B0] = multiCoilB0mapping(x, TE, MS, useUnwrapping)

%MULTICOILB0MAPPING   Performs B0estimation from multi-TE raw coil images.
%   * X is the array with different echoes stored in the 4th dimension and spatial dimensions in dimensions 1:3
%   * TE is a list with the corresponding echo times in ms
%   ** B0 is the estimated B0 in Hz.
%
%   Yannick Brackenier 2022-10-05

if nargin<3 || isempty(MS);MS = [1 1 1];end
if nargin<4 || isempty(useUnwrapping);useUnwrapping = 0;end

N = size(x);N=N(1:3);
NS = size(x,4);
NT = size(x,5);

%%% Extract first 2 echoes
if NT>2
    x = dynInd(x,1:2,5);%Use only first 2 echoes
    NT = size(x,5);
    TE = dynInd(TE(:),1:2,1);
    warning('multiCoilB0mapping:: Extracted the first 2 echoes.');
end
    
%%% Control input
assert(NT == length(TE(:)),'multiCoilB0mapping:: Number of input images does not correspond with echo times provided.');
assert(all(diff(TE)>0),'multiCoilB0mapping:: Echo times must be increasing.');

%%% Prepare 2 echo images
x1 = dynInd(x,1,5);
x2 = dynInd(x,2,5);

%%% Compute B0
dTE = TE(2)-TE(1);%ms
dTE = dTE/1000;%s
B0 = x2.*conj(sign(x1));%By taking only the sign of the second TE, this is weigted by the first coil sensitivity
if useUnwrapping
    for i=1:NS %For every coil individually
        fprintf('Unwrapping phase coil %d\n',i);
        %phaseOrig = angle(dynInd(B0,i,4));
        B0 = dynInd(B0,i,4,CNCGUnwrapping(dynInd(B0,i,4), MS, 'Magnitude','LSIt'));%MagnitudeGradient4
        %plotND([],cat(4, phaseOrig,B0, angle(exp(1i*B0))),[],[],0,{[],2},[],{'Original phase';'Unwrapped';'Unwrapped prediced phase'});
    end
    B0 = multDimMea(B0,4) ;%Average across the channel direction
else
    B0 = multDimMea(B0,4) ;%Average across the channel direction
    B0 = angle(B0);     
end
B0 = (1/dTE/2/pi)*B0; %Convert phase to Hz
