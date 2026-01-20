
function [flag, elFalse, elTrue] = isequalExt(x,y,buffer,makeAssert)

%ISEQUALEXT checks if two arrays are identical, taking into account machine precision.
%   [FLAG,ELFALST,ELTRUE]=ISEQUALEXT(X, Y, {BUFFER},{MAKEASSERT})
%   * X first array.
%   * Y second array.
%   * {BUFFER} a multiplicative buffer to the machine precision when comparing residuals.
%   * {MAKEASSERT} whther to make it an assert.
%   ** FLAG a logical flag to indicate the arrays are identical.
%   ** ELFALST the number of elements that do not fall within machine precision.
%   ** ELTRUE the number of elements that fall within machine precision.
%
%   Yannick Brackenier 2023-07-20

if nargin<3 || isempty(buffer);buffer=1;end
if nargin<4 || isempty(makeAssert);makeAssert=0;end

%%% CLASS HANDLING
classX = class(x);
classY = class(y);
epsClass = eps(classX);
%Cast if needed
if ~strcmp(classX, classY)
    warning('isequalExt:: x and y must have same class. y converted from %s to %s. Avoid this for computational overhead.',classY, classX);
    y = cast(y, classX);
end

%%% COMPUTE RELATIVE ERROR
err = abs(x-y)./ (multDimMea(abs(y)) + epsClass);
err = err(:);

%%% CHECK THRESHOLD
elFalse = multDimSum(err > buffer*epsClass);
elTrue = multDimSum(err <= buffer*epsClass);

flag = elFalse==0;
if makeAssert; assert(flag, sprintf('isequalExt:: %d/%d are not within machine precision.',elFalse,elTrue+elFalse));end

