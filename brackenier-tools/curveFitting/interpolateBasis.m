
function [fInt, c, fRes] = interpolateBasis (f, B, M, domain)

%INTERPOLATEBASIS interpolates an array f with a given basis B using a weighted Least-Squares (wLS). f = B*c
%   [FINT, C, RES]=INTERPOLATEBASIS(F,B,{M},{DOMAIN})
%   * F is the array to interpolate.
%   * B is the basis in which to interpolate f. Different basis functions are stored in the columns (or the end+1 dimensions of f).
%   * {M} is the mask to use for a weighted LS.
%   * {DOMAIN} is the domain to use. This reduces memory load. Equivalent to binary wLS, although the latter would not reduce memory load.
%   ** FINT is the interpolated array.
%   ** {C} are the coefficients of the respective basis functions in B.
%   ** {FRES} are the residuals between the original array f and the interpolated one fInt.
%
%   Yannick Brackenier

if nargin <2 || isempty(B); error('interpolate_basis: No basis provided.');end
if nargin <3 || isempty(M); M = onesL(f);end
if nargin <4 || isempty(domain); domain = onesL(f);end

N = size(f);ND = ndims(f);
if size(B,ND+1)>1; B = resSub(B, 1:ND);end

%%% Only retain domain to save memory
f(domain==0) = [];
M(domain==0) = [];
B(domain==0, :) = [];

%%% Reshape to column vector
f = resSub(f, 1:ND);
M = resSub(M, 1:ND);

%%% Weight basis and f to have weighted LS
Bw = bsxfun( @times, sqrt(M), B);
fw =  sqrt(M).*f;

%%% LS interpolation
c = Bw\fw;
fInt = B*c;

%%% Reshape and add domain
temp = fInt; fInt = zeros(N,'like',fInt) ; 
fInt(domain~=0) = temp; clear temp;

%%% Report coefficients and residuals
if nargout>2
    temp = f; f = zeros(N,'like',f) ; 
    f(domain~=0) = temp; clear temp;

    fRes = f - fInt;
end
