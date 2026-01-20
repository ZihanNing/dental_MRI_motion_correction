function x=multDimStd(x,dim)

%MULTDIMSTD   Takes the standard deviation of the elements of a multidimensional array 
%along a set of dimensions
%   X=MULTDIMSTD(X,{DIM})
%   * X is an array
%   * {DIM} are the dimensions over which to take the std of the elements 
%   of the array, defaults to all
%   ** X is the contracted array
%

if nargin<2 || isempty(dim);dim=1:numDims(x);end

for n=1:length(dim);x=std(x,[],dim(n));end
