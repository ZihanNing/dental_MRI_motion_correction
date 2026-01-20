
function on = oneStruct(n, type)

%ONESTRUCT creates a cell array where element i consist of an array ones([1 i]). 
%
%   [ON]=ONESTRUCT({N})
%   * {N} is the length of the cell array.
%   ** ON the cell array.
%
%   Yannick Brackenier 2022-01-30

if nargin < 1 || isempty(n); n=5;end
if nargin < 2 || isempty(type); type='single';end
    
on = cell(1,n);
for i = 1:n
    on{i} = ones([1 i],type);
end