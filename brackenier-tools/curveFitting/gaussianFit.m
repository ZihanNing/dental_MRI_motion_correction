
function [yFit, A, mu, sigma] = gaussianFit(y,x,dimFit)

% y = A*exp(-(x-mu)^2/(2*sigma^2) )
% A Fast, Accurate, and Separable Method for Fitting a Gaussian Function  Ibrahim Al-Nahhal, Octavia A. Dobre, Ertugrul Basar, Cecilia Moloney,  and Salama Ikki

if nargin <3 || isempty(dimFit);dimFit=numDims(y);end
if nargin <2 || isempty(x);x=1:size(y,dimFit);end

%%% RESHAPE
perm = 1:max(2,numDims(y));perm(dimFit)=[];perm=[perm dimFit];
y = permute(y,perm);
NYPerm = size(y);
y = resSub(y,1:(length(NYPerm)-1));
x = reshape(x,[1 size(y,2)]);
NG = size(y,2);%Number of elements along dimension on which to fit curve
NS = size(y,1);%Number of instances (separable fits)

%%% CREATE MATRICES TO SOLVE LINEAR SYSTEM X * Q = Y
%X matrix
X = zeros([3,3],'like', real(y));
for i=1:3
    for j=1:3
        if j>=i; X(i,j) = multDimSum(x.^(i+j-2));end%Avoid computing symmetric components and fill in later by enforcing symmetry
    end
end
%Enforce symmetry
X = (X+X.');
for i=1:3; X(i,i)=X(i,i)/2;end

%Y vector
Y = zeros([3,NS],'like', real(y));
for i=1:3; Y(i,:) = multDimSum(x.^(i-1).*log(y),2);end

%%% COMPUTE QUADRATIC PARAMETERS
Q = X\Y;
a = dynInd(Q,1,1);
b = dynInd(Q,2,1);
c = dynInd(Q,3,1);

%%% EXTRACT GAUSSIAN PARAMS
A = exp(a-b.^2./(4*c));
mu = -b./(2.*c);
sigma = sqrt(-1./(2.*c));

%%% FIT GAUSSIAN
yFit = A.'.*exp(-(x-mu.').^2./(2*sigma.'.^2) );

%%% RESHAPE INTO ORIGINAL SIZE
yFit = resSub(yFit, 1, NYPerm(1:(length(NYPerm)-1)));
yFit = ipermute(yFit,perm);
A = resSub(A(:), 1, NYPerm(1:(length(NYPerm)-1)));A = ipermute(A,perm);
mu = resSub(mu(:), 1, NYPerm(1:(length(NYPerm)-1)));mu = ipermute(mu,perm);
sigma = resSub(sigma(:), 1, NYPerm(1:(length(NYPerm)-1)));sigma = ipermute(sigma,perm);
