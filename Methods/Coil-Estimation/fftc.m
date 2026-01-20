
function res = fftc(x,dim)
% res = fftc(x,dim)
if nargin < 2
    dim = 2;
end

x = ifftshift(x,dim);
x = fft(x,[],dim);
x = fftshift(x,dim);
res = 1/sqrt(size(x,dim))*x;