
function res = ifftc(x,dim)

%res = ifftc(x,dim)

x = ifftshift(x,dim);
x = ifft(x,[],dim);
x = fftshift(x,dim);
res = sqrt(size(x,dim))*x;

