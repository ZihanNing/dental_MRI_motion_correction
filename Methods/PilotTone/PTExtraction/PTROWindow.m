
function [hh] = PTROWindow(N, factor, shift, gibbs)

if nargin < 3 || isempty(shift); shift=0;end
if nargin < 4 || isempty(gibbs); gibbs=.5;end

dimShift = find(N>1);

h=fftshift(buildFilter(N,'tukey',factor,0,gibbs));
startId = find(h(:)~=0);
L = length(startId);
startId=startId(1);

hh = circshift(h,-(startId-1));

shift=max(0,round(shift*N(dimShift) - L/2));
hh = circshift(hh,min(shift, N(dimShift)-L));


% figure(); 
% plot(hh(:));
% axis tight



