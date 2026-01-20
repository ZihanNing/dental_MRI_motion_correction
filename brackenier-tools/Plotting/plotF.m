

function [] = plotF(x,y,xlab,ylab,TIT)

if nargin<1 || isempty(x);x=[];end
if nargin<3 || isempty(xlab);xlab='';end
if nargin<4 || isempty(ylab);ylab='';end
if nargin<5 || isempty(TIT);TIT='';end

%figure('color','w');

N =size(y);N(end+1:2)=1;
if N(2)==1; y=y(:).';end
if isempty(x);x=1:size(y,2);end

plot(x,y);
grid on
axis tight
xlabel(xlab);
ylabel(ylab);

title(TIT);