
function [] = plotBox(x, y, xStart, xStop, yStart, yStop)

if nargin<3||isempty(xStart);xStart=[];end
if nargin<4||isempty(xStop);xStart=[];end
if nargin<5||isempty(yStart);xStart=[];end
if nargin<6||isempty(yStop);xStart=[];end

%Set parameters
factBox=1;
if ~isempty(x);dx=x(2)-x(1);end

%Set x boundaries
if isempty(xStart);xStart=x(1)-dx/2;end
if isempty(xStop);xStop=x(end)+dx/2;end

boxx=[xStart  xStart xStop xStop]; %+.5*dx because otherwise neighbourhing boxes have a small white space in between- only for visual purposes

%Set y boundaries
if isempty(yStart);yStart=gather(min(y(:))*factBox);end
if isempty(yStop);yStop=gather(max(y(:))*factBox);end
boxy = zerosL(boxx);
boxy([1,4])= yStart;
boxy([2,3])= yStop;

%Plot
patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none');%'EdgeColor',[.2 0 0])% 'LineWidth',0.00001,   

