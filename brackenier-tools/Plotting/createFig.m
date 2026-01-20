
function [h] = createFig(numFig,figName)

if nargin<1 || isempty(numFig);numFig=[];end
if nargin<2 || isempty(figName);figName='';end

%%% Test if we need to plot
if ~newFigFlag(numFig);h=gcf;return;end

%%% Create figure
if isempty(numFig);h = figure;else;h = figure(numFig);end

%%% Set parameters
clf; 
%set(h,'Name',figName);
set(h,'color', 'w');
set(h,'Position',get(0,'ScreenSize'));

%Pause - seemed to help with plotting the CURRENT figure
figure(h);
%pause(.01);
