
function [h,p] = histPlot(data, xLabel, yLabel, legendNames, Title, numFig, FontSize, xRange, yRange, binWidth, co)

%%% HISTPLOT is creates a histogram plot of multiple distributions.
%
%   [H, P] = HISTPLOT(DATA, XLABEL, YLABEL, LEGENDNAMES, TITLE, NUMFIG, FONTSIZE, XRANGE, YRANGE, BINWIDTH,CO)
%   * DATA contains the array with size [NGroups NSamples].
%   * {XLABEL} the label on the x-axis.
%   * {YLABEL} the label on the y-axis.
%   * {LEGENDNAMES} contains names of the legend.
%   * {TITLE} the title string.
%   * {NUMFIG} the figure number. If negative, no new figure will be created.
%   * {FONTSIZE} the font size to be used.
%   * {XRANGE} the range of the x-axis.
%   * {YRANGE} the range of the y-axis.
%   * {BINWIDTH} the width of the bins of the histogram. If negative, it is the number of bins.
%   * {CO} the colors of the histrograms as a matrix of size [NColors 3].
%   ** H the figure handle
%   ** P the individual handles of the histogram line and histrogram filling.
%
%   Yannick Brackenier 2022-11-30

if nargin < 2 || isempty(xLabel); xLabel = '';end
if nargin < 3 || isempty(yLabel); yLabel = '';end
if nargin < 4 || isempty(legendNames); legendNames = [];end
if nargin < 5 || isempty(Title); Title = '';end
if nargin < 6 || isempty(numFig); numFig = 999;end
if nargin < 7 || isempty(FontSize); FontSize = 16;end
if nargin < 8 || isempty(xRange); xRange=[];end
if nargin < 9 || isempty(yRange); yRange=[];end
if nargin < 10 || isempty(binWidth); binWidth = -1000;end%If negative, it's the number of bins
if nargin < 11 || isempty(co); co = defaultColors();end

if binWidth<0
    rangeData  = multDimMax(data) - multDimMin(data);
    rangeData = rangeData + eps;
    binWidth = rangeData/abs(binWidth);
end

numGroups = size(data,1);
numSamp = size(data,2);
 
if numGroups>size(co,1); co = repmat(co,[ceil(numGroups/size(co,1)) 1]);end
if isempty(legendNames);legendNames=cell(1,numGroups);legendNames(:) = {''};end

%%% Figure design
LineWidth = 8;
alphaLine = 1;
alphaFill = .4;
ModXAxis = .4*FontSize;
ModYAxis = .6*FontSize;
ModLeg = - 2;

%%% Create figure
if newFigFlag(numFig)
    h = figure(numFig);clf;
    set(h,'Name','','color', 'w');
    set(h,'color','w','Position',get(0,'ScreenSize'));
end

p= {};
listPlotsBar = zeros (1, numGroups);
listPlotsStair= zeros (1, numGroups);

%%% Add histograms
hold on;
for i =1:numGroups
        h1 = histogram(data(i,:),'DisplayStyle','bar','BinWidth',binWidth,'LineWidth',eps/100000,'FaceColor',co(i,:),'EdgeColor',co(i,:),'LineStyle','none','EdgeAlpha',alphaLine,'FaceAlpha',alphaFill);
        h2 = histogram(data(i,:),'DisplayStyle','stair','BinWidth',binWidth,'LineWidth',LineWidth,'EdgeAlpha',alphaLine,'EdgeColor',co(i,:));
        p{i} = {h1,h2};
        listPlotsBar(i) = h1;
        listPlotsStair(i) = h2;
end
hold off

%%% Title
title(Title, 'interpreter','latex','FontSize', FontSize)

%%% Axes
set(gca, 'FontSize',FontSize, 'TickLabelInterpreter','latex');
ylabel(yLabel, 'interpreter', 'latex', 'FontSize', FontSize + ModYAxis)
xlabel(xLabel, 'interpreter', 'latex', 'FontSize', FontSize + ModXAxis)

%%% Legend
if any( cellfun(@(x) ~isempty(x),legendNames) )
    %AX=legend(listPlotsBar,legendNames,'Location','NorthEast','FontSize',FontSizeLegend,'interpreter','latex');  
    AX=legend(listPlotsStair,legendNames,'Location','NorthEast','FontSize',FontSize+ModLeg,'interpreter','latex'); 
    legend boxoff
    set(AX,'Position',get(AX,'Position')+[0 0 0 0])%To the right and to the top   
end

%%% Set range limits
if isempty(xRange)
    xRange=gather( [min(multDimMin(data))  max(multDimMax(data)) ]    );
    xRange = xRange + diff(xRange)/20*[-1 1];
    xRange = xRange.*[.999 1.001];
end
if isempty(yRange);yRange=gather(1.1*[0 max(h1.BinCounts)]);end%getRange(data,[],1);

xlim([xRange(1) xRange(2)])  
if ~isempty(yRange);ylim([yRange(1) yRange(2)])  ;end

