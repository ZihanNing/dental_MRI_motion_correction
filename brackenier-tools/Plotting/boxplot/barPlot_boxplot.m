
function [h,p] = barPlot_boxplot(data, xLabel, xTickLabel, yLabel, legendNames, Title, numFig, FontSize, yRange, co)

%%% BARPLOT creates a barplot of the data 
%
%   [H, P] = HISTPLOT(DATA, XLABEL, XTICKLABEL,YLABEL, LEGENDNAMES, TITLE, NUMFIG, FONTSIZE, YRANGE,CO)
%   * DATA contains the array with size [experiments(blocks) groups 2] where the last dimension contains the mean and standard deviation. In case
%          size(DATA,3)>2, it is the number of samples and the mean and standard deviation will be computed along this dimension.
%   * {XLABEL} the label on the x-axis.
%   * {XTICKLABEL} the ticks on the x-axis.
%   * {YLABEL} the label on the y-axis.
%   * {LEGENDNAMES} contains names of the legend.
%   * {TITLE} the title string.
%   * {NUMFIG} the figure number.
%   * {FONTSIZE} the font size to be used.
%   * {YRANGE} the range of the y-axis.
%   * {CO} the colors of the histrograms as a matrix of size [NColors 3].
%   ** H the figure handle
%
%   Yannick Brackenier 2022-11-30

if nargin < 1 || length(size(data))<2; error('barPlot:: Data provided not suitable.');end
numExp = size(data,1);
numGroups = size(data,2);

%if size(data,3)>2;data = cat(3, multDimMea(data,3) , multDimStd(data,3) ) ; end %Compute mean and standard devitation
if size(data,3)<2;data = cat(3, data , zeros(size(data))) ; end %Append standard devitation if not provided

if nargin < 2 || isempty(xLabel); xLabel = '';end
if nargin < 3 || isempty(xTickLabel); xTickLabel = cell(1,numExp);xTickLabel(:) = {''};end
if nargin < 4 || isempty(yLabel); yLabel = '';end
if nargin < 5 || isempty(legendNames); legendNames = cell(1,numGroups);legendNames(:) = {''};end
if nargin < 6 || isempty(Title); Title = '';end
if nargin < 7 || isempty(numFig); numFig = [];end
if nargin < 8 || isempty(FontSize); FontSize = 16;end
if nargin < 9 || isempty(yRange); yRange=[];end
if nargin < 10 || isempty(co); co = defaultColors(numGroups);end
assert(size(co,1)>=numGroups,'barPlot:: Colors provided not compatible with data.');
FontSizeLegend = FontSize - 2;

idxList = 0:800; % Large enough if numbers of levels increases
barWidth = 1;% Width of an individual bar
expSpacing = 1;% Distance between center of neighbouring experiment blocks
groupInterleave = numGroups*barWidth+ expSpacing;% Distance between edges of neighbouring groups

h = createFig(numFig);hold on
p= {};
listPlots = zeros (1, numGroups);

for i = 1:numGroups
       %Prep data
        mean = data(:,i,1);
        stdv = data(:,i,2); 
        
        %Barplot
        idxPlot =  idxList(1:numExp)*(numGroups*barWidth + expSpacing) +i; %First element of each group
        if numExp==1
           width=1;
           width = 0.7
        else
           width = 1/(numGroups+1);
        end
        for jj=1:size(data,1)
            [p{i}] = plot_ind_box_yb(idxPlot(jj),vec(data(jj,i,:)),width, co(i,:))
            %boxplot(vec(data(jj,i,:)), 'Positions',idxPlot(jj))
        end
        
        %p{i} = bar(idxPlot,mean,width, 'FaceColor',co(i,:));
        %listPlots(i) = p{i};
        
        %Errorplot
        %errorbar(idxPlot, mean,stdv,'LineStyle','none','Color',[0 0 0],'Linewidth',2);
end
hold off

%%% Title
title( Title, 'interpreter','latex','FontSize', FontSize)

%%% Axes
xTicks = idxList(1:numExp)*(numGroups*barWidth + expSpacing) + (numGroups*barWidth+1)/2 ;

set(gca,'XTick',xTicks,'XTickLabel',xTickLabel);
set(gca, 'FontSize',FontSize, 'TickLabelInterpreter','latex');
ylabel(yLabel, 'interpreter', 'latex', 'FontSize', FontSize)
xlabel(xLabel, 'interpreter', 'latex', 'FontSize', FontSize)

%%% Legend
if any( cellfun(@(x) ~isempty(x),legendNames) )
    AX=legend(listPlots,legendNames,'Location','NorthEast','FontSize',FontSizeLegend,'interpreter','latex');    
    %LEG = findobj(AX,'type','text');
    set(AX,'Position',get(AX,'Position')+[0 0 0 0])%To the right and to the top         
end
%%% Set ranges
rangeLow = data(:,:,1)-data(:,:,2);
rangeUp = data(:,:,1)+data(:,:,2);
if isempty(yRange)
   range = gather(1.1* [ min(min(rangeLow(:)),0)  max(rangeUp(:)) ]    );
else
   range = yRange; 
end

xlim([0 numExp*(numGroups*barWidth + expSpacing) ])
ylim([range(1) range(2)])  

set(h,'color','w','Position',get(0,'ScreenSize'));
