
function [h, coilSortInfToSup, ccm] = displayCoilGeom(S,MT,y,rr,titleName, printValue)

%DISPLAYCOILGEOM displays the coil geometry based on the sensitivity profiles. 
%   [H, COILSORTINFTOSUP, CCM]=DISPLAYCOILGEOM(S,MT,{Y},{RR})
%   * S are the sensitivity profiles (stored in the 4th dimension).
%   * MT is the s-form of S.
%   * {Y} are the values to be plotted as color coding on each coil element.
%   * {RR} is range to plot for Y.
%   ** H is the figure handle.
%   ** COILSORTINFTOSUP is a list with the coils listed from inferior to superior direction.
%   ** CCM is the coil compression matrix associated to the array COILSORTINFTOSUP.
%
%   Yannick Brackenier 2023-03-21

if nargin<3 || isempty(y);y=[];end
if nargin<4 || isempty(rr);rr=[];end
if nargin<5 || isempty(titleName);titleName='';end
if nargin<6 || isempty(printValue);printValue=0;end

%%% COMPUTE CERTROILS OF COIL ELEMENTS
[RAS] = coilCentroids(S,MT);

%%% SORT IN INFERIOR-SUPERIOR DIRECTION
Sup = RAS(3,:);
[~,coilSortInfToSup] = sort(Sup);
ccm = eye(size(S,4));
ccm = ccm( flip(coilSortInfToSup) ,:);

%%% ASSIGN DATA FOR COLOR CODING
if ~isempty(y)
    if isempty(rr)
        rangeAssign = [min(y(:)),max(y(:))];
    else
        rangeAssign = rr;
    end
    [color, cMap] = mapToColormap(y(:), 'jet', rangeAssign );
else
    color = [0 0 1];
end

%%% PLOT
%Figure properties
FontSize = 20;

%Create figures
h = figure('color','w');
set(h,'Position',get(0,'screenSize'));

%Plot string with coil indices
offset = [10 10 10]/2;%To not overlap with markers
if printValue
    N=0;
    str = string(round(y,N));
else
    str = string(1:size(S,4));
end
textscatter3(RAS(3,:)+offset(3),RAS(1,:)+offset(1), RAS(2,:)+offset(2),str,'MarkerSize',20,'TextDensityPercentage',100);

%Plot the color data (markers)
hold on
offset = [0 0 0];
textscatter3(RAS(3,:)+offset(3),RAS(1,:)+offset(1), RAS(2,:)+offset(2),str,'MarkerSize',20,'TextDensityPercentage',0,'ColorData',color);
if ~isempty(y)
    colorbar; colormap(cMap);%If you do this after the text, it won't be applied
    caxis(gather(rangeAssign));
end

%Labels
xlabel('Inferior-Superior [mm]','FontSize',FontSize,'Interpreter','latex');
ylabel('Left-Right [mm]','FontSize',FontSize,'Interpreter','latex');
zlabel('Posterior-Anterior [mm]','FontSize',FontSize,'Interpreter','latex');

%Title
if ~isempty(titleName)
   sgtitle(titleName,'FontSize',1.5*FontSize,'Interpreter','latex');
end
