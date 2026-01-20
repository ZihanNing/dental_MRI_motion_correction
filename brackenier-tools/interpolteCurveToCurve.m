

function [cInt, c1New, cInsert, c2New] = interpolteCurveToCurve(c1, c2, NInsert, NFix1 , NFix2, offsetFix1, offsetFix2, typeInterp,deb)

%INTERPOLATECURCETOCURVE interpolates 2 ND-curves by independently interpolating each dimension.
%   [CINT,C1NEW,CINSERT,C2NEW] = INTERPOLATECURCETOCURVE(C1,C2,{NINSERT},{NFIX1},{NFIX2},{OFFSETFIX1},{OFFSETFIX2},{TYPEINTERP},{DEB})
%   * C1 is the first curve as a [NSamples x NDimensions] array
%   * C1 is the second curve
%   * {NINSERT} is the number of samples to insert between curves
%   * {NFIX1} is the number of samples from curve 1 to use (and fix) for interpolation 
%   * {NFIX2} is the number of samples from curve 2 to use (and fix) for interpolation 
%   * {OFFSETFIX1} is the offset w.r.t. the edge of curve 1 to set the samples to use for interpolation
%   * {OFFSETFIX2} is the offset w.r.t. the edge of curve 2 to set the samples to use for interpolation
%   * {TYPEINTERP} is the type of interpolation to use. Defaults to 'spline'
%   * {DEB} is a debug flag to plot the interpolation
%   ** CINT is the combined interpolated curve
%   ** C1NEW is the part that originally belonged to curve 1.
%   ** CINSERT is the inserted interpolated curve.
%   ** C2NEW is the part that originally belonged to curve 2.
%
%   Yannick Brackenier 2022-08-22

if nargin<3 || isempty(NInsert); NInsert=round((size(c1,1)+size(c1,1))/4);warning('interpolteCurveToCurve:: No number of interpolation points set. Default value used.');end
if nargin<4 || isempty(NFix1); NFix1=round(size(c1,1)/10);end
if nargin<5 || isempty(NFix2); NFix2=round(size(c2,1)/10);end
if nargin<6 || isempty(offsetFix1); offsetFix1=0;end
if nargin<7 || isempty(offsetFix2); offsetFix2=0;end
if nargin<8 || isempty(typeInterp); typeInterp='spline';end
if nargin<9 || isempty(deb); deb=0;end

%%% SET PARAMETERS
N1 = size(c1,1);
N2 = size(c2,1);
NTot = N1+NInsert+N2;
NDim = size(c1,2);
assert(size(c1,2)==size(c2,2),'interpolteCurveToCurve:: Dimensions of curves not consistent.');

%%% CHECK FOR INCORRECT PARAMETERS
assert(mod(NInsert,1)==0 && NInsert>=0,'interpolteCurveToCurve:: NInsert shoud be positive integer.');
assert(mod(NFix1,1)==0 && NFix1>0 && mod(NFix2,1)==0 && NFix2>0,'interpolteCurveToCurve:: NFix shoud be positive integer.');
assert(mod(offsetFix1,1)==0 && offsetFix1>=0 && mod(offsetFix2,1)==0 && offsetFix2>=0,'interpolteCurveToCurve:: offsetFix shoud be positive integer.');
if N1==0; assert(NFix1==1 && offsetFix1==0,'interpolteCurveToCurve:: If frist curve is a single start point, NFix1 should be 1 (was %d) and offsetFix1 should be 0 (was %d).',NFix1,offsetFix1);end
if N2==0; assert(NFix2==1 && offsetFix2==0,'interpolteCurveToCurve:: If second curve is a single end point, NFix2 should be 1 (was %d) and offsetFix2 should be 0 (was %d).',NFix2,offsetFix2);end

%Make sure offset is not larger than the array size
if offsetFix1>N1;offsetFix1=N1;warning('interpolteCurveToCurve:: offsetFix1 was set too large and was adapted (to %d).',offsetFix1);end
if offsetFix2>N2;offsetFix2=N2;warning('interpolteCurveToCurve:: offsetFix1 was set too large and was adapted (to %d).',offsetFix2);end

%Adjust NFix if too large
if (N1-NFix1-offsetFix1+1)<1; NFix1=N1-offsetFix1-1+1;warning('interpolteCurveToCurve:: NFix1 was set too large and was adapted (to %d) to start at first element of c1.',NFix1);end
if (N1+NInsert+offsetFix2+NFix2)>NTot; NFix2=NTot - (N1+NInsert+offsetFix2 );warning('interpolteCurveToCurve:: NFix2 was set too large and was adapted (to %d) to end at the last element of c2.',NFix2);end


%%% MAKE PARAMETERS COMPATIBLE
cInsert = zeros([NInsert, NDim]);
cInt = cat(1,c1,cInsert,c2);

%%% SET INDICES
%Indices of the different parts
idxInt = 1:NTot;
idx1 = 1:N1;
idx2 = (N1+1+NInsert):NTot;
idxInsert = (N1+1):N1+NInsert;
assert(isequal(idxInt, cat(2,idx1,idxInsert,idx2)),'interpolteCurveToCurve:: Indices do not add up.');

%Indecis of the parts to use as lookup values during interpolation
idxFix1 = [(N1-NFix1+1):N1 ] - offsetFix1;
idxFix2 = (N1+1+NInsert+offsetFix2):(N1+NInsert+offsetFix2+NFix2);
idFix = cat(2, idxFix1,idxFix2);

%Indecis of the quary points during interpolation
idxQuery = idxFix1(end)+1:(idxFix2(1)-1);
NQuery=length(idxQuery);
assert(numel(idxQuery)==NQuery,'interpolteCurveToCurve:: Indices do not add up.');

%%% PERFORM INTERPOLATION
valQuery = zeros([NQuery NDim],'like',c1);
for i=1:NDim
    %Assign values
    valFix1 =  cInt(idxFix1,i);
    valFix2 =  cInt(idxFix2,i);
    valFix = cat(1,valFix1,valFix2);

    %Interpolate
    valQuery(:,i)=interp1(idFix,valFix,idxQuery,typeInterp);
end
cInt(idxQuery,:)=valQuery;

%%% ASSIGN
c1New = cInt(idx1,:);
c2New = cInt(idx2,:);
cInsert = cInt(idxInsert,:);

%%% REPORT
if deb
    cOrig = cat(1,c1,zerosL(cInsert),c2);
    %%% Set parameters
    FontSize=10;%Baseline fontsize
    ModLab=7;%Axes labels
    ModTit=5;%Subtitles
    ModSupTit=13;%Suptitle
    ModTick=5;%Ticks of the axes
    ModLeg=5;%Legends
    LineWidth=1.5;%Linewidth of the plots
    invCol=1;%1 generates white background
    useScatter=1;
    
    %%% Create figure
    h=figure();clf;set(h,'color',[0 0 0]+invCol);%,'Position',get(0,'ScreenSize'))
    for i=1:NDim
        subplot(NDim,1,i);hold on;
        p2=plot_helper(idx1,cInt(idx1,i),LineWidth,'r',useScatter);
        p1=plot_helper(idx1,cOrig(idx1,i),LineWidth,'b',useScatter);
        plot_helper(idx2,cInt(idx2,i),LineWidth,'r',useScatter);
        plot_helper(idx2,cOrig(idx2,i),LineWidth,'b',useScatter);
        p3=plot_helper(idxInsert,cInt(idxInsert,i),LineWidth,'g',useScatter);
        axis([-inf inf -inf inf])
        if i==1;hL=legend([p1, p2,p3],{'$c_{orig}$', '$c_{interp}$', '$c_{insert}$'}, 'Location', 'NorthEast','Interpreter','Latex','FontSize',FontSize+ModLeg);end
        title(sprintf('$\\textbf{Dimension %d}$',i),'Interpreter','Latex','FontSize',FontSize+ModTit)
        
        xlabel(sprintf('$x_%d$ [a.u.]',i),'Interpreter','Latex','FontSize',FontSize+ModLab)
        ylabel(sprintf('$c_%d$ [a.u.]',i),'Interpreter','Latex','FontSize',FontSize+ModLab)
        set(gca,'TickLabelInterpreter','Latex','FontSize',FontSize+ModTick);
    end
    sgtitle('\textbf{Interpolation performance}','Interpreter','Latex','FontSize',FontSize+ModSupTit)
    
end
end

%%% HELPER FUNCTION
function [p] = plot_helper(idx,c,LineWidth,col,useScatter)
    if useScatter
      p = scatter(idx,c,[],col,'filled');
      p.SizeData = 10;
    else
      p=plot(idx,c,'Linewidth',LineWidth,'color',col);  
    end
   
end



