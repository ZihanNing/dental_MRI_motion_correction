
function visResidualsExt(We,outlWe,t,Residuals,ktraj,pau,folderName,fileName,numFig,titleName, mSt, stateSample)

%VISRESIDUALS   Visualizes reconstruction residuals
%   VISRESIDUALS({WE},{OUTLWE},{T},{RESIDUALS},{KTRAJ},PAU,{FOLDERNAME},{FILENAME})
%   * {WE} is the inverse of the median error per motion state normalized 
%   to the median among states
%   * {OUTLWE} are the detected outliers in the motion states
%   * {T} is the time of the motion states
%   * {RESIDUALS} are the residuals in the spectrum
%   * {KTRAJ} are the trajectories
%   * {PAU} indicates whether to pause the execution, it defaults to 1
%   * {FOLDERNAME} gives a folder where to write the results
%   * {FILENAME} gives a file where to write the results
%
%   Yannick Brackenier - 2023/10/20 - modified from visResiduals.m

if nargin<1;We=[];end
if nargin<2;outlWe=[];end
if nargin<3;t=[];end
if nargin<4;Residuals=[];end
if nargin<5;ktraj=[];end
if nargin<6 || isempty(pau);pau=1;end
if nargin<7;folderName=[];end
if nargin<8;fileName=[];end
if nargin<9 || isempty(numFig);numFig=400;end
if nargin<10 || isempty(titleName);titleName='';end
if nargin<11 || isempty(mSt);mSt=[];end
if nargin<12 || isempty(stateSample);stateSample=[];end

if ~iscell(folderName);folderName{1}=folderName;folderName{2}=folderName{1};end

co=[     0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330];

FontSizeA=30;
FontSizeB=24;
LineWidth=3;
MarkerSize=20;

if ~isempty(We)
    h = createFig(numFig);
    if ~isempty(t);xlab='$t$ (s)';else; xlab='$t$ (au)';end
    if isempty(t);t=1:NS;end
    plot(t,We(:),'Color',co(1,:),'LineWidth',LineWidth,'LineStyle','-')
    hold on
    plot(t(outlWe),We(outlWe),'*','Color',co(2,:),'MarkerSize',MarkerSize,'LineWidth',LineWidth);
    if t(end)>t(1);xlim([t(1) t(end)]);end
    xlabel(xlab,'Interpreter','latex','Color',[1 1 1],'FontSize',FontSizeA+6)
    grid on        
    set(gca,'Color','none','XColor',[1 1 1],'YColor',[1 1 1],'FontSize',FontSizeA)
    dire{1}='Weights';
    if any(outlWe(:)~=0);dire{2}='Outliers';end
    AX=legend(dire);
    LEG = findobj(AX);
    set(LEG,'Color','none','TextColor',[1 1 1],'Location','NorthWest','FontSize',FontSizeA)          
    set(gcf,'Color',[0 0 0])
    set(gcf, 'Position', get(0,'Screensize'))  
    if ~isempty(titleName);sgtitle(titleName,'Interpreter', 'latex', 'FontSize', FontSizeA,'Color',[1 1 1]);end
    if pau==2 && ~isempty(folderName{1}) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
        if ~exist(folderName{1},'dir');mkdir(folderName{1});end
        %print(strcat(folderName,filesep,fileName),'-dpng');
        export_fig(strcat(folderName{1},filesep,fileName,'.png'));
        close(h);
    end
end
if ~isempty(Residuals)
    for m=1:2;Residuals=fftshift(Residuals,m);end
    Nrep = max(multDimSize(Residuals,3:16));
    NN=size(Residuals);
    Residuals = Residuals(:,:);
    outlWeMat = Residuals(:,:);
    outlWeMatInv = zerosL(outlWeMat);
    ResidualTime = zeros([1 size(ktraj,1)],'like',Residuals);
    mStT = unique(mSt);
    for i=1:length(mStT)
       idx = find(stateSample==mStT(i));
       for j=1:length(idx)
           repId = ceil(idx(j)/length(stateSample)*Nrep);
           xx = ktraj(idx(j),1); 
           yy = ktraj(idx(j),2) + (repId-1)*NN(2);
           if outlWe(i);outlWeMat(xx,yy)=multDimMax(Residuals);end
           if outlWe(i);outlWeMatInv(xx,yy)=Residuals(xx,yy);end
           ResidualTime(idx(j)) = Residuals(xx,yy);
       end
    end
    h = createFig(numFig+1);
    subtightplot(2,3,1,.1,.04,.02)
    if ~isempty(ktraj)
        kMin=min(ktraj,[],1);
        kMax=max(ktraj,[],1);        
        for n=1:2;kGrid{n}=kMin(n):kMax(n);end
        kGrid{1}=kGrid{1}';
        kGrid{2}=repmat(kGrid{2},[1 Nrep]);    
        %for m=1:2;Residuals=fftshift(Residuals,m);end
        LResiduals=log(1+Residuals(:,:));
        imagesc('XData',kGrid{2},'YData',kGrid{1},'CData',mapToColormap(LResiduals,'jet'));
    else
        imagesc('CData',Residuals);
    end
    colormap(jet);%colorbar
    axis image
    title('log(1+$r$)','Interpreter','latex','FontSize',FontSizeB)
    xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeB)
    ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeB,'Rotation',0)
    
    subtightplot(2,3,2,.1,.04,.02)
    imagesc('CData',log(1+outlWeMat));
    colormap jet
    axis image
    title(sprintf('%s - Outliers',titleName),'Interpreter','latex','FontSize',FontSizeB)
    xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeB)
    ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeB,'Rotation',0)
    
    subtightplot(2,3,3,.1,.04,.02)
    imagesc('CData',log(1+outlWeMatInv));
    colormap jet
    axis image
    title(sprintf('%s - Outliers',titleName),'Interpreter','latex','FontSize',FontSizeB)
    xlabel('$k_3$','Interpreter','latex','FontSize',FontSizeB)
    ylabel('$k_2$','Interpreter','latex','FontSize',FontSizeB,'Rotation',0)
    
    subtightplot(2,3,4:6,.1,.04,.02)
    plot(log(ResidualTime))
    if ~isempty(outlWe) && ~isempty(stateSample)
       for i=find(gather(outlWe(:).'))
            x = find(stateSample==mStT(i));

            factBox=1;
            boxx=[x(1)-.5  x(1)-.5 x(end)+.5 x(end)+.5]; %+.5 because otherwise neighbourhing boxes have a small white space in between- only for visual purposes
            boxy=[0    1    1      0     ]*max(log(1+ResidualTime))*factBox;
            patch(boxx,boxy,[1 0 0],'FaceAlpha',0.2,'LineStyle','none');%'EdgeColor',[.2 0 0])% 'LineWidth',0.00001,   
       end
    end

    if pau==2 && ~isempty(folderName{2}) && ~isempty(fileName)%WE SIMPLY WRITE TO FILE
        if ~exist(folderName{2},'dir');mkdir(folderName{2});end
        %print(strcat(folderName,filesep,fileName),'-dpng');
        gcf;
        export_fig(strcat(folderName{2},filesep,fileName,'.png'));
        %close(h);
    end
end
if pau==1;pause;end


