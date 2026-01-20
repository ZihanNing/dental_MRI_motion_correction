
function [f, S, SFiltered, t, sFiltered] = freqAnal(s, deltaT, freqRem, windowRem, useFCT, numFig, titleName)

%FREQANAL analysese the signal in the frequency domain and allows to filter out hormonics.
%   [SFILTERED, T, SFILTERED, F, S] = freqAnal(S, DELTAT, FREQREM, WINDOWREM, USEFCT)
%   * S is the signal with dimensions [ NCha x NSamples ] 
%   * {DELTAT} is the sample duraction in seconds
%   * {FREQREM} is the interval of frequencies to remove
%   * {WINDOWREM} is the window to use when filtering in the frequency domain
%   * {USEFCT} is a flag whether to use the FCT instead of the FFT
%   ** SFILTERED the filtered signal in the temporal domain
%   ** T the time axis
%   ** SFILTERED the filtered signal in the frequency domain
%   ** F the frequency axis
%   ** S the original signal in the frequency domain
%
% Yannick Brackenier 

%p has dimensions [ NCha x NSamples ] 
%deltaT in seconds

if nargin<2 || isempty(deltaT);deltaT=1;warning('freqAnal:: No deltaT provided, so frequency axis will not be appropriate.');end
if nargin<3 || isempty(freqRem);freqRem=[];end
if nargin<4 || isempty(windowRem);windowRem=[];end
if nargin<5 || isempty(useFCT);useFCT=1;end
if nargin<6 || isempty(numFig);numFig=1;end
if nargin<7 || isempty(titleName);titleName='';end

%%% HANDLE DATA SIZES
NSamples= size(s,2);
NCha= size(s,1);

Fs = 1/deltaT;%Sampling frequency
t = (0:NSamples-1)*deltaT;        % Time vector

dimTime = 2;
shiftFFT = 1;

%%% TRANFORM TO FREQUENCY DOMAIN
if useFCT
    S = fct(s,NSamples,dimTime);
    f = ((0:1/NSamples:1-1/NSamples)*Fs)/2;
else
    S = fft(s,NSamples,dimTime);
    f = ((0:1/NSamples:1-1/NSamples)*Fs);
    if shiftFFT
       S = fftshift(S);
       f = f- f(ceil((NSamples+1)/2) );
    end
end

%%% FILTER (IN PROGRESS!!)
if nargout>2
    sFiltered = s;
    
    if ~isempty(freqRem)
        idxZero = abs(f)>=freqRem(1) & abs(f)<=freqRem(2);%Make sure negative frequencies are also removed
        SFiltered = S;
        SFiltered(:,idxZero) = 0;%Need to add windowing here

        if useFCT
            sFiltered = ifct(SFiltered,NSamples,dimTime);
        else
            if shiftFFT
                sFiltered = ifft(ifftshift(SFiltered),NSamples,dimTime);
            else
                sFiltered = ifft(SFiltered,NSamples,dimTime);
            end
        end
    end   
end

%%% PLOT
FontSize=10;%Baseline fontsize
ModLab=7;%Axes labels
ModTit=7;%Subtitles
ModSupTit=10;%Suptitle
ModTick=3;%Ticks of the axes
LineWidth=1;
invCol=1;%1 generates white background

%%% PLOT
NPlot = [2 1];
xLabel = 'Frequency [Hz]';
if ~isempty(numFig);h=figure(numFig);clf;else; h = figure();end
set(h,'Color',[0 0 0]+invCol,'Position', get(0,'Screensize'));

%Magnitude
subplot(NPlot(1),NPlot(2),1); plot(f,abs(S),'LineWidth',LineWidth,'LineStyle','-')
set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
ylabel('$|S_{PT}|$ [a.u.]','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
axis tight;
title('\textbf{Magnitude}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)
%Phase
subplot(NPlot(1),NPlot(2),2); plot(f, angle(S),'LineWidth',LineWidth,'LineStyle','-')
set(gca,'Color','none','XColor',[1 1 1]*(1-invCol),'YColor',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTick,'TickLabelInterpreter','latex')
xlabel(xLabel,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
ylabel('\angle$S_{PT}$ [rad]','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModLab)
axis tight;
title('\textbf{Phase}','Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModTit)

%%% SUPTITLE
if ~strcmp(titleName,''); sgtitle(titleName,'Interpreter','latex','Color',[1 1 1]*(1-invCol),'FontSize',FontSize+ModSupTit);end
   
