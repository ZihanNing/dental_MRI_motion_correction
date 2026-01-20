
function [rec] = preProcessPTInRec(rec, type, timeWindow, blockLength, deb)

if nargin<4 || isempty(blockLength);blockLength=length(rec.Assign.z{2});end
if nargin<5 || isempty(deb);deb=0;end

%timeWindow in ms

%%% PARAMETERS
pTime = rec.PT.pSliceImage;
N = size(pTime);N(end+1:16)=1;
[~, idx] = PESamplesFromRecPTTest(rec);
perm = [4 2 1 3 5:16];

filtKernelWidthidx = round(timeWindow./rec.Par.Labels.RepetitionTime);
filtKernelWidthidx = filtKernelWidthidx-mod(filtKernelWidthidx,2)-1;%make odd
plotType = 1;
coilPlot = 1:5;

%%% RESPAPE INTO TIME PROFILES
for n=2:3
    pTime=fftshiftOperator(pTime,2,1,n);
    pTime=fftGPU(pTime,n)/N(n);
    pTime=fftshiftGPU(pTime,n);
end

pTime = resPop(pTime, 2:3,[],2);
pTimeFiltFinal = zerosL(pTime);

%%% ORDER IN TIME
pTimeTemp=[];
for i=1:N(5)%Multiple repeats
    idxToExtract = (i-1)*(length(idx)/N(5))+1:(i)*(length(idx)/N(5)) ;
    pTimeTemp=cat(2,  pTimeTemp, ...
                  dynInd( dynInd(pTime,i,5),...
                          dynInd(idx,idxToExtract,2) ,2));
end
pTime = pTimeTemp; pTimeTemp = [];

pTime = permute(pTime,perm);%sorted
phaseRem = sign(pTime(1,:));
pTime = pTime./phaseRem;
pTime = pTime./sqrt(normm(pTime,[],1));

visPTSignal(dynInd(pTime,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 200, replace(sprintf('\\textbf{Original PT signal for first %d coils:}\n %s',length(coilPlot),''),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    

%pTime = pTime./sqrt(normm(pTime,[],1));
%visPTSignal(dynInd(pTime,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 201, replace(sprintf('\\textbf{Normalised  PT signal for first %d coils:}\n %s',length(coilPlot),''),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    

%%% FILTER
TFEfactor = rec.Par.Labels.TFEfactor;
if TFEfactor<length(rec.Assign.z{2});blockLength = TFEfactor;end
filtKernelWidthidx = min(filtKernelWidthidx, blockLength);
filtKernelWidthidx = filtKernelWidthidx-mod(filtKernelWidthidx,2)-1;%make odd

shotId = 0;

for i=1:blockLength:size(pTime,2)
    shotId=shotId+1;
    idxShot = i:min(i+blockLength-1,size(pTime,2));
        
    if deb
        createFig(100);
        subplot 211; plot(abs(dynInd(pTime,idxShot,2)).','b') ; hold on
        subplot 212; plot(angle(dynInd(pTime,idxShot,2)).','b') ; hold on
    end
    
    %Median filtering
    if strcmp(type,'med') && filtKernelWidthidx>1
        pTime = dynInd(pTime,idxShot,2,   cdfFilt(real(dynInd(pTime,idxShot,2)),'med',[1 filtKernelWidthidx 1 1 ],'replicate') + ...%real channel
                                   1i*cdfFilt(imag(dynInd(pTime,idxShot,2)),'med',[1 filtKernelWidthidx 1 1 ],'replicate')...%imaginary channel
                                   );
    end   
    %Savitzky-Golay filtering
    if strcmp(type,'golay') &&  filtKernelWidthidx>1
        order = 2;
        pTime = dynInd(pTime,idxShot,2, single(sgolayfilt(dynInd(double(pTime),idxShot,2),order,double(filtKernelWidthidx),[],2)) );
    end    
    
    if deb
        subplot 211; plot(abs(dynInd(pTime,idxShot,2)).','r') 
        subplot 212; plot(angle(dynInd(pTime,idxShot,2)).','r') 
        legend({'Original','Filtered'})
        sgtitle(sprintf('Shot %d',shotId))
        if deb>1;pause();else;pause(.1);end
    end
end

pTimeFilt = pTime;
visPTSignal(dynInd(pTimeFilt,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 202, replace(sprintf('\\textbf{Filtered PT signal for first %d coils:}\n %s',length(coilPlot),''),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    

%pTimeFiltTest = pTimeFilt./sqrt(normm(pTimeFilt,[],1));
%visPTSignal(dynInd(pTimeFiltTest,{coilPlot,rec.PT.idxMB},[1 3]), [], [],[], plotType, [], [], 203, replace(sprintf('\\textbf{Normalised filtered PT signal for first %d coils:}\n %s',length(coilPlot),''),'_',' '), [], [],'Time [s]',[],rec.Par.Labels.RepetitionTime/1000);    

pTimeFilt = pTimeFilt.*phaseRem;

rec.PT.pTimeTest = pTimeFilt;

%%% Put back
pTimeFilt = ipermute(pTimeFilt,perm);%sorted

for i=1:N(5)%Multiple repeats
    idxToExtract = (i-1)*(length(idx)/N(5))+1:(i)*(length(idx)/N(5)) ;
    pTimeFiltFinal = dynInd(pTimeFiltFinal, {idx(idxToExtract), i}, [2 5],  dynInd(pTimeFilt,idxToExtract,2) );
end

pTimeFilt = pTimeFiltFinal;pTimeFiltFinal=[];
pTimeFilt = resSub(pTimeFilt, 2:3,N(2:3));

for n=2:3
    pTimeFilt=ifftshiftGPU(pTimeFilt,n);%fftshift since in solveXT there is an iffthift on timeIndex
    pTimeFilt=ifftGPU(pTimeFilt,n)*N(n);
    pTimeFilt=fftshiftOperator(pTimeFilt,2,0,n);
end

rec.PT.pSliceImage = pTimeFilt;
