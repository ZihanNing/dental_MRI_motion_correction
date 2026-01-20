
function [PT, parPT] = preProcessPT(PT, parPT, NY, kIndex, deb, parPreProcessing, TFEfactor)

if nargin<3 || isempty(NY);NY=[];end
if nargin<4 || isempty(kIndex);kIndex=[];end
if nargin<5 || isempty(deb);deb=1;end
if nargin<6 || isempty(parPreProcessing);parPreProcessing=[];end
if nargin<7 || isempty(TFEfactor);TFEfactor=[];end

%%% SET ASIDE VARIABLES
pSlice = PT.pSlice;
reOrderTemporally = any(parPT.signalUsage.eigTh<1) || parPT.PTFilter.medFiltKernelWidthidx>1 || parPT.PTFilter.golayFiltKernelWidthidx>1; 
reOrderTemporally = reOrderTemporally || ~isempty(parPT.signalUsage.orderPreProcessing) ;%only for debugging during development
reOrderTemporally = reOrderTemporally || parPT.Calibration.calibrationOrder>1;
fprintf('<strong>Pilot Tone</strong>: Pre-processing\n')

%%% LOG THE PRE-PROCESSING STEPS
PT.preProcessing = [];
forceProcessing = parPT.Calibration.externalFit.Flag || ~isempty(parPT.signalUsage.preProcessing) ;
if forceProcessing && parPT.Calibration.externalFit.Flag
    parPTpreProcessing = parPT.Calibration.externalFit.preProcessing;
    fprintf('	Enforcing pre-processing steps from the pre-calibrated acquisition: %s.\n',parPTpreProcessing.preProcessName);
    if ~isempty(parPT.signalUsage.preProcessing);warning('preProcessPT:: Taking pre-processing parameters from the pre-calibration scan, although pre-processing parameters from a separate acquisition were also provided in parPT.signalUsage.preProcessing.\n');end
elseif forceProcessing && ~isempty(parPT.signalUsage.preProcessing)
    parPTpreProcessing = parPT.signalUsage.preProcessing;
    fprintf('	Enforcing pre-processing steps from separate acquisition: %s.\n',parPTpreProcessing.preProcessName);
else
    fprintf('	Not enforcing pre-processing steps from pre-calibrated acquisition.\n'); 
end
%if forceProcessing; assert(isequal(parPT.signalUsage.orderPreProcessing,[]),'Pre-processing steps should be the same for pre-calibration use.');end

%%% REMOVE COILS FROM INFERIOR-SUPERIOR DIRECTION
PT.NChaPTAcq=size(pSlice,4);% Number of non-compressed channels
if ~isempty(parPT.signalUsage.SupportIS) && ~isempty(parPreProcessing) && isfield(parPreProcessing,'coilGeom')
    assert(length(parPreProcessing.coilGeom.idxIS )==PT.NChaPTAcq,'preProcessPT:: Coil geometry information not consistent with provided PT signal.');
    coilIdx = max(1,round(parPT.signalUsage.SupportIS(1)*PT.NChaPTAcq)):round(parPT.signalUsage.SupportIS(2)*PT.NChaPTAcq);
    coilIdx = parPreProcessing.coilGeom.idxIS (coilIdx);
    if forceProcessing
        pSlice = dynInd(pSlice,parPTpreProcessing.coilIdx,4);
    else
        PT.preProcessing.coilIdx = coilIdx;%Track
        pSlice = dynInd(pSlice,coilIdx,4);
    end  
end

%%% BROADBAND HANDLING
if ~isempty(PT.idxMB)
    if parPT.signalUsage.combineMB==1%Extract peak in RO direction
        pSlice = dynInd(pSlice, PT.idxMB,3);
    elseif parPT.signalUsage.combineMB==2%Average frequencies
        pSlice = multDimMea(pSlice, 3);
    else
       assert(size(pSlice,3)==1, 'Reconstruction not able to work with multi-band Pilot Tone signal.');
    end
    fprintf('	Removing multi-band PT signal.\n')
end

%%% RE-ORDER IN TEMPORAL DOMAIN IF NEEDED
if reOrderTemporally
    %%% Move the PT signal to the hybrid space where samples have meaning
    for n=1:2
        %pSlice=fftshiftGPU(pSlice,n);%Already done in solveXTB_PT script.
        pSlice=fftGPU(pSlice,n)/(NY(n));
        pSlice=fftshiftGPU(pSlice,n);
    end
    %%% Order in time
    idx = sub2ind(NY(1:2),kIndex(:,1),kIndex(:,2))';%%% Get the index for the sampling
    pTime=[];
    for i=1:NY(5)%Multiple repeats
        idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
        pTime=cat(1,  pTime, ...
                      dynInd( resPop(dynInd(pSlice,i,5) ,1:2,[],1),...
                              dynInd(idx,idxToExtract,2) ,1));
    end
else
    pTime = pSlice;
end

%%% PHASE REFERENCING
if parPT.signalUsage.referencePhase>0
   pTime = pTime./sign(dynInd(pTime,parPT.signalUsage.referencePhase,4));
   fprintf('	Referencing phase of PT signal w.r.t coil %d.\n',parPT.signalUsage.referencePhase)
end

if deb && reOrderTemporally
    tt=permute(pTime,[1 4 2:3]).';
    visPTSignal(tt,'pTime',[],[],2,[],[],89,'Original Pilot Tone signal');
end

%%% TIME-INDEPENDENT PRE-PROCESSING
for step = parPT.signalUsage.orderPreProcessing
    if step == 1 %Removing complex mean values along PE directions
        if forceProcessing
            PT.preProcessing.step1_mean = parPTpreProcessing.step1_mean;
            pTime = bsxfun(@minus, pTime, parPTpreProcessing.step1_mean);
        else
            PT.preProcessing.step1_mean = multDimMea((pTime),1:2); %Track
            pTime = bsxfun(@minus, pTime, multDimMea((pTime),1:2));
        end
        fprintf('	Removing complex mean values along PE directions.\n')
    elseif step == 2 %Removing mean phase along PE directions
        if forceProcessing
            PT.preProcessing.step2_meanPhase = parPTpreProcessing.step2_meanPhase;
            pTime = bsxfun(@rdivide, pTime, parPTpreProcessing.step2_meanPhase);
        else
            PT.preProcessing.step2_meanPhase = sign(multDimMea(pTime,1:2)); %Track
            pTime = bsxfun(@rdivide, pTime, sign(multDimMea(pTime,1:2)));
        end   
        fprintf('	Removing mean phase along PE directions.\n')
    elseif step == 3 %Normalising along coil dimension
        %PT.preProcessing.step3_norm = sqrt(normm(pTime,[],4)); %Track
        pTime = bsxfun(@rdivide, pTime , sqrt(normm(pTime,[],4)) );
        fprintf('	Normalising along coil dimension.\n')
    elseif step == 4 %Whitening complex signal
        %permWhit = [1 4 2:3];
        %percSampes=[];
        %pTime = ipermute(whiten(permute(pTime,permWhit),percSampes),permWhit);
        %fprintf('	Whitening complex signal.\n')
        error('Do not use this.')
    elseif step == 5 %Making unit-variance signal
        %pTime = bsxfun(@rdivide, pTime, std((pTime),[],4)+0*std(imag(pTime),[],4));
        %fprintf('	Making unit-variance signal.\n')
        error('Do not use this.')
    elseif step == 6 %Concatenate real and imag component
        pTime = cat(4,real(pTime),imag(pTime));
        fprintf('	Concatenate real and imag component.\n')
    elseif step == 7 %Magnitude
        pTime = abs(pTime);
        fprintf('	Taking magnitude.\n')
    end
end

if deb && reOrderTemporally
    tt=permute(pTime,[1 4 2:3]).';
    visPTSignal(tt,'pTime',[],[],2,[],[],90,'After time-independent pre-processing');
    %[~,pReshaped] = reshapeP(rand([6 32]),permute(pTime,[1 4 2:3]).','mat2vec');
end
                
%%% DIMENSIONALITY REDUCTION USING SINGULAR VALUE DECOMPOSITION
if any(parPT.signalUsage.eigTh~=1)  
    %%% Perform SVD
    permSVD = [1 4 2:3];
    pTime = permute(pTime, permSVD);
    if parPT.signalUsage.useRealImagSVD && ~isreal(pTime); pTime = cat(2,real(pTime),imag(pTime)); end
    pTimeOrig = pTime;
    [U,S,V] = svd(pTime,'econ');
    
    %%% Select #singular vectors to keep
    sv = diag(S);
    nv = zerosL(parPT.signalUsage.eigTh);
    for l=1:length(nv)
        if parPT.signalUsage.eigTh(l)==0 %Automatic/tuned detection
            [~,nv(l)]=screePoint(sv,.7);
        elseif parPT.signalUsage.eigTh(l) < 0 && mod(parPT.signalUsage.eigTh(l),1)==0 %-eigTh is the pre-defined number of components
            nv(l) = -parPT.signalUsage.eigTh(l);
        else %eigTh is fraction of sv to discard
            nv(l)= gather( find(sv>=sv(1)*parPT.signalUsage.eigTh(l),1,'last'));
        end
    end 
    if forceProcessing
       PT.preProcessing.S = parPTpreProcessing.S;
       PT.preProcessing.V = parPTpreProcessing.V;
       pTime = pTime*parPTpreProcessing.V(:,1:max(nv));
       if ~parPT.signalUsage.includeSVScaling;pTime=pTime*inv(parPTpreProcessing.S(1:max(nv),1:max(nv))); end
    else
       PT.preProcessing.S = S; %Track
       PT.preProcessing.V = V; %Track
       pTimeTest = pTime;
       pTime = U(:,1:max(nv));
       if parPT.signalUsage.includeSVScaling;pTime=pTime*S(1:max(nv),1:max(nv)); end
    end
    pTime = ipermute(pTime, permSVD);
    
    
    %%% Add a reporting stage with the number of componenents in each level
    if 0%deb
        pTimeProj = U(:,1:max(nv))*S(1:max(nv),1:max(nv))*V(:,1:max(nv))';
        pTimeProjError = pTimeOrig-pTimeProj;
        visSVDThresholding(U,S,V,max(nv),91,[],0);
        visPTSignal(pTimeOrig.','pTime',[],[],[],[],[],93,'Original signal');
        visPTSignal(pTimeProj.','pTime',[],[],[],[],[],94,'Projected signal');        
        visPTSignal(pTimeProjError.','pTime',[],[],[],[],[],95,'Error');        
    end
    fprintf('	Low rank respresentation: using %s components.\n',sprintf('%d-',nv));
else
   nv = size(pTime,4)*onesL(parPT.signalUsage.eigTh);
end

%%% INFER DATA TYPE AND RE-ARRANGE
if isreal(pTime)
    nvSave = ceil(max(nv)/2); %Needs to be even when we store it as complex numbers again
    if size(pTime,4)<2*nvSave; pTime = cat(4, pTime , dynInd(zerosL(pTime),1:(2*nvSave-size(pTime,4)),4));end
    pTime = dynInd(pTime,1:nvSave,4) + 1i*dynInd(pTime,nvSave+1:2*nvSave,4);%Store as complex array for fft resampling functionality
    PT.useRealImag = 1;
    PT.NChaPTRec = nv;
else
    PT.useRealImag=0; 
    PT.NChaPTRec = nv;
end
    
%%% FILTER IN TIME DOMAIN
medFiltKernelWidth = round(parPT.PTFilter.medFiltKernelWidthidx);if medFiltKernelWidth>1; fprintf('	Filtering with %d-element median kernel.\n',medFiltKernelWidth);end
golayFiltKernelWidth = round(parPT.PTFilter.golayFiltKernelWidthidx);if golayFiltKernelWidth>1; fprintf('	Filtering with %d-element Savitzky-Golay kernel.\n',golayFiltKernelWidth);end

if isempty(TFEfactor); TFEfactor = size(pTime,1);end
assert(mod( size(pTime,1)/TFEfactor,1)==0,'preProcessPT: Shot ordering for TFE sequence not correct.')

for shot = 1:(size(pTime,1)/TFEfactor)
    idxShot = (shot-1)*TFEfactor + [1:TFEfactor];
    %Median filtering
    if medFiltKernelWidth>1
        pTime = dynInd(pTime,idxShot,1,   cdfFilt(real(dynInd(pTime,idxShot,1)),'med',[medFiltKernelWidth 1 1 1 ],'replicate') + ...%real channel
                                   1i*cdfFilt(imag(dynInd(pTime,idxShot,1)),'med',[medFiltKernelWidth 1 1 1 ],'replicate')...%imaginary channel
                                   );
    end   
    %Savitzky-Golay filtering
    if golayFiltKernelWidth>1
        order = 2;
        pTime = dynInd(pTime,idxShot,1, single(sgolayfilt(dynInd(double(real(pTime)),idxShot,1),order,golayFiltKernelWidth,[],1)) + ...real channel
                                    1i* single(sgolayfilt(dynInd(double(imag(pTime)),idxShot,1),order,golayFiltKernelWidth,[],1))) ;
    end    
end
if deb;visPTSignal(squeeze(pTime).','pTime',[],[],2,[],[],96,'After filtering');end

%%% UPDATE CALIBRATION ORDER
if parPT.Calibration.calibrationOrder>1
    nvOrig = nv;
    %Undo complex storage
    if PT.useRealImag; pTime = cat(4, real(pTime) , imag(pTime)); end
    %Concatenate higher-order terms
    pTimeOrig = dynInd(pTime,1:max(nv),4);
    pTime = pTimeOrig;
    for order = 2:parPT.Calibration.calibrationOrder
        pTime = cat(4, pTime, pTimeOrig.^order);
    end
    pTimeOrig=[];
    %Store as complex numbers
    nv = size(pTime,4)*onesL(nv);
    if PT.useRealImag
        nvSave = ceil(max(nv)/2); %Needs to be even when we store it as complex numbers again
        if size(pTime,4)<2*nvSave; pTime = cat(4, pTime , dynInd(zerosL(pTime),1:(2*nvSave-size(pTime,4)),4));end
        pTime = dynInd(pTime,1:nvSave,4) + 1i*dynInd(pTime,nvSave+1:2*nvSave,4);%Store as complex array for fft resampling functionality
    end
    PT.NChaPTRec = nv;
	fprintf('	Increasing order of the PT signal to %d. From %d --> %d (virtual) channels.\n',parPT.Calibration.calibrationOrder, max(nvOrig),max(nv));
end

%%% CONVERT TO SLICE DOMAIN
if reOrderTemporally
    %%% Create k-space array
    pSlice = zerosL(pSlice);
    pSlice = resPop(pSlice,1:2,[],1);%Flatten PE dimensions
    pSlice = dynInd(pSlice, 1:size(pTime,4),4);%Make coil dimensions compatible with pre-processing
    
    %%% Deal with multiple averages
    for i=1:NY(5)%Multiple repeats
        idxToExtract = (i-1)*(size(kIndex,1)/NY(5))+1:(i)*(size(kIndex,1)/NY(5)) ;
        pSlice = dynInd(pSlice, {idx(idxToExtract), i} , [1 5] , dynInd(pTime, idxToExtract,1)  );
    end

    pSlice = resPop(pSlice,1,NY(1:2),1:2);%Reshape to PE dimensions
        
    %%% Move the PT signal back to the image domain
    for n=1:2
        pSlice=ifftshiftGPU(pSlice,n);
        pSlice=ifftGPU(pSlice,n)*(NY(n));
        %pSlice=fftshiftGPU(pSlice,n);
    end
else
    pSlice = pTime;
end

%%% ASSIGN NEW VARIABLES
PT.pSlice = pSlice;

