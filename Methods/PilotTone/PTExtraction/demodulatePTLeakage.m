
function [y, pPeakPad, idxPeakPad, idxPeakPadAverage] = demodulatePTLeakage(y, dim, padFac, windowCent, windowWidth, gibbs, facFOVTh, isImageDomain,averageFreqCoils,Lin,Par)

if nargin<2 || isempty(dim);dim=numDims(y);end 
if nargin<3 || isempty(padFac);padFac=20;end
if nargin<4 || isempty(windowCent);windowCent=[];end
if nargin<5 || isempty(windowWidth);windowWidth=[];end
if nargin<6 || isempty(gibbs);gibbs=0;end%No windowing as interpolation using Fourier domain (~= zero-padding)
if nargin<7 || isempty(facFOVTh);facFOVTh=.45;end
if nargin<8 || isempty(isImageDomain);isImageDomain=0;end%If the readout is in the image domain
if nargin<9 || isempty(averageFreqCoils);averageFreqCoils=0;end%Per TR, enforce that every coil has the same index
if nargin<10 || isempty(Lin);Lin=[];end
if nargin<11 || isempty(Par);Par=[];end

forward = @(x) fftshift(ifft(ifftshift(x)));
backward= @(x) fftshift(fft(ifftshift(x)));
deb = 1;

%%% SHUTTER HANDLING
M = zeros(multDimSize(y,2:3));
if isempty(Lin) || isempty(Par)
    M = onesL(M);
else
    for i=1:length(Lin); M(Lin(i),Par(i)) = 1;end
end
M = resSub(M,1:2);
idShutter = M==0;

%%% PERMUTE AND FLATTEN ARRAY
nD = max(5,numDims(y));
perm = 1:nD; perm(perm==dim)=[]; perm=[dim perm];
perm = [perm(1:3) nD perm(4:end-1)];
y = permute(y,perm);
NY = size(y);
y = resSub(y,2:3);
y = dynInd(y,~idShutter,2);

NY2 = size(y);
y = resSub(y,2:3);

%y = resSub(y,2:nD-1);
dimNew = 1;

%%% ZERO-PADDING DIMENSIONS
N = size(y,1);
NPad = padFac*N;

%%% MOVE TO FOURIER DOMAIN (specific to TWIX2rec.m convention)
if isImageDomain
    y = fftshiftOperator(y, 2, 1, dimNew);
    y = ifft(y, [], dimNew);
    y = fftshiftOperator(y, 1, 1, dimNew);
end

%%% CREATE NEW ARRAYS
NP = size(y); NP(1)=1; %Only one element extracted
pPeakPad = zeros(NP,'like',y);
idxPeakPad = zerosL(pPeakPad);

%%% EXTRACT PART OF THE ARRAY IN THE DIMENSION OF INTEREST (ECC etc.)
croppOffset = [0 0];
if ~isempty(windowWidth) && ~isempty(windowCent)
    NH = onesL(size(y));
    NH(dimNew)=NY(dimNew);
    h=PTROWindow(NH, windowWidth, windowCent, gibbs);
    temp = find(h);
    croppOffset = [temp(1)-1 N-temp(end)];
    croppOffset = max(croppOffset)*onesL(croppOffset);%works better when same pre- and post-cropping
    if deb;fprintf('    Cropping offset: %d-%d\n',croppOffset);end
end

%%% LOOP OVER BLOCKS OF ARRAY IN OTHER DIMENSIONS
pPadAccum = 0;
blSz = 500;
reverseStr=[];
for i=1:blSz:size(y,2)
        
    %Block index
    vI=i:min(i+blSz-1,size(y,2));
    if deb
        msg = sprintf('    Performing sinc-interpolation on samples: %d/%d\n',vI(end),size(y,2));
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    yBlock = y(:,vI,:);
    
    %Crop
    yBlockPad = dynInd(yBlock, (croppOffset(1)+1):(N-croppOffset(2)),1);
    NCrop = size(yBlockPad,1);
    NNew = NCrop+NPad;

    %Zero-pad in Fourier domain
    yBlockPad = resampling(yBlockPad, [NNew, length(vI)], 2);%padArrayND(yBlock,padDim,1,0,'both');

    %%% Move to image domain (specific to TWIX2rec.m convention)
    yBlockPad = fftshiftOperator(yBlockPad, 1, 0, dimNew);
    yBlockPad = fft(yBlockPad, [], dimNew);
    yBlockPad = fftshiftOperator(yBlockPad, 2, 0, dimNew);

    %%% Find peak
    considerAll = 0;
    if ~considerAll
        yBlockPadZero = fillCentPT(yBlockPad, dimNew, 0, facFOVTh);
        orig=round(size(yBlockPadZero,1)*(1+(facFOVTh-.15))/2);
        yBlockPadZero = dynInd(yBlockPadZero,1:orig,1,0);
    else
        yBlockPadZero = yBlockPad;
    end
    [~, idxPeakPad(1,vI,:) ] = max(abs(yBlockPadZero),[],dimNew);
    if nargout>2; pPadAccum = pPadAccum + multDimSum(abs(yBlockPad),2)/NP(2); end
    
    if averageFreqCoils; idxPeakPad(1,vI,:) = repmat( round(multDimMed(idxPeakPad(1,vI,:),3)),[1 1 size(idxPeakPad,3)]);end
    
    for j=1:length(vI)
        for k = 1:size(idxPeakPad,3)
            pPeakPad(1,vI(j),k) = dynInd(yBlockPadZero,{idxPeakPad(1,vI(j),k),j,k},1:3)/N*(N/NCrop);
        end
    end
    %De-modulate
    t = [0:(N-1)].';
    offset = centerIdx(NCrop)-1 + (croppOffset(1));
    
    ampDet = abs(pPeakPad(1,vI,:));
    freqDet = 1 + (idxPeakPad(1,vI,:)-1)/(NCrop+NPad)*N;
    phaseDet = -angle(pPeakPad(1,vI,:));
    
    yDet = ampDet.*exp(-1i.*   (-2*pi*(freqDet-centerIdx(N))./N.*(t-offset) + phaseDet)   );   

    %Fill original k-space    
    y = dynInd(y,vI,2,yBlock-yDet);
    %Enforce zero
    y(1:croppOffset(1),:) = 0;
    y(end-(croppOffset(2)-1):end,:) = 0;
    
    if deb>1
        %close all
        An = yBlock(:,1,1);
        an = backward(An);
        AnFilt = y(:,1,1);
        anFilt = backward(AnFilt);
        pDet = yDet(:,1,1);
        plotDemodFittingTest(An,an,AnFilt,anFilt,pDet);pause(eps)
    end
end    

if nargout>2;[~,idxPeakPadAverage] = max(pPadAccum,[],dimNew);end

%%% RESHAPE IN CORRECT ARRAY SIZES
NY(end+1:16)=1;
NY2(end+1:16)=1;
%Reshape
%pPeakPad = dynInd(zeros([1 prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(pPeakPad,2,[multDimSum(~idShutter) NY(4:end)]) );
pPeakPad = dynInd(zeros([1 prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(pPeakPad,2,NY2(2:3)));
pPeakPad = resSub(pPeakPad,2,NY(2:3));

%idxPeakPad = dynInd(zeros([1 prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(idxPeakPad,2,[multDimSum(~idShutter) NY(4:end)]) );
idxPeakPad = dynInd(zeros([1 prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(idxPeakPad,2,NY2(2:3)));
idxPeakPad = resSub(idxPeakPad,2,NY(2:3));

%y = dynInd(zeros([N prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(y,2,[multDimSum(~idShutter) NY(4:end)]) );
y = dynInd(zeros([N prod(NY(2:3)) NY(4:end) ],'like',y), ~idShutter,2, resSub(y,2,NY2(2:3)));
y = resSub(y,2,NY(2:3));

%Re-order
pPeakPad = ipermute(pPeakPad,perm);
idxPeakPad = ipermute(idxPeakPad,perm);
y = ipermute(y,perm);

%%% MOVE TO IMAGE FOURIER DOMAIN AGAIN
if isImageDomain
    y = fftshiftOperator(y, 1, 0, dimNew);
    y = fft(y, [], dimNew);
    y = fftshiftOperator(y, 2, 0, dimNew);
end
