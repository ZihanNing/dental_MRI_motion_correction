
function [pPeakPad, idxPeakPad, idxPeakPadAverage] = sincInterpolatePT(p, dim, padFac, windowCent, windowWidth, gibbs, facFOV, isImageDomain)

%SINCINTERPOLATEPT interpolates k-space to detect Pilot Tone (PT) signals as the peak amplitude in the image domain in the readout direction.
%   [PPEAKPAD,IDXPEAKPAD,IDXPEAKPADAVERAGE] = SINCINTERPOLATEPT(P,{DIM},{PADFAC},{WINDOWCENT},{WINDOWWIDTH},{GIBBS},{ENFORCEUPPERFOV},{ISIMAGEDOMAIN})
%   * P is the (multi-channel) k-space signal that is contaminated with PT signal.
%   * {DIM} is the readout dimension where to detect PT signal.
%   * {PADFAC} the padding factor w.r.t. the array size to zero-fill k-space.
%   * {WINDOWCENT} the window center that is extracted from the k-space signal (e.g. to extract the readout signal least affected by eddy currents).
%   * {WINDOWWIDTH} the width of the window in relative units (w.r.t the array size).
%   * {GIBBS} the gibbs factor of the window in relative units (w.r.t the array size).
%   * {ENFORCEUPPERFOV} wether to enforce peak detection in the upper part of the array (FOV).
%   * {ISIMAGEDOMAIN} is to indicate if the k-space readout dimension is in the image domain.
%   ** PPEAKPAD is the PT signal as the peak of the interpolated sinc. PPEAKPAD has the same shape as P (apart from the readout dimension).
%   ** IDXPEAKPAD is the detected peak index in the readout dimension.
%   ** IDXPEAKPADAVERAGE is the detected peak index in case all the readouts are averaged.
%
%   Yannick Brackenier 2023-03-20

if nargin<2 || isempty(dim);dim=numDims(p);end 
if nargin<3 || isempty(padFac);padFac=20;end
if nargin<4 || isempty(windowCent);windowCent=[];end
if nargin<5 || isempty(windowWidth);windowWidth=[];end
if nargin<6 || isempty(gibbs);gibbs=0;end%No windowing as interpolation using Fourier domain (~= zero-padding)
if nargin<7 || isempty(facFOV);facFOV=1;end
if nargin<8 || isempty(isImageDomain);isImageDomain=0;end%If the readout is in the image domain (TWIX2rec.m)

%%% PERMUTE AND FLATTEN ARRAY
nD = numDims(p);
perm = 1:nD; perm(perm==dim)=[]; perm=[dim perm];
p = permute(p,perm);
NP = size(p);
p = resSub(p,2:nD);
dimNew = 1;

%Zero-pad parameters
padDim = zerosL(size(p));
padDim(dimNew)=padFac*NP(dimNew);

%%% MOVE TO FOURIER DOMAIN (specific to TWIX2rec.m convention)
if isImageDomain
    p = fftshiftOperator(p, 2, 1, dimNew);
    p = ifft(p, [], dimNew);
    p = fftshiftOperator(p, 1, 1, dimNew);
end

%%% CREATE NEW ARRAYS
NPNew = size(p); NPNew(1)=1; %Only one element extracted
pPeakPad = zeros(NPNew,'like',p);
idxPeakPad = zerosL(pPeakPad);

%%% EXTRACT PART OF THE ARRAY IN THE DIMENSION OF INTEREST (ECC etc.)
if ~isempty(windowWidth) && ~isempty(windowCent)
    NH = onesL(size(p));
    NH(dimNew)=NP(dimNew);
    h=PTROWindow(NH, windowWidth, windowCent, gibbs);
    p=bsxfun(@times, p, h);
	%temp = find(h);
    %croppOffset = [temp(1)-1 NP(1)-temp(end)];
    %croppOffset = max(croppOffset)*onesL(croppOffset);%works better when same pre- and post-cropping
    %p = dynInd(p, (croppOffset(1)+1):(NP(1)-croppOffset(2)),1);
end

%%% LOOP OVER BLOCKS OF ARRAY IN OTHER DIMENSIONS
pPadAccum = 0;
blSz = 5000;
for i=1:blSz:size(p,2)
    %Block index
    vI=i:min(i+blSz-1,size(p,2));
    pBlock = p(:,vI);
    
    fprintf('    Performing sinc-interpolation on samples: %d/%d\n',vI(end),size(p,2));
    
    %Zero-pad in Fourier domain
    pBlockPad = padArrayND(pBlock,round(padDim/2),1,0,'both');

    %%% Move to image domain (specific to TWIX2rec.m convention)
    pBlockPad = fftshiftOperator(pBlockPad, 1, 0, dimNew);
    pBlockPad = fft(pBlockPad, [], dimNew);
    pBlockPad = fftshiftOperator(pBlockPad, 2, 0, dimNew);

    %%% Find peak
    pBlockPadZero = fillCentPT(pBlockPad, dimNew, 0, facFOV);

    [pPeakPad(1,vI), idxPeakPad(1,vI) ] = max(pBlockPadZero,[],dimNew);

    if nargout>2
        pPadAccum = pPadAccum + multDimSum(abs(pBlockPad),2)/NPNew(2);
    end
end

if nargout>2
    [~,idxPeakPadAverage] = max(pPadAccum,[],dimNew);
end

%%% RESHAPE IN CORRECT ARRAY SIZES
%Reshape
pPeakPad = resSub(pPeakPad,2,NP(2:nD));
idxPeakPad = resSub(idxPeakPad,2,NP(2:nD));
%Re-order
pPeakPad = ipermute(pPeakPad,perm);
idxPeakPad = ipermute(idxPeakPad,perm);

