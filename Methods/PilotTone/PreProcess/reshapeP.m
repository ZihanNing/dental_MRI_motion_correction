
function [AReshaped, pReshaped, idxToReshape, pTemp] = reshapeP(A, p, type, idx, isOffset)

%RESHAPEP reshapes the Pilot Tone (PT) signal and calibartion data between system models T=pA and T=Ap where T/A/p are respectively the motion parameters/calibration data/PT signal.
%   [ARESHAPED, PRESHAPED, IDXTORESHAPE, PTEMP]=RESHAPEP(A, P,{TYPE},{IDX},ISOFFSET)
%   * A is the calibration data as either a matrix with size [6 NCha] (used in T=Ap) or a "row-flattened" vector with size [NCha*6 1] (used in T=pA).
%   * P is the PT signal as either a matrix with size [NCha NSamples] (used in T=Ap) or a matrix with size [6 NCha*6 NSamples] (used in T=pA with the NSamples dimensions decoupled).
%   * {TYPE} is the type of conversion (w.r.t. the calibration data). Defaults to 'mat2vec', which if the forward conversion T=Ap --> T=pA.
%   * {IDX} are the row indices of A in vector form.
%   * {ISOFFSET} a flag to indicate A includes an offset (see PilotToneAlgorithm.m).
%   ** ARESHAPED is the reshaped calibration data.
%   ** PRESHAPED is the reshaped PT signal.
%   ** IDXTORESHAPE are the row indices used for reshaping A.
%   ** PTEMP is a temporary reshaped version of p to compute the conditioning number.
%
%   Yannick Brackenier 2023-04-06

if nargin<2 || isempty(p);p=[];end
if nargin<3 || isempty(type);type='mat2vec';end
if nargin<4 || isempty(idx);idx=[];end
if nargin<5 || isempty(isOffset);isOffset=0;end%By default

if strcmp(type,'vec2mat') && isempty(idx);error('reshapeP:: Need indices for reshaping from vector to matrix');end

%%% DEDUCE PARAMETER SIZE
if strcmp(type,'mat2vec')
   ndG = size(A,1);
   if ~isempty(p) && ~isempty(A);isOffset= size(A,2)~=size(p,1);end
   if ~isempty(p); NCha = size(p,1)-isOffset;else;NCha = size(A,2)-isOffset;end
elseif strcmp(type,'vec2mat')
   ndG = max(idx);
   NCha = length(idx)/ndG-isOffset;%isOffset must be provided in the input
else
    error('reshapeP:: Type not supported.')   
end

if ~isempty(p)
    if strcmp(type,'mat2vec')
       ndS = size(p,2);
    elseif strcmp(type,'vec2mat')
       ndS = size(p,3);
    end
end

%%% RESHAPE
%Calibration data
if strcmp(type,'mat2vec')
    AReshaped=[];idxToReshape=[];
    for i=1:ndG
        Atemp = A(i,:).';
        AReshaped=cat(1,AReshaped,Atemp);
        idxToReshape = cat(1,idxToReshape, i*onesL(Atemp) );
    end
elseif strcmp(type,'vec2mat')
   AReshaped = zeros([ndG NCha+isOffset],'like',A);
   for i=1:ndG
       Atemp = dynInd(gather(A),idx==i,1).';
       AReshaped = dynInd(AReshaped,i,1,Atemp);
   end 
   idxToReshape=idx;
end

%Pilot Tone signal
if  ~isempty(p)
    if strcmp(type,'mat2vec') 
        pReshaped=zeros([ndG ndG*(NCha+isOffset) ndS], 'like', p);
        if isOffset; p = padArrayND(p,[1 0],[],1,'post');end%Offset in calibration
        for i=1:ndG
            for ii=1:ndS
                pReshaped = dynInd(pReshaped,{i,idxToReshape==i, ii},1:3,dynInd(p,{ii},2).');
            end 
        end
    elseif strcmp(type,'vec2mat') 
        pReshaped=zeros([NCha ndS], 'like', p);
        for i=1:ndS
            pReshaped = dynInd(pReshaped, i, 2, dynInd(p, {1,1:NCha,i},1:3));
        end
    end
else
    pReshaped=[];
end

%------------- temp
if ~isempty(pReshaped) &&  strcmp(type,'mat2vec') 
    assert(ndG==size(pReshaped,1),'temp error')
    pTemp =zeros([size(pReshaped,3)*ndG size(pReshaped,2)],'like',p);
    for i=1:size(pReshaped,3); pTemp= dynInd(pTemp, (1+(i-1)*ndG):i*ndG,1, dynInd(pReshaped,i,3));end
else
    pTemp=[];
end
    
    
    
