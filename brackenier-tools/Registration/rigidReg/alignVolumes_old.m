
function [x,T] = alignVolumes_old(x, res, ref, exclDim, excludeRef, padding, debug)

%ALIGNVOLUMES   Aligns a set of volumes x to the reference volume x{ref} taking into account a difference in resolution and/or field of view (FOV). 
%   [X,T]=ALIGNVOL(X,{RES},{EXCLDIM},{EXCLUDEREF},{PADDIN},{DEBUG})
%   * X is a cell array with the different volumes to align (w.r.t. x{ref}). 
%     The second row of the cell array are the volumes to apply the transformation to (e.g. in fMRI). Defaults to the original volumes used for registration.
%   * {RES} is a cell array with the respective resolutions.
%   * {REF} the reference volume to align to. Defaults to the 1st.
%   * {EXCLDIM} the dimension to use for FOV extraction during groupwise volume registration
%   * {EXCLUDEREF} to remove the reference volume from the output to save memory. Defaults to 1.
%   * {PADDIN} padding to add to the common space to avoid edges of FOV appearing at other side.
%   * {DEBUG} to show intermediate results and the region to be excluded  during registration.
%   ** X are the aligned volumes, arranged in a multi-dimensional array with the same FOV and resolution as x{ref}. Volume x{2,:} are stacked in the 5th dimension
%   ** T cell of motion parameters to tranfsorm from x{vol} to x{ref} in the common image space (see below).
%

applyToOther = size(x,1)>1; %whether to apply the transformation to additional images as well
if nargin < 2 ; res = []; end 
if nargin < 3 || isempty(ref); ref = 1; end 
if nargin < 4 || isempty(exclDim); exclDim=[];end %Region exclusion will not be used
if nargin < 5 || isempty(excludeRef); excludeRef=0;end
if nargin < 6 || isempty(padding); padding=zeros(1,3);end
if nargin < 7 || isempty(debug); debug=0;end

%%% Make resolutions compatible
if isempty(res); res = cell(1,length(x));end
if length(res)~=size(x,2); res = dynInd(res, length(res)+1:size(x,2),2,{[]});end%cell(1, length(length(res)+1:size(x,2))));end 
if isempty(res{ref}); res{ref}=ones(1,3);end
idx = cellfun(@isempty, res);
if any(idx); res(idx)= { res{ref} }; end %replace empty resolution with the resolution from the 1st volume

%%% Check whether transformation to be applied to other volumes
if applyToOther
    xToTransform = dynInd(x,2,1); 
    x = dynInd(x,1,1); %take out the volumes from x that are used for regisrtation
    %%% Fill empty cell arrays with the volumes in x
    idx = cellfun(@isempty, xToTransform);
    xToTransform = dynInd( xToTransform, find(idx),2, dynInd(x,find(idx),2) );  
end

%%% Make same resolution
N = cellfun(@size, x, 'UniformOutput', false);
nVolumes = size(N,2);
assert( isequal( cellfun(@length, dynInd(N,1,1)), 3*ones(1,nVolumes) ) ,'alignVol:: volumes are expected for all arrays in x.')

for vol = 1:nVolumes 
    if ~isequal(res{ref},res{vol})
        Nv=round(N{vol}.*res{vol}./res{ref});
        x{vol}=resampling(x{vol},Nv);%x{vol} has same resolution as x{1}
        if applyToOther;  xToTransform{vol} = resampling(xToTransform{vol},Nv);end %assumed that xToTransform{vol} has same resolution as x{vol}
    end     
end
Nnew = cellfun(@size, x, 'UniformOutput', false);

%%% Bring to common space
Ncom = multDimMax( cat(4,Nnew{:}), 4 );
Ncom = Ncom + 2*padding;%factor 2 since at both sides
centCom = ceil((Ncom+1)/2);%centre of common field

for vol = 1:nVolumes
    grid = generateGrid(Nnew{vol}, [],Nnew{vol},ceil((Nnew{vol}+1)/2));
    ind = {centCom(1) + grid{1}, centCom(2) + permute(grid{2},[2 1 3]),centCom(3) + permute(grid{3},[3 1 2]) };
    if vol==ref; indRef = ind;end %save to later extract FOV of first volume

    x{vol} = dynInd(zeros(Ncom,'like', x{1}), ind,1:3,x{vol});%Used to be called xcom (common space) but this way memory is saves
    if applyToOther; xToTransform{vol} = dynInd(zeros(Ncom,'like', x{1}), ind,1:3,xToTransform{vol}); end
end

%%% Registration settings
pyr=[8 4 2 1];
tolAccel=0.01;
fraction_order = 0.5;
w = cat(4,x{:});
if applyToOther; wToTransform = cat(4,xToTransform{:});end
    
%%% Region exclusion
if length(exclDim)>1; mask_factor = exclDim(2);else mask_factor = 3;end
NW=size(w);W=w; W=W(:,:,:,1);W(:)=1; %all ones
if ~isempty(exclDim)
    W = dynInd(W, 1:floor(NW(exclDim(1))/mask_factor), exclDim(1),0); % set (neck) to 0 - neck in lower part of third dimension 
    if debug; tt = x{1}.*W; plotND([],abs(tt),[],[],0,{0:2,2}); title('alignVol.m:: Region excluded from rigid registration for ref volume.');end
end            

%%% Register
if debug; plot_(abs(dynInd(w,1,4)),[],[],abs(dynInd(w,1,4)),abs(dynInd(w,2,4)),[],[],[],[],[],200); sgtitle('alignVol.m:: Before alignment for pair 1-2');end
[~,Treg,~,w]=groupwiseVolumeRegistration(abs(w),W,[],0,pyr,fraction_order,[],w,tolAccel); % export Treg for motion estimate
if debug; plot_(abs(dynInd(w,1,4)),[],[],abs(dynInd(w,1,4)),abs(dynInd(w,2,4)),[],[],[],[],[],201);sgtitle('alignVol.m:: After alignment for pair 1-2');end
W=[];

%%%% Transform back to reference volume x{ref}
[~,kGrid,rkGrid,~,~] = generateTransformGrids(size(dynInd(w,1,4)));
[etRef] = precomputeFactorsSincRigidTransform(kGrid,rkGrid,dynInd(Treg,ref,4), 1);       
w = sincRigidTransform(w, etRef, 1, [],[],0); 

if applyToOther
    for vol=1:nVolumes
        di = [0 1];
        Tpair = permute(dynInd(Treg,[vol,ref],4), [1:3 6 4 5]); %compositeTransform.m expects states and parameters in 5th and 6th dimension
        Ttemp{vol} = compositeTransform(gather(Tpair) , di);
    end
    Ttemp = cat(5, Ttemp{:}); Ttemp = permute(Ttemp, [1:3 5 6 4]);
    [etToTr] = precomputeFactorsSincRigidTransform(kGrid,rkGrid,Ttemp, 1); 
    wToTransform = sincRigidTransform(wToTransform, etToTr, 1, [],[],0); 
end

%%% Stack together
if applyToOther; x = cat(5,w, wToTransform);wToTransform=[]; else x = w; end;  w=[];

%%% Extract original FOV x{ref}
x = dynInd(x,indRef,1:3);

%%% Extract reference volume
if excludeRef; x = dynInd(x, setdiff(1:size(x,4),ref) ,4); end

%%% Output transformation parameters
if nargout>1
    T = Ttemp;
    if excludeRef; T = dynInd(T, setdiff(1:size(T,4),ref) ,4); end
    warning('Motion parameters returned might have offset due to difference in FOV.');
end

end
            