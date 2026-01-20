
function [E, x] = LMsolverBasis(y,E,x,C,deb)

%LMSOLVER   Performs a Levenger-Marquardt iteration on a series of 
%   parameters of the encoding matrix to minimize a least squares 
%   backprojection of the reconstruction to the measured data
%   E=LMsolver(Y,E,X,{C},{DEB})
%   * Y is the measured data
%   * E is the encoding structure
%   * X is the reconstructed data
%   * {C} is a constrain structure
%   * {DEB} indicates whether to print information about convergence
%   ** E is the modified encoding structure
%
%   Yannick Brackenier 2023-07-14 (modified from original script by Lucilio Cordero-Grande)

if nargin<4;C=[];end
if nargin<5 || isempty(deb);deb=2;end

%GENERAL PARAMETERS AND ARRAYS
gpu=isa(x,'gpuArray');
T = E.Tr;E=rmfield(E,'Tr'); %YB: need to remove it otherwise the encode will perform the tranfsormation twice
if isfield(E,'Dl'); TD=E.Dl.Tr;end
if isfield(E,'Ds'); TD=E.Ds.Tr;end
if isfield(E,'Db'); Db = E.Db; E = rmfield(E,'Db');end%YB: need to remove it otherwise the encode will dephase it twice
if isfield(E,'Dc'); Dc = E.Dc; E = rmfield(E,'Dc');end%YB: need to remove it otherwise the encode will dephase it twice
if isfield(E,'Dl'); Dl = E.Dl; E = rmfield(E,'Dl');end%YB: need to remove it otherwise the encode will dephase it twice
if isfield(E,'Ds'); Ds = E.Ds; E = rmfield(E,'Ds');end%YB: need to remove it otherwise the encode will dephase it twice
if isfield(E,'B0Exp'); B0Exp = E.B0Exp; E = rmfield(E,'B0Exp');end%YB: need to remove it otherwise the encode will dephase it twice
if isfield(E,'Dbs'); Dbs = E.Dbs; E = rmfield(E,'Dbs');end%Could include this part in encode.m, but it will be skipped there as it in in the if E.Tr....end
if isfield(E,'B1m'); B1m = E.B1m; E = rmfield(E,'B1m');end
if isfield(E,'Shim'); Shim = E.Shim; E = rmfield(E,'Shim');end

NT=size(T);
NCr=multDimSize(Db.cr,1:6);

NX=size(x);    
dimG=max(6,length(NCr));
dimS=find(NCr(1:dimG-1)~=1);if isempty(dimS);dimS=dimG-1;end%Dimensions of motion states and parameters of the transform 
ndS=NCr(dimS);
ndG=NCr(dimG);NTS=length(dimS);
ha=horzcat(repmat(1:ndG,[2 1]),nchoosek(1:ndG,2)');%Combinations of derivatives with repetitions to approximate the Hessian
ndH=size(ha,2);   
flagw=zeros(NCr(1:dimS(end)));
En=single(zeros([ndS 1]));EnPrev=En;EnPrevF=EnPrev;EnF=En;   
dH=single(zeros([ndS ndH]));dG=single(zeros([ndS ndG]));dGEff=dG;    %YB: For every state you have 21 Hessian elements and 6 gradient terms    
multA=1.2;multB=2;%Factors to divide/multiply the weight that regularizes the Hessian matrix when E(end)<E(end-1)    
NElY=numel(y);
if isfield(E,'nF');E.Sf=dynInd(E.Sf,E.nF,3);y=dynInd(y,E.nF,3);end%Extract ROI in the readout direction
if isfield(E,'Ps') && ~isempty(E.Ps);E.Sf=bsxfun(@times,E.Sf,E.Ps);y=bsxfun(@times,y,E.Ps);end%Preconditioner of the coils, deprecated

%BLOCK SIZES AND CONVERGENCE VALUES FOR THE TRANSFORM
ET.bS=E.bS;ET.dS=E.dS;ET.oS=E.oS;ET.kG=E.kG;ET.rkG=E.rkG;   
if isfield(E,'Bm');ET.Bm=E.Bm;ET.bS(1)=1;end%YB: in encode.m = slice mask
if isfield(E,'Fms');ET.Sf=E.Sf;ET.Fms=E.Fms;ET.bS(1)=NX(3);E.dS(1)=1;end%YB: in encode.m = Fourier encoding multislice
if isfield(E,'ZSl');ET.ZSl=E.ZSl;E.ZSl=-E.ZSl;end%YB: in encode.m = Slab extraction
ET.cTFull=Db.cT;%YB:B0-modification

if NTS==1;NB=1;else; NB=NT(dimS(1));end
for b=1:NB
    if NB==1;ET.cT=ET.cTFull;else; ET.cT=dynInd(ET.cTFull,b,dimS(1));end
    if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

    %COMPUTE THE JACOBIAN AND STORE VALUES
    ind2EstOr=find(~ET.cT(:));NEst=length(ind2EstOr); %YB: indeces which state to update (so not converged)
    if isfield(E,'Fs');ry=zeros([size(y,1) 1],'like',real(y));else; ry=zeros([NT(6) 1],'like',real(y));end

    if exist('Ds','var') && isfield(Ds,'B0')% Susceptibility model with voxel basis
        x = bsxfun(@times, x , exp(1i*2*pi*Ds.TE*Ds.B0));
    end

    for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2EstOr(vA);
        if isfield(E,'Fms')
           [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
           E.Sf=dynInd(ET.Sf,vA,3);            
        end                               
        if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
        if NB==1;vT=vA;else; vT={b,vA};end
        Ti=dynInd(T,vT,dimS);
        [Tf,~]=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,1,1,1);%Transform factors 

        xT=x;
        if exist('B0Exp','var');xT=bsxfun(@times,xT, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end  
        if exist('Db','var');xT=bsxfun(@times,xT,dephaseBasis( Db.B, dynInd(Db.cr,vA(vA<=E.NMs),5), size(xT),Db.TE));end
        if exist('Dc','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
        if exist('Dl','var');TiD = dynInd(TD,vT,dimS);xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Dl.f,Dl.D,Dl.TE));end% Linear model with voxel model
        if exist('Ds','var') && isfield(Ds,'D');TiD = dynInd(TD,vT,dimS);xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Ds.f,Ds.D,Ds.TE));end% Susceptibility model with voxel basis
        if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xT=bsxfun(@times,xT,intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xT),Db,Dl));end

        if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
        if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                       
        xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);

        if exist('Dbs','var');xT=bsxfun(@times,xT,dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xT),Dbs.TE));end
        if exist('B1m','var');xT=bsxfun(@times,xT,exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xT),[],1)));end   
        if exist('Shim','var');xT=bsxfun(@times,xT,exp(+1i*2*pi*Shim.TE *Shim.B0) );end   
        
        if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end%YB: Taking out E.nF only after transformation, otherwise transformations not consistent with rest of algorithm
        if isfield(E,'Fs')
            E.vA=vA;
            E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
            xT=encode(xT,E)-dynInd(y,indY,1); % YB:This is to create the term w in equation 11 of Lucilio's first paper
        elseif isfield(E,'Bm')
            E.Bm=dynInd(ET.Bm,vA,6);
            xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
        else
            xT=encode(xT,E)-dynInd(y,vT,[5 3]);
        end

        if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
        if isfield(E,'Fs')
            ry(indY)=normm(xT,[],3:4)/NElY;
            EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
        else
            ry(vA)=normm(xT)/NElY;
            EnPrevF=dynInd(EnPrevF,vT,1:NTS,gather(ry(vA)));
        end
        if isfield(E,'Fd') && ~isempty(E.Fd)%Filtering
            xT=fftGPU(xT,3);
            xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
            EnPrev=dynInd(EnPrev,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
        else
            EnPrev=dynInd(EnPrev,vT,1:NTS,dynInd(EnPrevF,vT,1:NTS));
        end

        %xTG=sincRigidTransformGradient(xG,Tf,Tg,E.Fof,E.Fob); % This creates outer 2 operations of w_l in Lucilios 1st paper (still needs encoding!)
        xG=x;
        if exist('B0Exp','var');xG=bsxfun(@times,xG, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end  
        xG=phaseGradient(xG,E,Db,vA,1,NX);   %YB: if you start using vS, thes 2 lines will conflict with each other
        
        for g=1:ndG %YB: Building Gradient (with the use of the Jacobian)
            %Phase in the head frame
            if exist('Dc','var');xG{g}=bsxfun(@times,xG{g},dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
            if exist('Dl','var');TiD = dynInd(TD,vT,dimS);xG{g}=bsxfun(@times,xG{g},dephaseRotationTaylor(TiD, Dl.f,Dl.D,Dl.TE));end% Linear model with voxel model
            if exist('Ds','var') && isfield(Ds,'D');TiD = dynInd(TD,vT,dimS);xG{g}=bsxfun(@times,xG{g},dephaseRotationTaylor(TiD, Ds.f,Ds.D,Ds.TE));end% Susceptibility model with voxel basis
            if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xG{g}=bsxfun(@times,xG{g},intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xG{g}),Db,Dl));end
            
            %Transformation
            xG{g}=sincRigidTransform(xG{g},Tf,1,E.Fof,E.Fob);
            
            %Apply phase in scanner reference to Gradient as well
            if exist('Dbs','var');xG{g}=bsxfun(@times,xG{g},dephaseBasis(Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xG{g}),Dbs.TE));end
            if exist('B1m','var');xG{g}=bsxfun(@times,xG{g},exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xG{g}),[],1)));end   
            if exist('Shim','var');xG{g}=bsxfun(@times,xG{g},exp(+1i*2*pi*Shim.TE*Shim.B0) );end 

            if isfield(E,'nF');xG{g}=dynInd(xG{g},E.nF,3);end
            xG{g}=encode(xG{g},E); % YB: This is the encoding step still needed to become term w_l in 11th equation in 1st paper
            if isfield(E,'Pf') && ~isempty(E.Pf);xG{g}=bsxfun(@times,xG{g},dynInd(E.Pf,indY,1));end  %Precondition using all the residuals, deprecated
            if isfield(E,'Fd') && ~isempty(E.Fd)
                xG{g}=fftGPU(xG{g},3);
                xG{g}=bsxfun(@times,xG{g},dynInd(E.Fd,indY,1));
            end

            if isfield(E,'Fs')
                dG(vA,g)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xG{g}.*conj(xT)),3:4)),vA,1));
            elseif isfield(E,'Bm')
                dG(vA,g)=gather(multDimSum(real(xG{g}.*conj(xT)),1:6));
            else
                dG=dynInd(dG,[vT g],1:NTS+1,gather(permute(multDimSum(real(xG{g}.*conj(xT)),[1:2 4]),[5 3 1 2 4])));
            end
        end
        for h=1:ndH %YB: Building Hessian - approximated by JHJ = different from AlignedRecon paper
            if isfield(E,'Fs')
                dH(vA,h)=gather(dynInd(accumarray(E.mSt(indY),multDimSum(real(xG{ha(1,h)}.*conj(xG{ha(2,h)})),3:4)),vA,1));
            elseif isfield(E,'Bm')
                dH(vA,h)=gather(multDimSum(real(xG{ha(1,h)}.*conj(xG{ha(2,h)})),1:6));
            else
                dH=dynInd(dH,[vT h],1:NTS+1,gather(permute(multDimSum(real(xG{ha(1,h)}.*conj(xG{ha(2,h)})),[1:2 4]),[5 3 1 2 4])));
            end
        end
    end;xG=[];
end

%UPDATE
En=EnPrev;EnF=EnPrevF;
MH=single(eye(ndG));% YB: This is the actual matrix of the Hessian with symmetry
fina=0;
while ~fina   
    %BUILD HESSIAN MATRIX AND POTENTIAL UPDATE
    ry(:)=0;
    for a=find(~ET.cTFull(:) & Db.w(:)<1e10)'            
        aa=ind2subV(NCr(dimS),a);
        for h=1:ndH
            if ha(1,h)==ha(2,h)
                MH(ha(1,h),ha(2,h))=(1+Db.w(a))*dynInd(dH,[aa h],1:NTS+1);
            else
                MH(ha(1,h),ha(2,h))=dynInd(dH,[aa h],1:NTS+1);MH(ha(2,h),ha(1,h))=dynInd(dH,[aa h],1:NTS+1);
            end
        end%YB: Here, the LM parameter w is already included, so no need for identity matrix      
        GJ=double(dynInd(dG,aa,1:NTS));
        GJ=single(double(MH)\GJ(:));% YB: This is where you take the inverse of the Hessian
        GJ=resPop(GJ,1,[],NTS+1);
        dGEff=dynInd(dGEff,aa,1:NTS,-(Db.winit/Db.w(a))*GJ);%YB: this is initialised with 0, so converged shots will have Tupr=0
    end
    cupr=shiftdim(dGEff,-(dimS(1)-1));

    cup = Db.cr+cupr;

    %CHECK ENERGY REDUCTION
    Db.cTFull=(ET.cTFull | flagw);
    for b=1:NB
        if NB==1;Db.cT=Db.cTFull;else; Db.cT=dynInd(Db.cTFull,b,dimS(1));end
        if isfield(E,'Fms');E.Fms=dynInd(ET.Fms,b,5);end

        ind2Est=find(~Db.cT(:));NEst=length(ind2Est);
        for a=1:ET.bS(1):NEst;vA=a:min(a+ET.bS(1)-1,NEst);vA=ind2Est(vA);
            if isfield(E,'Fms')
                [~,E.kG,E.rkG]=generateTransformGrids(NX,gpu,NX,ceil((NX+1)/2),1,[],abs(E.ZSl),vA);
                E.Sf=dynInd(ET.Sf,vA,3);                   
            end                
            if isfield(E,'Fs');indY=[];for c=1:length(vA);indY=horzcat(indY,E.nSt(vA(c))+1:E.nSt(vA(c)+1));end;end
            if NB==1;vT=vA;else vT={b,vA};end
            Ti=dynInd(T,vT,dimS);
            Tf=precomputeFactorsSincRigidTransform(E.kG,E.rkG,Ti,1,0,1,1);

            xT=x;  %For Ds.B0, already computed early on     
            if exist('B0Exp','var');xT=bsxfun(@times,xT, exp(1i* 2*pi*B0Exp.TE * dynInd(B0Exp.B, mod(vA(vA<=E.NMs)-1,8)+1,5)));end
            if exist('Db','var');xT=bsxfun(@times,xT,dephaseBasis( Db.B, dynInd(cup,vA(vA<=E.NMs),5), NX,Db.TE));end
            if exist('Dc','var');xT=bsxfun(@times,xT,dephaseRotation(dynInd(Ti,Dc.d,6),Dc.D));end
            if exist('Dl','var');TiD = dynInd(TupD,vT,dimS);xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Dl.f,Dl.D,Dl.TE)); end% Linear model with voxel model
            if exist('Ds','var') && isfield(Ds,'D'); TiD = dynInd(TD,vT,dimS);xT=bsxfun(@times,xT,dephaseRotationTaylor(TiD, Ds.f,Ds.D,Ds.TE));end% Susceptibility model with voxel basis                
            if isfield(E,'intraVoxDeph') && E.intraVoxDeph && exist('Db','var')&& exist('Dl','var');xT=bsxfun(@times,xT,intravoxDephasing(Ti,dynInd(Db.cr,vA(vA<=E.NMs),5),size(xT),Db,Dl));end

            if isfield(E,'ZSl');xT=extractSlabs(xT,abs(E.ZSl),1,1);end
            if isfield(E,'Fms');xT=dynInd(xT,vA,6);end                
            xT=sincRigidTransform(xT,Tf,1,E.Fof,E.Fob);

            if exist('Dbs','var');xT=bsxfun(@times,xT,dephaseBasis( Dbs.B, dynInd(Dbs.cr,vA(vA<=E.NMs),5), size(xT),Dbs.TE));end
            if exist('B1m','var');xT=bsxfun(@times,xT,exp(dephaseBasis(B1m.B, dynInd(B1m.cr,vA(vA<=E.NMs),5), size(xT),[],1)));end   
            if exist('Shim','var');xT=bsxfun(@times,xT,exp(+1i*2*pi*Shim.TE *Shim.B0) );end   

            if isfield(E,'nF');xT=dynInd(xT,E.nF,3);end            
            if isfield(E,'Fs')
                E.bS(1)=vA(end)-vA(1)+1;E.oS(1)=vA(1);E.dS(1)=vA(end);
                E.vA=vA;
                xT=encode(xT,E)-dynInd(y,indY,1);
            elseif isfield(E,'Bm')
                E.Bm=dynInd(ET.Bm,vA,6);
                xT=encode(xT,E)-bsxfun(@times,E.Bm,y);
            else
                xT=encode(xT,E)-dynInd(y,vT,[5 3]);
            end
            if isfield(E,'Pf') && ~isempty(E.Pf);xT=bsxfun(@times,xT,dynInd(E.Pf,indY,1));end%Precondition using all the residuals, deprecated
            if isfield(E,'Fs')
                ry(indY)=normm(xT,[],3:4)/NElY;
                EnF=dynInd(EnF,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),ry(indY)),vA,1)));
            else
                ry(vA)=normm(xT)/NElY;
                EnF=dynInd(EnF,vT,1:NTS,gather(ry(vA)));                                 
            end
            if isfield(E,'Fd') && ~isempty(E.Fd)
                xT=fftGPU(xT,3);
                xT=bsxfun(@times,xT,dynInd(E.Fd,indY,1));
                En=dynInd(En,vT,1:NTS,gather(dynInd(accumarray(E.mSt(indY),normm(xT,[],3:4)/NElY),vA,1)));
            else
                En=dynInd(En,vT,1:NTS,dynInd(EnF,vT,1:NTS));
            end   
        end
    end
    En(Db.w(:)>1e10)=EnPrev(Db.w(:)>1e10);         
    flagw(En(:)<=EnPrev(:))=1; % YB: flags whether energy reduction in that state, when this is all ones, you quit (fina=1)

    if deb>=2;fprintf('Energy before: %0.6g / Energy after: %0.6g\n',sum(EnPrevF(:)),sum(EnF(:)));end  
    if any(~flagw(:)) % YB: not all converged yet
        Db.w(En(:)>EnPrev(:) & ~ET.cTFull(:))=Db.w(En(:)>EnPrev(:) & ~ET.cTFull(:))*multB;
    else
        Db.w(~ET.cTFull(:))=Db.w(~ET.cTFull(:))/multA;%YB: ~ET.cTFull(:) since this parameter comes from the caller and indicates these states have not been touced
        Db.w(Db.w<1e-8)=multA*Db.w(Db.w<1e-8);%To avoid numeric instabilities 
        Db.cr=cup;
        %if exist('Dl','var');TD = TupD; end%YB
        %if exist('Ds','var') && isfield(Ds, 'D');TD = TupD; end%YB
        fina=1;
        %traDiff=abs(dynInd(cupr,1:3,dimG));
        %rotDiff=abs(convertRotation(dynInd(cupr,4:6,dimG),'rad','deg'));  
        %traDiffMax=multDimMax(traDiff,dimS);rotDiffMax=multDimMax(rotDiff,dimS);    
        %ET.cTFull(max(traDiff,[],dimG)<E.tL(1) & max(rotDiff,[],dimG)<Db.tL(2))=1; %YB: at this point you update the convergence of the motion states
        if deb>=2
        %     fprintf('Maximum change in translation (vox): %s/ rotation (deg): %s\n',sprintf('%0.3f ',traDiffMax(:)),sprintf('%0.3f ',rotDiffMax(:)));
              fprintf('Not converged states: %d of %d\n',prod(ndS)-sum(single(ET.cTFull(:))),prod(ndS));
        end              
    end
end 

%T=constrain(T,C);
E.Tr=T;Db.cT=ET.cTFull;E.bS=ET.bS;E.dS=ET.dS;E.oS=ET.oS;E.kG=ET.kG;E.rkG=ET.rkG;
if exist('Dl','var');Dl.Tr = TD;end%Updated one
if exist('Ds','var') && isfield(Ds, 'D');Ds.Tr = TD;end%Updated one
if exist('B0Exp','var'); E.B0Exp = B0Exp;end%YB add again
if exist('Db','var'); E.Db = Db;end%YB add again
if exist('Dc','var'); E.Dc = Dc;end%YB add again
if exist('Dl','var'); E.Dl = Dl;end%YB add again
if exist('Ds','var'); E.Ds = Ds;end%YB add again
if exist('Dbs','var'); E.Dbs = Dbs;end%YB add again
if exist('B1m','var'); E.B1m = B1m;end%YB add again
if exist('Shim','var'); E.Shim = Shim;end%YB add again

if isfield(E,'Fms');E.Sf=ET.Sf;E.Fms=ET.Fms;end
if isfield(E,'ZSl');E.ZSl=ET.ZSl;end
if isfield(E,'Bm');E.Bm=ET.Bm;end
if isfield(E,'vA');E=rmfield(E,'vA');end
if isfield(E,'Fd');E=rmfield(E,'Fd');end   
Enmin=min(EnPrev,En);
if ~isfield(E,'En');E.En=Enmin;else E.En(ind2EstOr)=Enmin(ind2EstOr);end%To keep a record of the achieved energy

