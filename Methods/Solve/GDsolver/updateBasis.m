
function [E] = updateBasis(E, parXB, NX, MT)

%UPDATEBASIS updates the basis functions and headed in the recon structure for a certain resolution level.
%   [E]=UPDATEBASIS(E,PARXB,NX,MT)
%   * E the encoding structure.
%   * PARXB the reconstruction parameters for the B0-informed forward model.
%   * NX the image array size.
%   * MT the orientation information.
%   ** E the updated encoding structure containing the B0 information.
%
%   Yannick Brackenier 2023-07-04

assert(nargin==4,'updateBasis:: Need all inputs to be provided.');

%%% DEAL WITH SCANNER REF FRAME
if parXB.useSH==2; E.Db = E.Dbs;end

%%% DEAL WITH TRANSFORMED BASIS
transformBasis=0;
if isfield(E.Db,'transformed') %we have created basis function in past and have transformed it in LMsolver
    if E.Db.transformed==1; transformBasis=1;T=E.Db.transformHist; end
end

%%%CREATE NEW BASIS
E.Db.SHorder = parXB.SHorder;
E.Db.B = shimBasis(NX, MT, parXB.SHorder);%[E.Db.B, E.Db.Bidx] = SH_basis(NX, parXB.SHorder);
E.Db.B = resPop(E.Db.B,1:3,[],1);
E.Db.B = squeeze(E.Db.B);
E.Db.NX = NX;

%%%TRANSFORM BASIS - so rotation basis functions taken into account
if transformBasis
    Tf = precomputeFactorsSincRigidTransform(E.kG,E.rkG,T,1,0,1,1);%Transform factors 
    E.Db = TransformBasis(E.Db,T,Tf, E.Fof, E.Fob); 
end

%%%RESET TRANFORM PARAMETERS
E.Db.transformed = 0;
E.Db.transformHist = [];

%%% DEAL WITH SCANNER REF FRAME
if parXB.useSH==2; E.Dbs = E.Db;E = rmfield(E,'Db'); end


