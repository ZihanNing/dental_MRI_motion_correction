

function [xG] = phaseGradient(x,E, Db, vA, vS,NX)

NB = size(Db.B);

xG = [];
for i=1:NB(2)
    xG{i} = bsxfun(@times, x, dephaseBasis(Db.B, dynInd(Db.cr,vA(vA<=E.NMs),5),NX, dynInd(Db.TE,vS,7) )  );
    xG{i} = bsxfun(@times, xG{i}, reshape(dynInd(Db.B,i,2),NX));
    xG{i} = bsxfun(@times, xG{i}, 1i*2*pi*Db.TE);
end