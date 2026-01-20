
function [h] = makePTApod(N, facFOV, facFOVTh)

%%% SET DIMENSIONS OF BOTH EDGES
facFOV = rescaleND(facFOV,[0 1],[-1 1]);
facFOVTh = rescaleND(facFOVTh,[0 1],[-1 1]);
upperSide = facFOV>0.5;%If PT signal is at the upper side of the (O)FOV

if upperSide
    N1 = round(facFOVTh*N);
    N2 = N - round(facFOV*N); 
else
    N1 = round(facFOV*N);
    N2 = N - round(facFOVTh*N); 
end
Nint = N-N1-N2;

%%% MAKE BOTH ENDS
h1 = ones([N1 1],'single');
h2 = ones([N2 1],'single');

if upperSide
   h2 = 0*h2; 
else
   h1 = 0*h1;  
end
    

%%% INTERPOLATE CURVES
NFix1 = [];NFix2 = [];
offsetFix1=0;offsetFix2=0;
deb=0;
h = interpolteCurveToCurve(h1, h2, Nint, NFix1, NFix2, offsetFix1, offsetFix2, 'cubic',deb);

%%% PLOT
deb=1;
if deb
    figure; 
    plot(h);
end

