
function [M] = refineMaskExt(x, Otsu, nErode, nDilate, MS)

if nargin<2 || isempty(Otsu);Otsu=[0:.2:.9];end
if nargin<3 || isempty(nErode);nErode=0;end
if nargin<4 || isempty(nDilate);nDilate=0;end
if nargin<5 || isempty(MS);MS=[1 1 1];end

parS=[];

parS.maskNorm=1;%Norm of the body coil intensities for mask extraction%It was 2
parS.maskTh=1;%Threshold of the body coil intensities for mask extraction%It was 0.2
parS.Otsu=Otsu;%Binary vector of components for multilevel extraction (it picks those components with 1)
parS.nErode=nErode;%Erosion for masking (in mm)
parS.nDilate=nDilate;%Dilation for masking (in mm)
parS.conComp=2;%Whether to get the largest connected component after erosion (when >0). If set to 2, holes in masks will be filled, provided they are not the FOV boundaries. 

M = refineMask(x,parS,MS);
