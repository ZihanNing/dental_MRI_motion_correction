

function [fPT, fac, fPTOrig, facOrig]= setPTFreq (fCentre, baseResol, voxBW, facOrig, offsetROFOV, resolRO)

%SETPTFREQ determined the Pilot Tone frequency to set for a scanning protocol.
%   [FPT, FAC, FPTORIG, FACORIG]=SETPTFREQ(FCENTRE, BASERESOL, VOXBW, {FACORIG}, {OFFSETROFOV}, {RESOLRO} )
%   * FCENTRE the scanner cenre frequency in Hz.
%   * BASERESOL the base resolution (number of voxels in the FOV).
%   * voxBW the bandwidth (BW) per voxel in the readout in Hz/vox.
%   * {FACORIG} the factor where to put the PT signal in the over-sampled FOV. -1 and 1 for respectively both edges of the over-sampled FOV.
%   * {OFFSETROFOV} the offset in mm of the FOV w.r.t the scanner isocenre.
%   * {RESOLRO} theresolution in mm in the readout.
%   ** FPT is the the PT frequency to set.
%   ** FAC is the adjusted factor.
%   ** FPTORIG is the the PT frequency without re-adjusting FAC (see code).
%   ** FACORIG is the the PT frequency without re-adjusting FAC (see code).
%
%   Yannick Brackenier 2023-07-31

if nargin < 4 || isempty(facOrig); facOrig = .45;end
if nargin < 5 || isempty(offsetROFOV); offsetROFOV = 0;end
if nargin < 6 || isempty(resolRO); resolRO = inf; end

%%% CHECKS
setLowerPart = facOrig<0;
if setLowerPart; fprintf('Pilot Tone signal set to lower part of the RO FOV.\n');end

fac = round(baseResol * facOrig)/baseResol;
if fac ~= facOrig; warning('setPTFreq:: Desired FOV factor not ideal from a DFT leakage point of view.');end

%%% COMPUTE FREQUENCY
mmBW = voxBW/resolRO;%Bandwidth in Hz/mm

fPT = fCentre + ...%In Hz
      mmBW * offsetROFOV + ...
      baseResol * voxBW * fac ;

fPTOrig = fCentre + ...%In Hz
          mmBW * offsetROFOV + ...
          baseResol * voxBW * facOrig ;

