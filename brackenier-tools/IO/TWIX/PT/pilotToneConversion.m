
function [PT] = pilotToneConversion(PT)
 
%PILOTTONECONVERSION sets the Pilot Tone conversion parameters in the rec.PT structure used during reconstruction.
%   [PT]=PILOTTONECONVERSION(PT)
%   * PT is the PT data structure in rec.
%   ** PT is the PT data structure in rec with the filled parameters.
%
%   Yannick Brackenier 2023-03-17

assert(nargin==1,'pilotToneConversion:: PT struct needs to be provided.')

%%% IMAGE K-SPACE DE-MODULATION
PT.deModulate.Flag=1;
PT.deModulate.windowCent=.5;
PT.deModulate.windowWidth=.95;
PT.deModulate.padFac=1;
PT.deModulate.NRep=2;%3 of 20 should work as well on a test dataset
PT.deModulate.facFOVTh=-0.42;
PT.deModulate.applyApod=1;
PT.deModulate.averageFreqCoils=1;

%%% PEAK DETECTION
%Where peak is located
PT.extraction.facFOVTh=PT.deModulate.facFOVTh;

%Peak extraction
PT.extraction.sincInterp.Flag=1;
PT.extraction.sincInterp.windowCent=PT.deModulate.windowCent;
PT.extraction.sincInterp.windowWidth=PT.deModulate.windowWidth;
PT.extraction.sincInterp.padFac=PT.deModulate.padFac;
PT.extraction.sincInterp.averageFreqCoils=PT.deModulate.averageFreqCoils;

PT.extraction.sincInterp.useDemodParam = PT.deModulate.windowCent==PT.extraction.sincInterp.windowCent &&...
                                         PT.deModulate.windowWidth==PT.extraction.sincInterp.windowWidth &&...
                                         PT.deModulate.padFac==PT.extraction.sincInterp.padFac &&...
                                         PT.deModulate.averageFreqCoils==PT.extraction.sincInterp.averageFreqCoils &&...
                                         PT.deModulate.facFOVTh==PT.extraction.facFOVTh;

%Broadband extraction - only used if PT.sincInterp.Flag==0
PT.offsetMBHz = 0;%Multiband in Hz to extract - put this to inf to extract all PT data in the RO
PT.preProcessing.isAveragedMB = 0;

%%% PHASE REFERENCING
PT.preProcessing.isRelativePhase = 0;
PT.preProcessing.relativePhaseCoilIdx = 1;



