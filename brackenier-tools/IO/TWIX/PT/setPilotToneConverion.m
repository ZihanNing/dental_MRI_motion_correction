
function [flagDemod, flagSinc, facFOVTh, averageCoils] = setPilotToneConverion(flagDemod, flagSinc, facFOVTh, averageCoils)

if nargin<1 || isempty(flagDemod);flagDemod=[];end
if nargin<2 || isempty(flagSinc);flagSinc=[];end
if nargin<3 || isempty(facFOVTh);facFOVTh=[];end
if nargin<4 || isempty(averageCoils);averageCoils=[];end

%%% READ ORIGINAL FILE
A = readM('pilotToneConversion.m');

%%% SET WHAT YOU WANT
if ~isempty(flagDemod)
    id = find(contains(A,'PT.deModulate.Flag')); A{id(1)} = sprintf('PT.deModulate.Flag=%d;',flagDemod);
end
if ~isempty(flagSinc)
    id = find(contains(A,'PT.extraction.sincInterp.Flag')); A{id(1)} = sprintf('PT.extraction.sincInterp.Flag=%d;',flagSinc);
end
if ~isempty(facFOVTh)
    id = find(contains(A,'PT.deModulate.facFOVTh')); A{id(1)} = sprintf('PT.deModulate.facFOVTh=%.2f;',facFOVTh);
end
if ~isempty(averageCoils)
    id = find(contains(A,'PT.deModulate.averageFreqCoils')); A{id(1)} = sprintf('PT.deModulate.averageFreqCoils=%d;',averageCoils);
end

%%% WRITE
fileName = '/home/ybr19/Software/Utilities/IO/TWIX/PT/pilotToneConversion.m';
writeM(A, fileName);

%%% READ OUT FILES
id = find(contains(A,'PT.deModulate.Flag')); flagDemod = findAssignedVal(A{id(1)});
id = find(contains(A,'PT.extraction.sincInterp.Flag')); flagSinc=findAssignedVal(A{id(1)});
id = find(contains(A,'PT.deModulate.facFOVTh')); facFOVTh=findAssignedVal(A{id(1)});
id = find(contains(A,'PT.deModulate.averageFreqCoils')); averageCoils=findAssignedVal(A{id(1)});
