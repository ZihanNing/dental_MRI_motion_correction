
function [CoilSelectMap,fftFactorArray,rawCorrectionFactorArray] = CoilScalingFactors(raw_filename,twix_file)
%[CoilSelectMap] = CoilScalingFactors(raw_filename,Ncha)
%   Returns structure CoilSelectMap with coil scaling factors for the raw 
%   data coil combination as well as arrays fftFactorArray and 
%   rawCorrectionFactorArray that contains those factors sorted by the coil
%   order in twix.image
%   twix_file can be optionally passed: if so then this function won't need
%   to load in again, speeding it up
%   Code refactored from John Asutin Roberts Magnetom post
%   (https://www.magnetom.net/t/fftscale-and-rawdatacorrectionfactor/3330/5)

if nargin==1 || isempty(twix_file)
    twix_file = mapVBVD(raw_filename);
    if numel(twix_file)>1
        twix_file = twix_file{end};
    end
end
Ncha = twix_file.image.NCha;

fid = fopen(raw_filename);
% iMeasHeaderSize = fread(fid,1,'int32'); -> could not make this work, taking 1e6 samples from header:
textHeader = fread(fid, 1e6, 'uchar=>char')'; %iMeasHeaderSize-4 
findThis = '{\s*{\s*{\s*"[^"]+"[^\n]+';
pStart = regexp(textHeader, findThis, 'start');
if ~isempty(pStart)
   pEnd = regexp(textHeader(pStart(1):end), '}[\n\s]*}[\n\s]*}', 'end');
   if ~isempty(pEnd)
      allCoilsInHeaderCell = textHeader(pStart(1):pStart(1)+pEnd(1)+1);
      findThis = '{\s*{\s*"(?<name>[^"]+)"\s*}\s*{\s*(?<fft>[\d\.]+)\s*}\s*{\s*(?<re>[\d\.-]+)\s*}\s*{\s*(?<im>[\d\.-]+)\s*}\s*}';
      CoilStructArray = regexp(allCoilsInHeaderCell, findThis,'names');
      if length(CoilStructArray) == Ncha
         for c=1:Ncha
            tName = CoilStructArray(c).name;
            CoilSelect.fftScale = sscanf(CoilStructArray(c).fft,'%f');
            CoilSelect.rawDataCorrectionFactor = complex(sscanf(CoilStructArray(c).re,'%f'),sscanf(CoilStructArray(c).im,'%f'));
            CoilSelect.txtOrder = c-1;
            CoilSelectMap.(tName) = CoilSelect;
            clear CoilSelect
         end
      end
   end
end
fclose(fid);

%extract coil order in raw data
for c=1:Ncha
    CoilElement{c} = twix_file.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList{1,c}.sCoilElementID.tElement;
end

fftFactorArray = zeros(Ncha,1);
rawCorrectionFactorArray = zeros(Ncha,1);
for c=1:Ncha
    fftFactorArray(c) = CoilSelectMap.(CoilElement{c}).fftScale;
    rawCorrectionFactorArray(c) = CoilSelectMap.(CoilElement{c}).rawDataCorrectionFactor;
end



