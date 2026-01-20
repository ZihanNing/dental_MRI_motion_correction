
function [idOu] = extractCoilID(hdr)

id1=[];
id2=[];
id3=[];
id4=[];
id5 = [];
id6 = [];

hdr.Meas.atCoilSelectInfoText

ssPhoenix = hdr.Phoenix.sCoilSelectMeas.aRxCoilSelectData{1, 1}.asList;
ssMeas = hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1, 1}.asList;

for i=1:length(ssPhoenix)
    id1(end+1) = ssPhoenix{1, i}.lADCChannelConnected;
    %id2(end+1) = ssPhoenix{1, i}.lRxChannelConnected;%No good
    id3{end+1} = ssPhoenix{i}.sCoilElementID.tElement;
    id4(end+1) = ssPhoenix{i}.sCoilElementID.ulUniqueKey;
    
    id5(end+1) = ssMeas{i}.lADCChannelConnected;
    id6{end+1} = ssMeas{i}.sCoilElementID.tElement;
end

id3 = cellfun(@(x) str2double(x{1}(2:end)), id3);
id3 = int16(id3);
id6 = cellfun(@(x) str2double(x{1}(2:end)), id6);
id6 = int16(id6);

idOu = id1;