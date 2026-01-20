
function TW = dat2TWIX(fileName)

%%% NAME HANDLING
fileName= removeExt(fileName);

%%% CONVERT
if existsFileVar(strcat(fileName,'.dat'))~=1
    TW=[];
    warning('dat2TWIX: %s does not exist. Empty twix object returned.', strcat(fileName,'.dat'));
else
    evalc('TW = mapVBVD(strcat(fileName,''.dat''))');%Suppress command line output
end

