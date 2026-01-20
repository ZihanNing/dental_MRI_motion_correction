
function recS = dat2Se(refName, refBName, supportReadout, writeNIIFlag, writeRAWFlag, pathOu, isPilotTone)

%DAT2SE processes a Siemens .dat file from a low-resolution acquisition and calculates the sensitivity maps using the ESPIRiT algorithm. 
%   [RECS]=DAT2SE(FILENAME,{UPPORREADOUT},{WRITENIIFLAG},{WRITERAWFLAG},{PATHOU},{ISPILOTTONE})
%   * REFNAME is the name of the raw data (.dat) of the reference acquisition used for sensitivity estimation in the image reconstruction step.
%   * {SUPPORTREADOUT} is the ralative range of the readout FOV to be extracted (between 0 and 1).
%   * {WRITENIIFLAG} is a flag to write the B0 map to NIFTI files (1) or to write both maps and raw images (2). Defaults to 1.
%   * {WRITERAWFLAG} is a flag to write the reconstruction structure to a .mat file.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the multiple TE file is stored.
%   * {ISPILOTTONE} indicated whether Pilot Tone was used during this acquisition. This activates a median filter across k-space to detect outliers. Poor performance, so only use when applicable.
%   ** REC is the reconstruction object.
%
%   Yannick Brackenier 2023-03-24

if nargin<2 || isempty(refBName); refBName=[];end

[pathOuTemp,file,suff]=fileparts(refName);
if ~strcmp(refBName,'') && ~isempty(refBName); [~,fileB,~]=fileparts(refBName); else;fileB=[];end

fprintf('=====  Converting %s.dat file to coil sensitivities  ====\n',file);
if ~strcmp(refBName,'') && ~isempty(refBName); fprintf('Using body coil from acquisition: %s\n',fileB);end

if nargin<3 || isempty(supportReadout); supportReadout=[];end
if nargin<4 || isempty(writeNIIFlag); writeNIIFlag=0;end
if nargin<5 || isempty(writeRAWFlag); writeRAWFlag=1;end
if nargin<6 || isempty(pathOu); pathOu=[];end
if nargin<7 || isempty(isPilotTone); isPilotTone=0;end
if isempty(pathOu); pathOu = pathOuTemp;end
logFlag = 1;%~isempty(pathOu);

%%% INVERT SIEMENS DATA
writeNIIFlagRec = (writeNIIFlag>1);%Writing NIFTI of the raw coil images
writeRAWFlagRec = 1;
removeOversampling = 1;
%Array data
if existsFileVar(strcat(pathOuTemp,file), 'rec')==2
    ss=load(strcat(pathOuTemp,file,'.mat'),'rec');%Attention: when recS contains gpuArray variables and you re-use dat2Rec.m, this gives the following error Data no longer exists on the GPU.
    recS = ss.rec;clear ss
else
    recS = dat2Rec(refName, supportReadout, writeNIIFlagRec, writeRAWFlagRec, [], isPilotTone, removeOversampling); 
end
%Body coil data
if ~strcmp(fileB,'') && ~isempty(fileB)
     recB = dat2Rec(refBName, supportReadout, writeNIIFlagRec, writeRAWFlagRec, [], isPilotTone, removeOversampling); 
else
     recB = []; 
end

%%% ACTIVATE LOGGING
if logFlag
    if ~isempty(pathOu); logDir = strcat(pathOu,filesep,'Re-Se_Log');else;logDir='Re-Se_Log';end
    logName =  strcat( logDir,filesep, file, '.txt'); 
    if exist(logName,'file'); delete(logName) ;end; if ~isfolder(logDir); mkdir(logDir);end
    diary(logName); tStart=tic;
    c = clock; fprintf('Date of sensitivity estimation: %d/%d/%d \n', c(1:3));
end

% H=fftshift(buildFilter(size(recS.y,1),'tukeyIso',.8,[],.5));
% H(round(length(H)/2):end)=1;
% figure; plot(H)
% recS.y  = recS.y  .* H;

%%% DEAL WITH MULTI-ECHO DATA
nD = numDims(recS.y);
if nD>4; recS.y = permute(recS.y,[1:4 nD 5:(nD-1) ] );end%Move different echoes from last dimension to the 5th dimension
if size(recS.y,5)>1
    if nD==10;elementToUse=size(recS.y,5);else; elementToUse=1;end
    fprintf('Dealing with muliple images: Taking image %d.\n',elementToUse);
    recS.y = dynInd(recS.y,elementToUse,5);%Extract 1st echo to have minimal phase
end
%TO DO: if multiple echoes, use phase unwrapping and B0 mapping to remove B0 phase and only use transmit/receive phase in ESPIRiT

%%% SOLVE FOR COIL SENSITIVITIES
if (isfield(recS.Enc,'UnderSampling') && recS.Enc.UnderSampling.Flag) || (isfield(recB,'Enc') && isfield(recB.Enc,'UnderSampling') && recB.Enc.UnderSampling.Flag)
    error('dat2Se:: Cannot estimates sensitivities from undersampled data.');
end
fprintf('Estimating coil sensitivities.\n');
recS.Names.pathOu = pathOu;
recS.Names.Name = file;
recS.Plan.Suff=''; recS.Plan.SuffOu='';

doDivEstimation=0
if doDivEstimation
    recS=solveSensit7T_Div(recS, recB);
else
    recS=solveSensit7T(recS, recB);
end

%%% SAVE RAW DATA
if writeRAWFlag
    if existsFileVar(strcat(pathOuTemp,filesep,file),'rec')==2
        recS.y=[];
        %recS.N=[];
        recS.Assign=[];%Remove raw data to avoid storing twice
    end
    save(strcat(pathOuTemp,filesep,file,'.mat'),'recS','-append');%If this file is also used for something else, don't overwrite .mat file    
    fprintf('Raw file saved:\n   %s\n', file);
end

%%% SAVE NIFTI
if writeNIIFlag>0
    fprintf('Writing NIFTI files.\n');
    pathOuNII = strcat(pathOu, '/Re-Se/'); if ~isfolder( pathOuNII);mkdir(pathOuNII);end
    fileSave = strcat( pathOuNII,file);
    xW=[];xW{1} = recS.S; 
    MSW=[];MSW{1} = recS.Enc.AcqVoxelSize; 
    MTW =[]; MTW{1} = recS.Par.Mine.APhiRec;
    writeNII(fileSave, {'Se'},xW, MSW, MTW);
end
if logFlag;tStop = toc(tStart); fprintf('Sensitivities estimated in %.0fmin %.0fs.\n',floor(tStop/60),mod(tStop,60)); diary off;end

end