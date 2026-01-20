
function [rec, TW] = dat2Rec(fileName, supportReadout, writeNIIFlag, writeRAWFlag, pathOu, isPilotTone, ...
                             removeOversampling, noiseName, resRec, removeZeroSamples, targetR, logFlag, suff)
    
%DAT2REC converts a Siemens .dat file to a reconstruction structure with all the necessary data in it. 
%   [REC , TW] = DAT2REC(FILENAME,SUPPORTREADOUT,{WRITENIIFLAG},{WRITERAWFLAG},{PATHOU},{ISPILOTTONE},{REMOVEOVERSAMPLING},{NOISENAME},{RESREC},{REMOVEZEROSAMPLE}, {TARGETR}, {LOGFLAG}, {SUFF})
%   * FILENAME is the name of the raw data (.dat) of the multiple TE acquisition.
%   * {SUPPORTREADOUT} is the ralative range of the readout FOV to be extracted (between 0 and 1).
%   * {WRITENIIFLAG} s a flag to write the raw coil images to NIFTI files (1). Defaults to 0.
%   * {WRITERAWFLAG} is a flag to write the reconstruction structure to a .mat file.
%   * {PATHOU} is the output path to write all files. Defaults to the directory where the .dat file is stored.
%   * {ISPILOTTONE} indicated whether Pilot Tone was used during this acquisition. For implications, see TWIX2Rec.m.
%   * {REMOVEOVERSAMPLING} indicated whether to remove Siemens oversampling.
%   * {NOISENAME} is the name of the file from which to use noise samples for channel de-correlation.
%   * {RESREC} is the resolution at which to downsample the reconstruction structure.
%   * {REMOVEZEROSAMPLE} is a flag to remove the zero samples when accelerated scans are acquired.
%   * {LOGFLAG} a flag to activate logging.
%   * {SUFF} is the suffix to add to the original filename.
%   ** REC is the reconstruction object.
%   ** TW is the untouched TWIX object.
%       
%   Yannick Brackenier 2023-03-24

if nargin<1 || isempty(fileName); fileName=cd;end%Defaults to current folder
if nargin<2 || isempty(supportReadout); supportReadout=[];end
if nargin<3 || isempty(writeNIIFlag); writeNIIFlag=0;end
if nargin<4 || isempty(writeRAWFlag); writeRAWFlag=0;end
if nargin<5 || isempty(pathOu); pathOu=[];end
if nargin<6 || isempty(isPilotTone); isPilotTone=0;end
if nargin<7 || isempty(removeOversampling); removeOversampling=1;end
if nargin<8 || isempty(noiseName); noiseName='';end
if nargin<9 || isempty(resRec); resRec=[];end
if nargin<10 || isempty(removeZeroSamples); removeZeroSamples=1;end
if nargin<11 || isempty(targetR); targetR=[];end
if nargin<12 || isempty(logFlag); logFlag=1;end
if nargin<13 || isempty(suff); suff='';end

if isfolder(fileName)
    warning(sprintf('dat2Rec:: Function called on all .dat file in the folder %s.\n                   Input parameters are applied to all files!',fileName))
    
    %%% List all .dat files in the folder
    dirInfo = dir(fullfile(fileName,'*.dat') );%Only .dat files
    dirInfo( [dirInfo.isdir]) = []; %Remove directories

    %%% Initialise
    Names = {dirInfo.name};
    Folder = {dirInfo.folder};
    for i=1:length(Names)
        fileNameTemp = fullfile( Folder{i}, Names{i} );
        dat2Rec(fileNameTemp, supportReadout, writeNIIFlag, writeRAWFlag, pathOu, isPilotTone, removeOversampling, noiseName, resRec,removeZeroSamples,targetR,logFlag,suff);
    end
    
    %%% Assign empty arrays to rec and TW as not assigned during loop
    rec=[];
    TW=[];
    
else %We have a file
    %%% NAME HANDLING
    [pathOuTemp,file] = fileparts(fileName);   
    if ~isempty(suff) && ~strcmp(suff(1),'_');suff = strcat('_',suff);end
    file = strcat(file, suff);%Add suffix to the file - considered a new file and not assigned in rec.Names.suff
    if isempty(pathOu); pathOu = pathOuTemp;end
    
    %%% LOGGING
    logFolder = fullfile(pathOu,'Parsing_Log'); if ~isfolder(logFolder);mkdir(logFolder);end
    if logFlag
        logName = strcat(logFolder, filesep, file, '.txt'); 
        if exist(logName, 'file'); delete(logName) ;end
        diary(logName); tStart = tic;
    end

    %%% CONVERT .DAT FILE TO TWIX OBJECT
    fprintf('=====  Converting %s.dat file to a reconstruction structure  ====\n', file);
    c = clock; fprintf('Date of conversion: %d/%d/%d \n', c(1:3));
    fprintf('Converting .dat file to a TWIX object.\n');
    TW = dat2TWIX(fileName);

    %%% CONVERT TWIX TO A REC STRUCTURE
    fprintf('Converting TWIX object to a reconstruction structure.\n');
    if iscell(TW);TW=TW{end};warning('dat2Rec:: Multiple TWIX objects detected. Only last one extracted.');end%Take last one by default as previous ones might be adjustment volumes
    if isfield(TW.hdr.Meas,'tScanningSequence') && contains(TW.hdr.Meas.tScanningSequence,'SE')
        rec = TWIX2recSPACE(TW, supportReadout, isPilotTone, [], removeOversampling, noiseName, fullfile(pathOu,file), resRec,[], removeZeroSamples, targetR);
    else
        rec = TWIX2rec(TW, supportReadout, isPilotTone, [], removeOversampling, noiseName, fullfile(pathOu,file), resRec,[], removeZeroSamples, targetR);        
    end
    rec.Names.Name=file;
    rec.Names.pathOu=pathOu;

    %%% WRITE RAW
    if writeRAWFlag
        fprintf('Saving raw file:\n');
        rec = gatherStruct(rec);
        save(fullfile(pathOu,strcat(file,'.mat')),'rec','-v7.3');    
        fprintf('Raw file saved:\n   %s\n', fullfile(pathOu,strcat(file,'.mat')));
        
        if isfield(rec,'ACS')
           recS = rec2recS(rec);
           ss = [];ss.rec=recS;
           save(fullfile(pathOu,strcat(file,'_ACS.mat')),'-struct','ss','-v7.3'); 
           fprintf('Raw file of ACS data saved:\n   %s\n', fullfile(pathOu,strcat(file,'_ACS.mat')));
        end
    end

    %%% WRITE NIFTI
    if writeNIIFlag>0
        outDir = fullfile(pathOu,'Parsing'); if exist(outDir,'dir')~=2; mkdir(outDir);end
        xW=[];xW{1} = rec.y; 
        MSW=[];MSW{1} = rec.Enc.AcqVoxelSize; 
        MTW =[];MTW{1} = rec.Par.Mine.APhiRec;
        writeNII( strcat(outDir,filesep,rec.Names.Name) , {'y'},xW, MSW, MTW);
        fprintf('NIFTI file with raw coil images saved:\n   %s\n', strcat(outDir,filesep,rec.Names.Name));
        
        % added by ZN for the case of using ACS line
        if isfield(rec,'ACS') % ZN: save the ACS line reconstructed image, only for high-res image
            data = sqrt(sum(square(abs(rec.ACS)),4));
            xW=[];xW{1} = data;
            MSW=[];MSW{1} = rec.Enc.UnderSampling.ACSVoxelSize; 
            MTW =[];MTW{1} = rec.Par.Mine.APhiACS;
            writeNII( strcat(outDir,filesep,rec.Names.Name) , {'y_ACS'},xW, MSW, MTW);
            fprintf('NIFTI file with ACS line recontructed image saved:\n   %s\n', strcat(outDir,filesep,rec.Names.Name));
        end
    end

    %%% WRITE METADATA TO JSON FILE
    writeJSONFlag=1;%Hard-coded to 1
    if writeJSONFlag>0
        recJSON = rec;if nargout<1;rec=[];end
        recJSON.NY=size(recJSON.y);
        %Remove all big arrays
        recJSON.y=[];
        recJSON.S=[];recJSON.N=[];recJSON.Assign=[];
        recJSON.x=[];recJSON.M=[];
        if isfield(recJSON,'ACS');recJSON.ACS=[];end
        recJSON.Par.Labels.Shim.correctFact_mA2au=[];
        if isfield(recJSON.Par,'preProcessing') && isfield(recJSON.Par.preProcessing,'deCorrNoise');recJSON.Par.preProcessing.deCorrNoise.covMatrix=[];end
        if isfield(recJSON,'PT')
            recJSON.PT.yProjRO=[];
            recJSON.PT.pSliceImage=[];
            recJSON.PT.pTimeTest=[];
            recJSON.PT=rmfield(recJSON.PT,'hApod');
        end
        %Remove certain fields
        if isfield(recJSON.Par.Labels,'Shim') && isfield(recJSON.Par.Labels.Shim, 'crossTerms');recJSON.Par.Labels.Shim=rmfield(recJSON.Par.Labels.Shim,'crossTerms');end
        if isfield(recJSON,'Varia');recJSON=rmfield(recJSON,'Varia');end
        if isfield(recJSON.Par,'preProcessing')
            recJSON.Par.preProcessing=rmfield(recJSON.Par.preProcessing,'fftScale');
            recJSON.Par.preProcessing=rmfield(recJSON.Par.preProcessing,'rawCorrectionScale');
            recJSON.Par.preProcessing=rmfield(recJSON.Par.preProcessing,'deCorrNoise');
        end
        %Save JSON
        savejson('',recJSON,sprintf('%s.json',fullfile(pathOu,file)));%writeJSON has specific fields 
        fprintf('JSON file saved:\n   %s\n', fullfile(pathOu,file));
    end
    
    %%% END LOG
    if logFlag; tStop = toc(tStart); fprintf('Data converted in %.0fmin %.0fs.\n',floor(tStop/60),mod(tStop,60)); diary off; end
end
