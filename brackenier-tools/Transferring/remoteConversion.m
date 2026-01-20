
function [] = remoteConversion(strStudies, idFile, doCopy, doConvert, doRemove, typeToConvert, userName, remoteName, serverName, suff)

%REMOTECONVERSION converts Siemens raw data from a remote machine and cleans raw files afterwards to save memory.
%   []=REMOTECONVERSION(STRSTUDIES, {IDFILE},{DOCOPY},{DOCONVERT},{DOREMOVE},{TYPECONVERT},{USERNAME},{REMOTENAME},{SERVERNAME},{SUFF})
%   * STRSTUDIES is a string with the study list that will be evaluated within the function.
%   * {IDFILE} is the id of the file to convert from STRSTUDIES.
%   * {DOCOPY} is a flag to copy the files from the remote machine.
%   * {DOCONVERT} is a flag to convert the raw data (after copying).
%   * {DOREMOVE} is a flag to remove the raw files (after converting).
%   * {TYPECONVERT} is the type of conversion to apply. Currently supported for empty or 'ref' (to process reference acquisitions).
%   * {USERNAME} is the username under which login to transfer the data.
%   * {REMOTENAME} is the remote machine were to copy data from.
%   * {SERVERNAME} is the name of the machine/server running this function on.
%   * {SUFF} the suffix to add to the converted files.
%
%   Yannick Brackenier 2023-03-24

if nargin<3 || isempty(doCopy);doCopy=1;end
if nargin<4 || isempty(doConvert);doConvert=1;end
if nargin<5 || isempty(doRemove);doRemove=1;end
if nargin<6 || isempty(typeToConvert);typeToConvert='';end%'', 'ref'
if nargin<7 || isempty(userName);userName='ybr19';end
if nargin<8 || isempty(remoteName);remoteName='perinatal005-pc';end
if nargin<9 || isempty(serverName);serverName='gpubeastie05';end
if nargin<10 || isempty(suff);suff='';end

%%% LOAD THE STUDY FILES
evalc(strStudies);
p=1;%Path 1

%%% SELECT FILES TO 
if nargin < 2 || isempty(idFile);idFile=1:length(fileIn{p});end
idRef=onesL(idFile);
if strcmp(typeToConvert,'ref')
    idRef=idFile;
    idFile=onesL(idFile);
end
[pathIn, fileIn, refIn, refBIn, B0In, B1In,fileUnique, refUnique, refBUnique, B0Unique, B1Unique, isPTUnique, noiseFileUnique, supportReadoutUnique, resRecUnique, RDesiredUnique, facFOVThUnique] ...
 = extractStudies(pathIn, fileIn, refIn, refBIn, B0In, B1In, p, idFile, idRef,[],[], isPT, noiseFile, supportReadout, resRec, RDesired, facFOVTh);

%%% CHANGE TO REF FILES IF NEEDED
if strcmp(typeToConvert,'ref')
    fileUnique=refUnique;
    isPTUnique=fillCell(isPTUnique,0);%No PT signal for ref scans
    noiseFileUnique=fillCell(noiseFileUnique,refIn{1}{1});%First ref scan
    resRecUnique=fillCell(resRecUnique,[]);
    RDesiredUnique = fillCell(RDesiredUnique,[]);
    supportReadoutUnique = fillCell(supportReadoutUnique,[]);
end

%%% RUN OVER FILES
for f=1:length(fileUnique{p})
    %%% COPY RAW FILES FROM REMOTE MACHINE
    if doCopy
        cmd = {
            sprintf('ssh %s@%s "cd %s;',userName, remoteName, pathRemote{p})
            sprintf(' rsync -azP %s.dat %s@%s:%s;logout"',fileUnique{p}{f},userName, serverName,pathData{p})
        };
        cmd =cat(2,cmd{:});
        dos(cmd);
    end
    
    %%% SET GLOBAL PARAMETERS
    supportReadout = supportReadoutUnique{p}{f};
    resRec = resRecUnique{p}{f};
    noiseFile = strcat(pathIn{p},filesep,noiseFileUnique{p}{f});
    RDesired = RDesiredUnique{p}{f};%If empty, it is going to assume R is consistent with the Over-Sampling (OS) set in the protocol.
    facFOVTh = facFOVThUnique{p}{f};
    
    %%% HANDLE PT CONVERION
    [~, ~, facFOVThOrig] = setPilotToneConverion();
    setPilotToneConverion([], [], facFOVTh);
       
    %%% FILENAME
    fileName=strcat(pathIn{p},filesep,fileUnique{p}{f});
    
    %%% CONVERT
    if doConvert && ~strcmp(fileName,'')
        dat2Rec(fileName,supportReadout,[],1,[],isPTUnique{p}{f},[],noiseFile,resRec,[],RDesired,[],suff);
    end        

    %%% RE-SET PT CONVERSION
    setPilotToneConverion([], [], facFOVThOrig);
       
    %%% DELETE DAT FILE
    if doRemove
        cmd = sprintf('rm %s/%s.dat',pathData{p},fileUnique{p}{f});
        dos(cmd);
    end
end    

