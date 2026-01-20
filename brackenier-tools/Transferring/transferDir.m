
function [] = transferDir(type, sourceMach, destMach, sourceDataDir, destDataDir, userName, MBLim, zipRepo)

%TRANFERRDIR  transfers a repo from one machine to the other, with the %option to specify directories and types of data to be transfered. 
%   []=TRANFERRDIR(TYPE,{SOURCEMACH},{DESTMACH},{SOURCEDATADIR},{DESTDATADIR},{USERNAME},{MBLIM})
%   * TYPE is a cell structure with strings of data types to be transfered. Defaults to all.
%   * {SOURCEMACH} is the machine were to copy data from.
%   * {DESTMACH} is the machine were to copy data to.
%   * {SOURCEDATADIR} is the directory where the data is stored on the source machine.
%   * {DESTDATADIR} is the directory where the data is stored on the destiny machine.
%   * {USERNAME} is the username under which login to transfer the data.
%   * {MBLIM} is the limits in MB for which to transfer files.
%   * {ZIPREPO} is a flag to zip the directory.
%
%   Yannick Brackenier 2023-03-27

allTypes = {'Projects', 'Software','Projects-Stanford', 'Software-Stanford'};

if nargin<1 || isempty(type); type = allTypes;end%By default, copy everything
if nargin<2 || isempty(sourceMach); sourceMach = '';end%'gpubeastie05'
if nargin<3 || isempty(destMach); destMach = 'perinatal005-pc';end
if nargin<4 || isempty(sourceDataDir); sourceDataDir = '/home/ybr19/';end
if nargin<5 || isempty(destDataDir); destDataDir = '/home/ybr19/Data/motion-correction/ybr19/Backup';end
if nargin<6 || isempty(userName); userName = 'ybr19';end
if nargin<7 || isempty(MBLim); MBLim = 0.2;end
if nargin<8 || isempty(zipRepo); zipRepo = 1;end

addColonSource = ~strcmp(sourceMach,'');
addColonDest = ~strcmp(destMach,'');

if addColonSource; colonSource = ':';nameSource=strcat(userName,'@');else; colonSource='';nameSource='';end
if addColonDest; colonDest = ':';nameDest=strcat(userName,'@');else; colonDest='';nameDest='';end

if ~isempty(MBLim); limSuff = sprintf('--max-size %.1fmb ',MBLim);else;limSuff='';end

%%% Run over types and transfer
for i = 1:length(type)
    fprintf('\nTransferring %s:\n',type{i});
    suff = type{i};
    
    if zipRepo %%% Create zip version        
        %Create temporary folder
        tempDir = '/home/ybr19/tmp_for_transferring/'; 
        if isfolder(tempDir);warning('transferDir:: Temporary directory exists. Not touching it and aborting.');return;else;mkdir(tempDir);end
        
        %Copy the files to the temporary directory
        command = sprintf('rsync -azP %s%s%s%s%s/%s %s/',...
                   limSuff,...
                   nameSource,sourceMach,colonSource,sourceDataDir,suff,...
                   tempDir);
        dos(command);
        
        %Zip file
        suffDate = getDate();
        ZIPName = sprintf('%s_%s.zip',suffDate,suff);
        command = sprintf('cd %s; zip -r %s %s',tempDir,ZIPName,suff);
        dos(command);
        
        %Send zip file
        limSuffZIP = sprintf('--max-size %dmb ',1000);
        command = sprintf('rsync -azP %s%s%s%s%s/%s %s%s%s%s/',...
                   limSuffZIP,...
                   nameSource,sourceMach,colonSource,tempDir,ZIPName,...
                   nameDest,destMach,colonDest,destDataDir);
        dos(command);            
        
        %Remove temporary folder
        dos(sprintf('rm -r %s',tempDir))
        
    else %%% Copy original directory
        command = sprintf('rsync -azP %s%s%s%s%s/%s %s%s%s%s/',...
                           limSuff,...
                           nameSource,sourceMach,colonSource,sourceDataDir,suff,...
                           nameDest,destMach,colonDest,destDataDir);
        dos(command);            
    end
end
