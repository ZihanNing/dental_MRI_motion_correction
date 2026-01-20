
function [] = transferData(strContain, type, sourceMach, destMach, sourceDataDir,destDataDir, userName, MBLim)

%TRANFERDATA  transfers data from one machine to the other, with the
%option to specify directories and types of data to be transfered. Note
%that this uses the Linux command line so might not work on Windows and
%might require manual password entry (more than once if multiple types
%copied). For the latter issue, consider setting up ssh keys. 
%   []=TRANFERDATA(STRCONTAIN, {TYPE},{SOURCEMACH},{DESTMACH},{SOURCEDATADIR},{DESTDATADIR},{USERNAME},{MBLIM})
%   * STRDIR is a string to is contained in the directory (relative to the dataDir) to be transfered.
%   * {TYPE} is a cell structure with strings of data types to be transfered. Defaults to all.
%   * {SOURCEMACH} is the machine were to copy data from.
%   * {DESTMACH} is the machine were to copy data to.
%   * {SOURCEDATADIR} is the directory where the data is stored on the source machine.
%   * {DESTDATADIR} is the directory where the data is stored on the destiny machine.
%   * {USERNAME} is the username under which login to transfer the data.
%   * {MBLIM} is the limits in MB for which to transfer files.
%
%   Yannick Brackenier 2023-03-24

allTypes = {'An-Ve','An-Ve_Log','An-Ve_Sn',...
            'Parsing_Log','Parsing_Sn',...
            'Re-Se','Re-Se_Log','Re-Se_Sn',...
            'ISD-PACS',...
            'rawmat','rawdat','rawmeta'};

if nargin<2 || isempty(type); type = allTypes;end%By default, copy everything
if nargin<3 || isempty(sourceMach); sourceMach = '';end%'gpubeastie05'
if nargin<4 || isempty(destMach); destMach = 'perinatal005-pc';end
if nargin<5 || isempty(sourceDataDir); sourceDataDir = '/home/ybr19/Data/';end
if nargin<6 || isempty(destDataDir); destDataDir = '/home/ybr19/Data/motion-correction/ybr19/Data_Recon/';end
if nargin<7 || isempty(userName); userName = 'ybr19';end
if nargin<8 || isempty(MBLim); MBLim = [];end

addColonSource = ~strcmp(sourceMach,'');
addColonDest = ~strcmp(destMach,'');

if addColonSource; colonSource = ':';nameSource=strcat(userName,'@');else; colonSource='';nameSource='';end
if addColonDest; colonDest = ':';nameDest=strcat(userName,'@');else; colonDest='';nameDest='';end

flag = cellfun( @(x) strcmp(x,'raw'), type); 
if any(flag); type(flag==1)=[]; type(end+1:end+2) = {'rawmat','rawdat'};end

%%% Extract the directory to transfer
dirInfo = dir(fullfile( sourceDataDir) );
dirInfo( ~[dirInfo.isdir]) = []; %Keep directories

Names = {dirInfo.name};

flag = contains( Names, strContain) ;
if multDimSum(flag)>1
    warning('transferData:: Multiple directories detected. Aborted data transfer.');
    return;
elseif multDimSum(flag)<1
    warning('transferData:: No directories detected. Aborted data transfer.');
    return;
end

idx = find(flag);
dirName = Names{idx(1)}(1:end);%Only first to transfer

%%% Run over types and transfer
for i = 1:length(type)
    fprintf('\nTransferring %s:\n',type{i});
    
    if strcmp(type{i},'rawmat')%Matlab converted raw data
        suff = '*.mat';
    elseif strcmp(type{i},'rawdat') %Siemens raw data
        suff = '*.dat';
    elseif strcmp(type{i},'rawmeta')%JSON metedata
        suff = '*.json';
    else 
        suff = type{i};
    end

    if ~isempty(MBLim); limSuff = sprintf('--max-size %dmb ',MBLim);else;limSuff='';end
    
    command = sprintf('rsync -azP %s%s%s%s%s%s/%s %s%s%s%s%s/',...
                       limSuff,...
                       nameSource,sourceMach,colonSource,sourceDataDir,dirName,suff,...
                       nameDest,destMach,colonDest,destDataDir,dirName);
    dos(command);    
end

