
function [A] = createStudy(dirName, writeName, writeMFile)

%CREATESTUDY creates a study file based on files detected in a folder compatible with the reconstruction pipeline.
%   [A]=CREATESTUDY({DIRNAME},{WRITEMFILE},{WRITENAME})
%   * {DIRNAME} is the directory where the raw files are stored.
%   * {WRITENAME} is the name of the study file to write.
%   * {WRITEMFILE} whether to write the study list into a MATLAB script.
%   ** A is the cell array containing all the string lines.
%
%   Yannick Brackenier 2023-04-24

if nargin<1 || isempty(dirName);dirName = cd;end
if nargin<2 || isempty(writeName);writeName = '';end
if nargin<3 || isempty(writeMFile);writeMFile = 1;end

%% PREPARATION
%%% Determine the study name
studyName = strsplit(dirName,'/');
studyName = studyName{end};

%%% List .dat files in folder
dirInfo = dir(strcat(dirName,filesep,'*.dat') );%Only .dat files
dirInfo([dirInfo.isdir]) = []; %Remove directories
Names = {dirInfo.name};
containsRef = contains(Names,'ref') | contains(Names,'REF');
containsAcq = containsRef==0;
NamesRef = Names(containsRef);
NamesAcq = Names(containsAcq);

%% WRITING CELL ARRAY
%%% Set the path names
A = '';
A{end+1} = '';
A{end+1} = sprintf("pathData{1}='/%s/%s/';",dirName,studyName);
A{end+1} = sprintf("pathIn{1}='/%s/%s/';",dirName,studyName);
A{end+1} = sprintf("pathRef{1}='/%s/%s/Re-Se/';",dirName,studyName);
A{end+1} = sprintf("pathRemote{1}='/%s/%s/';",dirName,studyName);
A{end+1} = '';
A{end+1} = '%Initials of HV';
A{end+1} = '';

%%% Experiment details
A{end+1} = '%% EXPERIMENT DETAILS';
A{end+1} = '';

%%% Acquisition names
A{end+1} = '%% ACQUISITIONS DATA';
for i=1:length(NamesAcq)
    [~,NameTemp] = fileparts(NamesAcq{i});
    A{end+1} = sprintf("fileIn{1}{%d}='%s';",i,NameTemp);
end
A{end+1} = '';

%%% Reference names
A{end+1} = '%% REFERENCE DATA';
A{end+1} = "refIn = fillCell(fileIn, '');";
A{end+1} = "refBIn = fillCell(fileIn, '');";
if any(containsRef)
    A{end+1} = '';
    for i=1:length(NamesRef)
        [~,NameTemp] = fileparts(NamesRef{i});
        A{end+1} = sprintf("refIn{1}{%d}='%s';",i,NameTemp);
    end
end
A{end+1} = '';

%%% B0/B1 names
A{end+1} = '%% B0 and B1';
A{end+1} = "B0In = fillCell(fileIn, '');";
A{end+1} = "B1In = fillCell(fileIn, '');";
A{end+1} = '';

%%% Pilot Tone information
A{end+1} = '%% PILOT TONE FLAG';
A{end+1} = "isPT = fillCell(fileIn, 0);%None have PT signal by default.";
A{end+1} = "facFOVPT = fillCell(fileIn, .75);%Where PT signal is placed -1-->1 limit for whole over-sampled FOV.";
A{end+1} = '';

%%% Data conversion
A{end+1} = '%% DATA CONVERSION SPECIFICATIONS';
A{end+1} = 'noiseFile = fillCell(fileIn,refIn{1}{1});';
A{end+1} = 'supportReadout = fillCell(fileIn,[]);';
A{end+1} = 'resRec = fillCell(fileIn,[]);';
A{end+1} = 'RDesired = fillCell(fileIn,[]);';

%% WRITE INTO MATLAB FILE
if writeMFile
    if isempty(writeName)
        temp =  strsplit(studyName,'_');
        writeName = sprintf('/home/ybr19/Projects/Reconstruction/Studies/studies_%s%s%s_%s_AUTOMATED',temp{1},temp{2},temp{3},studyName(12:end));
    else %make sure it's in temp and with date
        [~,writeName] =  fileparts(writeName);
        [~,d,m,y] = getDate();
        writeName = sprintf('/home/ybr19/Projects/Reconstruction/Studies/temp/studies_%s',writeName);
    end
	writeM(A, writeName);
end
