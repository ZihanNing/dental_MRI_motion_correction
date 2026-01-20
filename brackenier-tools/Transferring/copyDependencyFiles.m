
function [fList] = copyDependencyFiles(funName, pathOu, dirName, deb)

%COPYDEPENDECNYFILES copies all the files needed to run a certain function/script.
%   [FLIST]=COPYDEPENDECNYFILES(FUNNAME, {DIRNAME}, {PATHOU})
%   * FUNNAME is the name of the function.
%   * {PATHOU} is the path where to create the {DIRNAME}.
%   * {DIRNAME} is the name of the directory where to store the dependency files. If empty, no files are copied.
%   * {DEB} a flag to print debugging information.
%   ** FLIST is the list of dependency functions.
%
%   Yannick Brackenier 2023-07-18

if nargin<2 || isempty(pathOu);pathOu=cd;end
if nargin<3 || isempty(dirName);dirName=[];end
if nargin<4 || isempty(deb);deb=1;end

%%% GET DEPENDENCY LIST
fList = matlab.codetools.requiredFilesAndProducts(funName);

%%% WRITE FILES
if ~isempty(dirName)
    %SET AND CREATE PATHS
    pathOuTemp = fileparts(strcat(filesep,dirName,filesep));   
    dirName = strsplit(pathOuTemp,'/');
    dirName = dirName{end};

    %COPY
    dirWrite = strcat(pathOu,filesep,dirName,filesep);
    mkdir(dirWrite);
    for i=1:length(fList)
       if deb;fprintf('Copying %s \n',fList{i});end
       cmd = sprintf('scp %s %s',fList{i},dirWrite);
       res = dos(cmd); 
       assert(res==0,'copyDependencyFiles:: copying failed for command %s.', cmd)
    end
end

