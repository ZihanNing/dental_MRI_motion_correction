
function [fileName] = removeExt(fileName)

%REMOVEEXT removes the extension from a filename. 
%   [FILENAME] = REMOVEEXT(FILENAME)
%   * FILENAME is the filename.
%   ** FILENAME is the filename without the extension.
%       
%   Yannick Brackenier 2023-03-24

[path,file,~] = fileparts(fileName);
fileName=fullfile(path,file);

