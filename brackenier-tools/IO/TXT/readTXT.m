
function A = readTXT(fileName)

%%% Make sure extension is correct
[path,file] = fileparts(fileName);  
if ~isempty(path)
    fileName = strcat(path,filesep,file,'.txt');
else
    fileName = strcat(file,'.txt');
end

%%% Open text file
fid = fopen(fileName,'r');

%%% Loop over lines and read out
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ~feof(fid)%ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end

%%% Close file
fclose(fid);


% formatSpec='%s';
% A = fscanf(fid,formatSpec);