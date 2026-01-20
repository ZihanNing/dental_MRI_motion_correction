function [] = writeM(A, fileName, readType)

if nargin<3 || isempty(readType);readType='w';end%'a' for append

%%% Make sure extension is correct
[path,file] = fileparts(fileName);  
if ~isempty(path)
    fileName = strcat(path,filesep,file,'.m');
else
    fileName = strcat(file,'.m');
end

%%% Open text file
fid = fopen(fileName, readType);

%%% Write
toWrite = sprintf('%s\n', A{:});
fwrite(fid, toWrite);

%%% Close text file
fclose(fid);


% for i = 1:numel(A)
%     if A{i+1} == -1
%         fprintf(fid,'%s', A{i});
%         break
%     else
%         fprintf(fid,'%s\n', A{i});
%     end
% end