
function A = readTXT(fileName)

fid = fopen(fileName,'r');
i = 1;
tline = fgetl(fid);
A  {i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A  {i} = tline;
end
fclose(fid);
