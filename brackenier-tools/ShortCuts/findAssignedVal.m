
function [val] = findAssignedVal(str)

id1 = strfind(str,'=');
id1 = id1(end)+1;

id2 = strfind(str,';');
id2 = id2(end)-1;

val = single(str2double(str(id1:id2)));