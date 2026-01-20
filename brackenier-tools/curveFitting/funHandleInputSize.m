
function [nelements] = funHandleInputSize(funHandle)

defstr = func2str(funHandle);
test = regexp(defstr, 'x\((\d+)\)', 'tokens');

if isempty(test)
    % Assume we have theta instead of x(1)
    nelements = 1;
else
    nelements = max(str2double([test{:}]));
end
end