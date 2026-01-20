function [] = findStringInStruct(x, str, printVar, inputPath)

if nargin<3;printVar=0;end
if nargin<4;inputPath = '';end

if isstruct(x)
    y=fieldnames(x);
    for n=1:length(y)
        if contains(y{n},str)
            fprintf([inputPath '.' y{n} '\n']);
            if printVar
                x.(y{n})
            end
        end
        inputPathTemp = strcat(inputPath, '.',y{n});
        findStringInStruct(x.(y{n}),str,printVar,inputPathTemp);
    end
    
elseif iscell(x)
        for n=1:length(x);findStringInStruct(x{n},str,printVar,inputPath);end
        
elseif isnumeric(x)
    %nothing
    
end
