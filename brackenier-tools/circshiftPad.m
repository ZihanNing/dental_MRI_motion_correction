

function [x] = circshiftPad(x, shift)

x = circshift(x,shift);

for i=1:length(shift)
    if shift(i)>0
        x = dynInd(x, 1:shift(i),i,0);
    else
        x = dynInd(x, (size(x,i)+shift(i)+1):size(x,i) ,i,0);
    end
end
%Test
