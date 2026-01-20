
function [d] = ddAbs(x,y)
%double-difference absolute
d = abs( abs(x) - abs(y));