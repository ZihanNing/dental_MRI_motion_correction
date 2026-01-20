
function [linsStyle, stars, starsExt] = getStatsVis(pVal)

linsStyle = "-";
if pVal<=1E-3
    stars='***'; 
elseif pVal<=1E-2
    stars='**';
elseif pVal<=0.05
    stars='*';
else 
    stars='';
    linsStyle = '--';
end

starsExt = strcat(stars,sprintf(' p=%1.0E',pVal));
