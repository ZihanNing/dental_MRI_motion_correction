
function [figFlag] = newFigFlag(numFig)

%%% NEWFIGFLAG determines to create a new figure based on the figure number.
%
%   [FIGFLAG] = NEWFIGFLAG(NUMFIG)
%   * NUMFIG is the figure number.
%   ** FIGFLAG the logical flag.
%
%   Yannick Brackenier 2023-04-03


figFlag = ~( isscalar(numFig) && numFig<0);