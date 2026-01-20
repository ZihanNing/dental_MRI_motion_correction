
function [MPS, sliceOrientation] = acquisitionOrder(MT)
    
%ACQUISITIONORDER infers the acquisition order in physical coordinates from an orientation matrix.
%   [MPS , SLICEORIENTATION]=ACQUISITIONORDER(MT)
%   * MT is the orientation matrix in logical coordinates (with the Siemens convention PRS --> RAS).
%   ** MPS is the M(easurement)-P(hase)-S(slice) string containing the physical acquired coordinates in each dimension.
%   ** SLICEORIENTATION is a string containing the slice orientation.
%
%   Yannick Brackenier 2023-04-20
    
%%% Take out translation as this function determines orientation
MT = MT(1:3,1:3);

%%% Define axes in PRS space
vecPE = ([ 1 0 0]');
vecRO = ([ 0 1 0]');
vecSL = ([ 0 0 1]');

%%% Multiply each PRS direction to the RAS frame and deduce the acquistion order
dirPE = MT*vecPE;
dirRO = MT*vecRO;
dirSL = MT*vecSL;

%%% Get MPS
MPS = strcat( dirOrient(dirPE) , '-',...
              dirOrient(dirRO) , '-',...
              dirOrient(dirSL) ) ;

%%% Get slice orientation
if any(contains( 'HF' , MPS(7:8))) || any(contains( 'FH' , MPS(7:8))) %slice information
    sliceOrientation = 'Transversal';
elseif any(contains( 'PA' , MPS(7:8))) || any(contains( 'AP' , MPS(7:8))) %slice information
    sliceOrientation = 'Coronal';
elseif any(contains( 'LR' , MPS(7:8))) || any(contains( 'RL' , MPS(7:8))) %slice information
    sliceOrientation = 'Sagittal';
else
    error('acquisitionOrder:: A slice schould be detected.')
end


end

%%%%%%%%%%%%%%%%%%% HELPER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = dirOrient(vec)

    [ ~ , idx ] = max(abs(vec(1:3)));
    isPos = vec(idx)>0;

    if idx == 1 %Points R(ight)
        str = 'LR';
        if ~isPos; str = fliplr(str);end

    elseif idx ==2 %Points A(nterior)
        str = 'PA';
        if ~isPos; str = fliplr(str);end

    elseif idx ==3 %Points S(uperior)
        str = 'FH';
        if ~isPos; str = fliplr(str);end
    else
        error('One direction should at least be identified.')
    end

end