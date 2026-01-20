
function MS = MT2MS(MT)

%MT2MS computes the resolution MS from the orientation information MT.
%
%   [MS]=MT2MS(MT)
%   * MT is the orientation information as the s-form.
%   ** MS is the resolution.
%
%   Yannick Brackenier 2023/07/20

MS = sqrt(sum(MT(1:3,1:3).^2,1));
