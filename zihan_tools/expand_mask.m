function mask_expanded = expand_mask(mask, r)
%EXPAND_MASK_FAST
% Fast, isotropic mask expansion using distance transform.
%
% mask : binary mask
% r    : expansion radius in voxels
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

mask = logical(mask);

% Distance to nearest foreground voxel
d = bwdist(mask);

% Expand
mask_expanded = d <= r;
end
