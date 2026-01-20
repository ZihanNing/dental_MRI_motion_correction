function softMask = soften_mask_distance(mask, width)
%SOFTEN_MASK_DISTANCE
% Create a soft-edged mask using distance transform.
%
% mask  : binary 3D mask
% width : transition width in voxels (e.g. 2–5)
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

mask = logical(mask);

% Distance inside and outside
d_in  = bwdist(~mask);
d_out = bwdist(mask);

% Signed distance (positive inside)
sd = d_in - d_out;

% Smooth step function
softMask = 0.5 * (1 + tanh(sd / width));

% Clamp numerical issues
softMask = max(min(softMask, 1), 0);
end
