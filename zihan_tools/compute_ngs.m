function ngs = compute_ngs(I, roiMask, voxSize)
% Compute Normalized Gradient Square (NGS) inside ROI, excluding ROI edges.
%
% ngs = sum(|∇I|^2) / sum(I^2) over voxels that are strictly inside ROI
% (i.e., 6-neighbourhood fully inside), so mask boundaries do not bias gradients.
%
% Inputs:
%   I       : 3D image (single/double)
%   roiMask : logical 3D mask (same size as I)
%   voxSize : [dx dy dz] in mm (used to scale gradient)
%
% Output:
%   ngs     : scalar

if nargin < 3 || isempty(voxSize)
    voxSize = [1 1 1];
end
dx = voxSize(1); dy = voxSize(2); dz = voxSize(3);

roiMask = logical(roiMask);
if ~any(roiMask(:))
    ngs = NaN;
    return;
end

% "Interior" voxels only: remove 1-voxel boundary of the ROI (6-connected erosion)
interior = imerode(roiMask, strel('sphere', 1));
if ~any(interior(:))
    ngs = NaN;
    return;
end

% Compute gradients (scaled by voxel spacing)
% gradient() supports spacing arguments: gradient(F, dx, dy, dz)
[Gx, Gy, Gz] = gradient(I, dx, dy, dz);

g2 = Gx.^2 + Gy.^2 + Gz.^2;

num = sum(g2(interior), 'omitnan');
den = sum((I(interior)).^2, 'omitnan');

if den <= 0 || ~isfinite(den)
    ngs = NaN;
else
    ngs = num / den;
end