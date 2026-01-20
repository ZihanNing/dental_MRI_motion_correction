function mask = dental_mask_edge_aware_bg(vol, varargin)
%DENTAL_MASK_EDGE_AWARE_BG
% Initial edge-aware background removal for dental MRI.
%
% Purpose:
%   - Remove noisy background regions (even if bright)
%   - Preserve anatomy with clear boundaries (mouth, teeth, soft tissue)
%   - Do NOT enforce topology or smoothness (handled by wrapper)
%
% Assumptions:
%   - Single object roughly central
%   - Background touches image corners in axial view
%
% Usage:
%   mask = dental_mask_edge_aware_bg(vol, ...
%       'p_low', 15, 'p_high', 55, ...
%       'thr_tuning', 0.8, 'edge_pct', 85);
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

% ---------------- parameters ----------------
p = inputParser;
addParameter(p, 'p_low', 15);          % lower percentile for background
addParameter(p, 'p_high', 55);         % upper percentile for tissue scale
addParameter(p, 'thr_tuning', 0.8);    % threshold interpolation weight
addParameter(p, 'edge_pct', 85);       % percentile for strong edges
parse(p, varargin{:});
opt = p.Results;

vol = single(vol);
sz  = size(vol);

% ---------------- mild smoothing ----------------
% Suppress voxel-scale noise, preserve anatomy-scale edges
sigma = max(0.6, min(sz)/40);
vol_s = imgaussfilt3(vol, sigma);

% ---------------- adaptive intensity threshold ----------------
nz = vol_s(vol_s > 0);
if isempty(nz)
    mask = true(sz);
    return;
end

pL = prctile(nz, opt.p_low);
pH = prctile(nz, opt.p_high);
thr_int = pL + opt.thr_tuning * (pH - pL);

% ---------------- edge strength (gradient magnitude) ----------------
[gx, gy, gz] = gradient(vol_s);
edge_mag = sqrt(gx.^2 + gy.^2 + gz.^2);

edge_thr = prctile(edge_mag(:), opt.edge_pct);
strong_edge = edge_mag > edge_thr;

% ---------------- background-allowed region ----------------
% Background can flood only where:
%   - intensity is low
%   - edge is weak (no clear anatomical boundary)
bg_allowed = (vol_s < thr_int) | ~strong_edge;

% ---------------- axial corner-seeded background flooding ----------------
bg = false(sz);

for k = 1:sz(3)   % axial slices
    slice = bg_allowed(:,:,k);

    % corner seeds (guaranteed background)
    seed = false(size(slice));
    seed(1,1)     = true;
    seed(1,end)   = true;
    seed(end,1)   = true;
    seed(end,end) = true;

    % morphological reconstruction (flooding)
    bg(:,:,k) = imreconstruct(seed, slice);
end

% ---------------- object = complement of flooded background ----------------
mask = ~bg;

end
