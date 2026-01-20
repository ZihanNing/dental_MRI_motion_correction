function mask = dental_mask_edge_aware_bg_dynamic(vol, varargin)
%DENTAL_MASK_EDGE_AWARE_BG_DYNAMIC
% Edge-aware background removal with slice-adaptive, smooth thresholds.
%
% Implements:
%   - per-slice intensity statistics
%   - inter-slice threshold continuity
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

% ---------------- parameters ----------------
p = inputParser;
addParameter(p, 'p_low', 15);
addParameter(p, 'p_high', 55);
addParameter(p, 'thr_tuning', 0.8);
addParameter(p, 'edge_pct', 85);
parse(p, varargin{:});
opt = p.Results;

vol = single(vol);
sz  = size(vol);

%% ---- mild smoothing (3D) ----
sigma = max(0.6, min(sz)/40);
vol_s = imgaussfilt3(vol, sigma);

%% ---- per-slice intensity anchors ----
pL = zeros(sz(3),1);
pH = zeros(sz(3),1);

for k = 1:sz(3)
    sl = vol_s(:,:,k);
    nz = sl(sl > 0);
    if isempty(nz)
        continue;
    end
    pL(k) = prctile(nz, opt.p_low);
    pH(k) = prctile(nz, opt.p_high);
end

%% ---- smooth thresholds along z (continuity prior) ----
pL_s = pL;
pH_s = pH;

for k = 2:sz(3)-1
    pL_s(k) = median([pL(k-1), pL(k), pL(k+1)]);
    pH_s(k) = median([pH(k-1), pH(k), pH(k+1)]);
end

%% ---- edge strength (3D) ----
[gx,gy,gz] = gradient(vol_s);
edge_mag = sqrt(gx.^2 + gy.^2 + gz.^2);
edge_thr = prctile(edge_mag(:), opt.edge_pct);
strong_edge = edge_mag > edge_thr;

%% ---- axial edge-aware background flooding (dynamic thr) ----
bg = false(sz);

for k = 1:sz(3)

    thr_k = pL_s(k) + opt.thr_tuning * (pH_s(k) - pL_s(k));

    bg_allowed = (vol_s(:,:,k) < thr_k) & ~strong_edge(:,:,k);

    seed = false(size(bg_allowed));
    seed(1,1) = true;
    seed(1,end) = true;
    seed(end,1) = true;
    seed(end,end) = true;

    bg(:,:,k) = imreconstruct(seed, bg_allowed);
end

mask = ~bg;

end
