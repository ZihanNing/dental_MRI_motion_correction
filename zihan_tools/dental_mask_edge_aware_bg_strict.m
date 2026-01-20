function mask = dental_mask_edge_aware_bg_strict(vol, varargin)
%DENTAL_MASK_EDGE_AWARE_BG_STRICT
% Strict background removal with:
%  - axial single-component + no holes
%  - smooth inter-slice evolution
%  - single 3D component
%  - no holes in axial / coronal / sagittal views
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

% ---------------- parameters ----------------
p = inputParser;
addParameter(p, 'p_low', 15);
addParameter(p, 'p_high', 55);
addParameter(p, 'thr_tuning', 0.75);
addParameter(p, 'edge_pct', 60);
addParameter(p, 'area_jump_ratio', 1.5);
addParameter(p, 'max_iter', 4);
parse(p, varargin{:});
opt = p.Results;

vol = single(vol);
sz  = size(vol);

% ---------------- initial edge-aware background removal ----------------
mask = dental_mask_edge_aware_bg_dynamic(vol, ...
    'p_low', opt.p_low, ...
    'p_high', opt.p_high, ...
    'thr_tuning', opt.thr_tuning, ...
    'edge_pct', opt.edge_pct);

% ---------------- iterative topology enforcement ----------------
for it = 1:opt.max_iter

    mask_prev = mask;

    % ---------- R1: axial single component + no holes ----------
    mask = enforce_axial_single_component(mask);

    % ---------- R2: inter-slice smoothness ----------
    mask = enforce_inter_slice_smoothness(mask, opt.area_jump_ratio);

    % ---------- R3: single 3D component ----------
    mask = keep_central_component(mask);

    % ---------- R4: no holes in coronal + sagittal ----------
    mask = fill_holes_other_views(mask);

    % ---------- re-enforce axial (consistency) ----------
    mask = enforce_axial_single_component(mask);

    % ---------- convergence check ----------
    if isequal(mask, mask_prev)
        break;
    end
end

end


function mask = enforce_axial_single_component(mask)

sz = size(mask);
cx = round(sz(1)/2);
cy = round(sz(2)/2);

for k = 1:sz(3)
    sl = mask(:,:,k);
    if ~any(sl(:)), continue; end

    cc = bwconncomp(sl, 8);

    % choose component closest to centre
    best_i = 1; best_d = inf;
    for i = 1:cc.NumObjects
        [x,y] = ind2sub(size(sl), cc.PixelIdxList{i});
        d = (mean(x)-cx)^2 + (mean(y)-cy)^2;
        if d < best_d
            best_d = d;
            best_i = i;
        end
    end

    sl_clean = false(size(sl));
    sl_clean(cc.PixelIdxList{best_i}) = true;

    % fill 2D holes
    sl_clean = imfill(sl_clean, 'holes');

    mask(:,:,k) = sl_clean;
end
end

function mask = enforce_inter_slice_smoothness(mask, ratio)

areas = squeeze(sum(sum(mask,1),2));
sz = size(mask);

for k = 2:sz(3)-1
    a_med = median([areas(k-1), areas(k), areas(k+1)]);
    if areas(k) > ratio*a_med || areas(k) < a_med/ratio
        mask(:,:,k) = mask(:,:,k) & ...
                      (mask(:,:,k-1) | mask(:,:,k+1));
        mask(:,:,k) = imfill(mask(:,:,k), 'holes');
    end
end
end

function mask = fill_holes_other_views(mask)

% coronal (XZ)
for j = 1:size(mask,2)
    sl = squeeze(mask(:,j,:));
    sl = imfill(sl, 'holes');
    mask(:,j,:) = reshape(sl, [size(mask,1), 1, size(mask,3)]);
end

% sagittal (YZ)
for i = 1:size(mask,1)
    sl = squeeze(mask(i,:,:));
    sl = imfill(sl, 'holes');
    mask(i,:,:) = reshape(sl, [1, size(mask,2), size(mask,3)]);
end
end

function mask_out = keep_central_component(mask)

sz = size(mask);
cx = round(sz(1)/2); cy = round(sz(2)/2); cz = round(sz(3)/2);
center_idx = sub2ind(sz, cx, cy, cz);

cc = bwconncomp(mask, 26);
mask_out = false(sz);

for i = 1:cc.NumObjects
    if any(cc.PixelIdxList{i} == center_idx)
        mask_out(cc.PixelIdxList{i}) = true;
        return;
    end
end

% fallback: centroid closest to centre
center = [cx cy cz];
best_i = 1; best_d = inf;
for i = 1:cc.NumObjects
    [x,y,z] = ind2sub(sz, cc.PixelIdxList{i});
    c = [mean(x) mean(y) mean(z)];
    d = sum((c - center).^2);
    if d < best_d
        best_d = d;
        best_i = i;
    end
end
mask_out(cc.PixelIdxList{best_i}) = true;
end

