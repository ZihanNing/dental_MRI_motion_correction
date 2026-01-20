function mask = dental_background_mask(vol, varagin)
%DENTAL_MASK_AXIAL_PRIOR
% Enforces:
%  Please make sure the vol that input is in axial view (as the processing
%  priorize the axial for connected region detection and hole filling)
%
%  1) per-axial-slice: single connected region, no holes
%  2) globally: single 3D region, no holes in any view
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06
if nargin<2 && isempty(varagin);varagin.p_low = 20; varagin.p_high = 55; varagin.thr_tuning = 0.8; end

vol = single(vol);
sz  = size(vol);

%% ---- Step 0: smoothing + adaptive threshold ----
sigma = max(0.6, min(sz)/40);
vol_s = imgaussfilt3(vol, sigma);

nz = vol_s(vol_s > 0);
p_low = prctile(nz, varagin.p_low);
p_high = prctile(nz, varagin.p_high);
thr = p_low + varagin.thr_tuning * (p_high - p_low);

bw = vol_s > thr;

%% ---- Step 1: axial slice enforcement ----
mask = false(sz);
cx = round(sz(1)/2);
cy = round(sz(2)/2);

for k = 1:sz(3)
    slice = bw(:,:,k);

    if ~any(slice(:))
        continue;
    end

    % connected components in this slice
    cc = bwconncomp(slice, 8);

    % pick component closest to slice centre
    best_i = 0;
    best_d = inf;

    for i = 1:cc.NumObjects
        [x,y] = ind2sub(size(slice), cc.PixelIdxList{i});
        d = (mean(x) - cx)^2 + (mean(y) - cy)^2;
        if d < best_d
            best_d = d;
            best_i = i;
        end
    end

    slice_mask = false(size(slice));
    slice_mask(cc.PixelIdxList{best_i}) = true;

    % fill all 2D holes in this slice
    slice_mask = imfill(slice_mask, 'holes');

    mask(:,:,k) = slice_mask;
end

%% ---- Step 2: keep single 3D component (central) ----
mask = keep_component_touching_center(mask);

%% ---- Step 3: remove holes in other views (consistency) ----
for it = 1:2
    mask = fill_holes_coronal_sagittal(mask);
    mask = keep_component_touching_center(mask);
end

end

function mask_out = keep_component_touching_center(bw)
sz = size(bw);
cx = round(sz(1)/2); cy = round(sz(2)/2); cz = round(sz(3)/2);
center_idx = sub2ind(sz, cx, cy, cz);

cc = bwconncomp(bw, 26);
mask_out = false(sz);

for i = 1:cc.NumObjects
    if any(cc.PixelIdxList{i} == center_idx)
        mask_out(cc.PixelIdxList{i}) = true;
        return;
    end
end

% fallback: centroid closest to center
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

function m = fill_holes_coronal_sagittal(m)

% Coronal (XZ)
for j = 1:size(m,2)
    sl = squeeze(m(:,j,:));
    sl = imfill(sl, 'holes');
    m(:,j,:) = reshape(sl, [size(m,1), 1, size(m,3)]);
end

% Sagittal (YZ)
for i = 1:size(m,1)
    sl = squeeze(m(i,:,:));
    sl = imfill(sl, 'holes');
    m(i,:,:) = reshape(sl, [1, size(m,2), size(m,3)]);
end

end

