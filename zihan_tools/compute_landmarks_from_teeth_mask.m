function [idxHeadTop, idxLipsMid, idxChinBottom, idxRLMid, locPath] = ...
    compute_landmarks_from_teeth_mask(teethMaskPath,locPath)
% Compute HF landmarks from a teeth fullres mask NIfTI (.nii or .nii.gz).
% Landmarks (HF voxel indices, same convention as your manual logging):
%   1) Top of head: always the largest HF index (top of image)
%   2) Middle of lips: in RL mid-slice, midpoint between lowest upper (1) and highest lower (2)
%   3) Bottom of chin: in RL mid-slice, lowest lower (2)
%
% Outputs are also written into location.txt in the same folder as the mask:
%   idxHeadTop
%   idxLipsMid
%   idxChinBottom
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-02-25

%% ---- checks ----
if nargin < 1 || isempty(teethMaskPath)
    error('Input teethMaskPath is empty.');
end
if ~isfile(teethMaskPath)
    error('File not found: %s', teethMaskPath);
end

[pathPrep, ~, ~] = fileparts(teethMaskPath);

%% ---- read nifti ----
info = niftiinfo(teethMaskPath);
V = niftiread(info);
V = uint8(V);

sz = size(V);
if numel(sz) ~= 3
    error('Expect a 3D mask. Got size: [%s]', num2str(sz));
end

%% ---- infer RL/AP/HF dims from affine ----
T = info.Transform.T;     % 4x4
A = T(1:3, 1:3);          % 3x3

worldAxisOfVoxel = zeros(1,3); % 1=X (RL), 2=Y (AP), 3=Z (HF)
for vdim = 1:3
    [~, w] = max(abs(A(:, vdim)));
    worldAxisOfVoxel(vdim) = w;
end

dimRL = find(worldAxisOfVoxel == 1, 1);
dimAP = find(worldAxisOfVoxel == 2, 1); %#ok<NASGU>
dimHF = find(worldAxisOfVoxel == 3, 1);

if isempty(dimRL) || isempty(dimHF)
    error('Failed to infer RL/HF dims from affine. worldAxisOfVoxel=%s', mat2str(worldAxisOfVoxel));
end

%% ---- your explicit convention: head is largest HF index, chin is smallest ----
idxHeadTop = sz(dimHF);

%% ---- upper and lower masks ----
upperMask = (V == 1);
lowerMask = (V == 2);

if ~any(upperMask(:))
    error('Upper teeth label (1) not found in mask.');
end
if ~any(lowerMask(:))
    error('Lower teeth label (2) not found in mask.');
end

%% ---- RL mid-slice based on edges of upper and lower volumes ----
Uproj = any(upperMask, setdiff(1:3, dimRL));
Uvec = squeeze(Uproj);
rlU = find(Uvec);
rlU_min = min(rlU);
rlU_max = max(rlU);

Lproj = any(lowerMask, setdiff(1:3, dimRL));
Lvec = squeeze(Lproj);
rlL = find(Lvec);
rlL_min = min(rlL);
rlL_max = max(rlL);

midU = (rlU_min + rlU_max) / 2;
midL = (rlL_min + rlL_max) / 2;
idxRLMid = round((midU + midL) / 2);
idxRLMid = max(1, min(sz(dimRL), idxRLMid));

%% ---- sagittal slice at RL mid, then compute HF indices ----
idx = {':', ':', ':'};
idx{dimRL} = idxRLMid;
Vsag = squeeze(V(idx{:}));

remainDims = setdiff(1:3, dimRL, 'stable');  % dims preserved in Vsag after squeeze

if remainDims(1) == dimHF
    hfAxisInSag = 1;
elseif remainDims(2) == dimHF
    hfAxisInSag = 2;
else
    error('Internal error: cannot map HF axis in sagittal slice.');
end

if hfAxisInSag == 1
    upperHF_has = any(Vsag == 1, 2);  % rows are HF
    lowerHF_has = any(Vsag == 2, 2);
else
    upperHF_has = any(Vsag == 1, 1);  % cols are HF
    lowerHF_has = any(Vsag == 2, 1);
    upperHF_has = upperHF_has(:);
    lowerHF_has = lowerHF_has(:);
end

hfUpperIdx = find(upperHF_has);

if isempty(hfUpperIdx) 
    error('In RL mid-slice (RL=%d), cannot find both upper label.', idxRLMid);
end

% ---- Lower jaw: recompute from full 3D mask (V), not from Vsag ----
linLowerAll = find(V == 2);
if isempty(linLowerAll)
    error('Cannot find lower-jaw label (==2) anywhere in the 3D mask.');
end

subAll = cell(1,3);
[subAll{:}] = ind2sub(sz, linLowerAll);
hfAllLower = subAll{dimHF};                 % all HF indices where V==2
hfLowerIdx = unique(hfAllLower);            % analogous to "find(lowerHF_has)"

% With your upside-down convention:
% - "towards chin" = smaller HF index
% - "towards head" = larger HF index
upper_lowest_towardsChin  = min(hfUpperIdx);
lower_lowest_towardsChin  = min(hfLowerIdx);
lower_highest_towardsHead = max(hfLowerIdx);


idxLipsMid    = round((double(upper_lowest_towardsChin) + double(lower_highest_towardsHead)) / 2);
idxChinBottom = lower_lowest_towardsChin;

% safety clamp
idxLipsMid    = max(1, min(sz(dimHF), idxLipsMid));
idxChinBottom = max(1, min(sz(dimHF), idxChinBottom));

%% ---- write location.txt ----
fid = fopen(locPath, 'w');
if fid < 0
    error('Cannot write: %s', locPath);
end
fprintf(fid, '%d\n%d\n%d\n', idxHeadTop, idxLipsMid, idxChinBottom);
fclose(fid);

end