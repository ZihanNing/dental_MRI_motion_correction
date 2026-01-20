function mask_ls = dental_levelset_refine(vol, initMask)
%DENTAL_LEVELSET_REFINE
% Level-set refinement of a dental MRI mask using Chan–Vese model.
%
% vol      : input 3D volume (single)
% initMask: initial binary mask (same size)
%
% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-01-06

vol = single(vol);
initMask = logical(initMask);

% Normalise volume (important)
v = vol;
v = v - prctile(v(:), 1);
v = v / prctile(v(:), 99);
v = max(min(v,1),0);

% Chan–Vese parameters
numIter   = 40;     % keep small
smoothFac = 2.0;    % strong curvature regularisation
contract  = 0;      % no bias to shrink/expand

mask_ls = activecontour( ...
    v, initMask, numIter, ...
    'Chan-Vese', ...
    'SmoothFactor', smoothFac, ...
    'ContractionBias', contract);

end
