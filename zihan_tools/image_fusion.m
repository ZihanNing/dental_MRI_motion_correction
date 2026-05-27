function outFile = image_fusion(caseFolder, seqName, varargin)
%IMAGE_FUSION Fuse MoCo upper/lower reconstructions for one case or many.
%   outFile = image_fusion(caseFolder, seqName) fuses one sequence for one
%   case. If called with no inputs, the legacy batch mode over numeric case
%   folders under rootDir is used.

% by Zihan Ning <zihan.1.ning@kcl.ac.uk>
% @King's College London
% 2026-03-09

parser = inputParser;
parser.addParameter('hfDim', 1);
parser.addParameter('overlapHalfWidth', 5);
parser.addParameter('saveAsSingle', true);
parser.addParameter('verbose', true);
parser.addParameter('rootDir', '/home/zn23/Data/ddMRI/');
parser.parse(varargin{:});
opts = parser.Results;

if nargin == 0
    run_batch_mode(opts);
    outFile = '';
    return;
end

if nargin < 2 || strlength(string(seqName)) == 0
    error('image_fusion requires both caseFolder and seqName for single-case mode.');
end

outFile = fuse_single_case(caseFolder, seqName, opts);
end

function run_batch_mode(opts)
rootDir = opts.rootDir;
d = dir(rootDir);
isCase = [d.isdir] & ~startsWith({d.name}, '.') & ...
         cellfun(@(x) ~isempty(regexp(x, '^\d+$', 'once')), {d.name});
caseDirs = d(isCase);

fprintf('Found %d numeric case folders under %s\n', numel(caseDirs), rootDir);

for iCase = 1:numel(caseDirs)
    caseName = caseDirs(iCase).name;
    caseFolder = fullfile(rootDir, caseName);

    fprintf('\n==============================\n');
    fprintf('Processing case: %s\n', caseName);
    fprintf('Folder: %s\n', caseFolder);

    try
        seqNames = find_available_sequences(caseFolder);
        if isempty(seqNames)
            warning('Case %s: no upper/lower MoCo pair found. Skipping.', caseName);
            continue;
        end

        for iSeq = 1:numel(seqNames)
            outFile = fuse_single_case(caseFolder, seqNames{iSeq}, opts);
            fprintf('Saved fused image: %s\n', outFile);
        end
    catch ME
        warning('Case %s failed: %s', caseName, ME.message);
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

fprintf('\nAll done.\n');
end

function outFile = fuse_single_case(caseFolder, seqName, opts)
anVeDir = fullfile(caseFolder, 'An-Ve');
if ~isfolder(anVeDir)
    error('Folder not found: %s', anVeDir);
end

upperFile = fullfile(anVeDir, [seqName, '_Di_MotCorr_upperjaw_.nii']);
lowerFile = fullfile(anVeDir, [seqName, '_Di_MotCorr_lowerjaw_.nii']);

if ~isfile(upperFile)
    error('MoCo-upper file not found: %s', upperFile);
end
if ~isfile(lowerFile)
    error('MoCo-lower file not found: %s', lowerFile);
end

if opts.verbose
    fprintf('I_U: %s\n', upperFile);
    fprintf('I_L: %s\n', lowerFile);
end

locFile = fullfile(caseFolder, 'location.txt');
if ~isfile(locFile)
    locFile = fullfile(caseFolder, 'location_.txt');
end

if ~isfile(locFile)
    teethMaskFile = findLatestTeethMask(caseFolder);
    if isempty(teethMaskFile)
        error('Neither location.txt nor a full-resolution teeth mask was found for %s.', caseFolder);
    end

    [~, ~, ~, ~, locFile] = compute_landmarks_from_teeth_mask(teethMaskFile, fullfile(caseFolder, 'location.txt'));
    fprintf('Landmarks saved to: %s\n', locFile);
end

infoU = niftiinfo(upperFile);
IU = double(niftiread(infoU));

infoL = niftiinfo(lowerFile);
IL = double(niftiread(infoL));

if ~isequal(size(IU), size(IL))
    error('I_U and I_L have different sizes: %s vs %s', mat2str(size(IU)), mat2str(size(IL)));
end

[idx1, idx2, idx3] = readLocationTxt(locFile);

volSize = size(IU);
maskUpper = false(volSize);
maskLower = false(volSize);

maskUpper = fillMaskBetween(maskUpper, opts.hfDim, idx1, idx2);
maskLower = fillMaskBetween(maskLower, opts.hfDim, idx2, idx3);

maskUpperExp = dilateAlongDim(maskUpper, opts.hfDim, opts.overlapHalfWidth);
maskLowerExp = dilateAlongDim(maskLower, opts.hfDim, opts.overlapHalfWidth);

fixedReg  = IL .* maskUpper;
movingReg = IU .* maskUpper;

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;
optimizer.MinimumStepLength = 1e-5;
optimizer.RelaxationFactor = 0.5;

fixedRegSm  = imgaussfilt3(fixedReg, 1.0);
movingRegSm = imgaussfilt3(movingReg, 1.0);

Rfixed  = imref3d(size(fixedRegSm));
Rmoving = imref3d(size(movingRegSm));

tform_U_to_L = imregtform(movingRegSm, Rmoving, fixedRegSm, Rfixed, ...
                          'rigid', optimizer, metric);

Tinv = inv(tform_U_to_L.T);
Tinv(1:3,4) = 0;
Tinv(4,4)   = 1;
tform_L_to_U = affine3d(Tinv);

if opts.verbose
    fprintf('Rigid registration finished for %s / %s.\n', caseFolder, seqName);
    fprintf('T_{L->U}:\n');
    disp(tform_L_to_U.T);
    fprintf('T_{U->L}:\n');
    disp(tform_U_to_L.T);
end

IU_hat = imwarp(IU, tform_U_to_L, 'OutputView', imref3d(size(IL)), ...
        'Interp', 'nearest', 'FillValues', 0);

upperOnly = maskUpperExp & ~maskLowerExp;
overlap   = maskUpperExp & maskLowerExp;

dU = bwdist(~maskUpperExp);
dL = bwdist(~maskLowerExp);
p = 2;
w = dU.^p ./ (dU.^p + dL.^p + eps);

Ifused = IL;
Ifused(upperOnly) = IU_hat(upperOnly);
Ifused(overlap) = w(overlap) .* IU_hat(overlap) + ...
                  (1 - w(overlap)) .* IL(overlap);

outFile = buildOutputName(caseFolder, upperFile);
outInfo = infoU;

if opts.saveAsSingle
    IfusedToSave = single(Ifused);
else
    IfusedToSave = Ifused;
end

niftiwrite(IfusedToSave, outFile, outInfo, 'Compressed', false);

save(fullfile(anVeDir, [seqName, '_Di_fused_transform.mat']), ...
     'tform_U_to_L', 'tform_L_to_U', ...
     'idx1', 'idx2', 'idx3', ...
     'opts');
end

function seqNames = find_available_sequences(caseFolder)
anVeDir = fullfile(caseFolder, 'An-Ve');
seqNames = {};
if ~isfolder(anVeDir)
    return;
end

upperFiles = dir(fullfile(anVeDir, '*_Di_MotCorr_upperjaw_.nii'));
for k = 1:numel(upperFiles)
    seqName = erase(upperFiles(k).name, '_Di_MotCorr_upperjaw_.nii');
    lowerFile = fullfile(anVeDir, [seqName, '_Di_MotCorr_lowerjaw_.nii']);
    if isfile(lowerFile)
        seqNames{end+1} = seqName; %#ok<AGROW>
    end
end
end

function filePath = findLatestTeethMask(caseFolder)
searchDirs = {
    fullfile(caseFolder, 'An-Aq', 'seg')
    caseFolder
    };

files = [];
for iDir = 1:numel(searchDirs)
    f = dir(fullfile(searchDirs{iDir}, '*msk_teeth_fullresol.nii.gz'));
    if ~isempty(f)
        for iFile = 1:numel(f)
            f(iFile).folder = searchDirs{iDir};
        end
        files = [files; f(:)]; %#ok<AGROW>
    end
end

if isempty(files)
    filePath = '';
    return;
end

[~, idx] = max([files.datenum]);
filePath = fullfile(files(idx).folder, files(idx).name);
end

function [idx1, idx2, idx3] = readLocationTxt(locFile)
vals = readmatrix(locFile, 'FileType', 'text');
vals = vals(~isnan(vals));
if numel(vals) < 3
    error('location.txt does not contain 3 valid indices: %s', locFile);
end
idx1 = round(vals(1));
idx2 = round(vals(2));
idx3 = round(vals(3));
end

function mask = fillMaskBetween(mask, dim, a, b)
lo = min(a, b);
hi = max(a, b);

sz = size(mask);
lo = max(1, min(sz(dim), lo));
hi = max(1, min(sz(dim), hi));

idx = repmat({':'}, 1, ndims(mask));
idx{dim} = lo:hi;
mask(idx{:}) = true;
end

function out = dilateAlongDim(mask, dim, radius)
if radius <= 0
    out = mask;
    return;
end

kernelSize = ones(1, ndims(mask));
kernelSize(dim) = 2 * radius + 1;
kernel = ones(kernelSize);

out = convn(double(mask), kernel, 'same') > 0;
end

function outFile = buildOutputName(caseFolder, upperFile)
[anVeDir, baseName, ~] = fileparts(upperFile);

if endsWith(baseName, '.nii', 'IgnoreCase', true)
    baseName = erase(baseName, '.nii');
end

if contains(baseName, '_Di_MotCorr_upperjaw_')
    outBase = strrep(baseName, '_Di_MotCorr_upperjaw_', '_Di_fused_nearest_');
else
    outBase = [baseName, '_Di_fused_nearest_'];
end

outFile = fullfile(anVeDir, [outBase, '.nii']);
end
