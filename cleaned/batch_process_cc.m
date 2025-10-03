function batch_process_cc(inDir, outDir,depth, varargin)
% Batch-process complex PS-OCT volumes to orientation & birefringence NIfTIs (parallel).

if nargin <4
    oriMethod = "";
else
    oriMethod = varargin{1}; % "" or "new"
end
if nargin <5
    birefMethod = "";
else
    birefMethod = varargin{2}; % "" or "new"
end
if nargin < 6
    unwrap = false;
else
    unwrap = varargin{3}; % true or false
end

%inDir  = "/vast/fiber/projects/20250920_CCtest4_30degrees_CW/";
%outDir = "/local_mount/space/megaera/1/users/kchai/project/NewBiref_20250920_CCtest4_30degrees_CW/processed/";

addpath('/local_mount/space/megaera/1/users/kchai/code/psoct-data-processing/vol_recon/');
% ---- Parameters ----
surfaceFile = "";           % If you have per-volume surfaces, make this dynamic per file.
% depth       = 100;          % pixels below surface (depth_lin_ret)
lambda_um   = 0.0013;
zSize_um = 2.5;

% ---- Prep ----
if ~exist(outDir,"dir"), mkdir(outDir); end

files = [ dir(fullfile(inDir, "mosaic_*_image_*_processed_cropped.nii")); ...
    dir(fullfile(inDir, "mosaic_*_image_*_processed_cropped.nii.gz")) ];

if isempty(files)
    error("No input files found in %s", inDir);
end

fprintf("Found %d files.\n", numel(files));

% Ensure a pool exists (adjust workers as you like)
if isempty(gcp('nocreate'))
    parpool('threads');  % or parpool('local', feature('numcores'))
end

% Build per-file paths up-front (parfor-friendly)
N = numel(files);
inPaths   = strings(N,1);
outOri    = strings(N,1);
outBiref  = strings(N,1);

for k = 1:N
    inPath = string(fullfile(files(k).folder, files(k).name));
    [~, base, ext] = fileparts(inPath);
    if strcmpi(ext,'.gz')
        % For .nii.gz, fileparts -> ext=".gz" and base="...nii"
        [base2, ~] = strtok(base, '.');  % base2="..._cropped"
        base = base2;
    end

    [ok, mosaicStr, imageStr] = parse_mosaic_image(base);
    if ~ok
        warning("Skipping unrecognized file name: %s", files(k).name);
        continue;
    end

    inPaths(k)  = inPath;
    outdBI(k)   = fullfile(outDir, sprintf("mosaic_%s_image_%s_processed_dbi.nii",   mosaicStr, imageStr));
    outOri(k)   = fullfile(outDir, sprintf("mosaic_%s_image_%s_processed_ori.nii",   mosaicStr, imageStr));
    outBiref(k) = fullfile(outDir, sprintf("mosaic_%s_image_%s_processed_biref.nii", mosaicStr, imageStr));
    outOriNew(k)   = fullfile(outDir, sprintf("mosaic_%s_image_%s_processed_oriNew.nii",   mosaicStr, imageStr));
    outBirefNew(k) = fullfile(outDir, sprintf("mosaic_%s_image_%s_processed_birefNew.nii", mosaicStr, imageStr));
end

% Remove any empties from bad names
valid = inPaths ~= "";
inPaths  = inPaths(valid);
outOri   = outOri(valid);
outBiref = outBiref(valid);
N = numel(inPaths);

fprintf("Processing %d files in parallel...\n", N);

% ---- Parallel loop ----
parfor k = 1:N
% for k = 1:N
    inPath = inPaths(k);
    % surfaceFile = replace(inPath,'cropped', 'surface_finding')
    % Only produce ori + biref; skip others by passing ""
    aip=""; mip=""; ret=""; O3D=""; R3D=""; dBI3D="";
    ori ="";
    biref = "";
    oriNew = "";
    birefNew = "";
    ori = outOri(k);
    biref = outBiref(k);
    % dBI3D=outdBI(k);
    %oriNew = outOriNew(k);
    birefNew = outBirefNew(k);
    Complex2Processed( ...
        inPath, surfaceFile, depth, zSize_um, ...
        aip, mip, ret, ori, biref,...
        O3D, R3D, dBI3D, oriMethod, birefMethod, unwrap, ...
        "WavelengthUm", lambda_um);

    fprintf("[OK] %s\n", inPath);

end

fprintf("Done. Outputs in: %s\n", outDir);
end

% ================= helpers =================

function [ok, mosaicStr, imageStr] = parse_mosaic_image(base)
% Expect base like "mosaic_002_image_064_processed_cropped"
ok = false; mosaicStr=""; imageStr="";
tokens = regexp(base, '^mosaic_(\d{3})_image_(\d{3})_processed_cropped$', 'tokens', 'once');
if ~isempty(tokens)
    mosaicStr = string(tokens{1});
    imageStr  = string(tokens{2});
    ok = true;
end
end

