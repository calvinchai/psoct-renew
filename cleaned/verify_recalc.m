function verify_recalc
% Recompute biref + orientation from a complex NIfTI, write to temp,
% and verify results match existing outputs.

% --------- Paths you provided ---------
in_complex = "/local_mount/space/megaera/1/users/kchai/code/psoct-renew/vol_recon_test2/rawdata/mosaic_001_image_003_processed_cropped.nii";

ref_biref = "/local_mount/space/megaera/1/users/kchai/code/psoct-renew/vol_recon_test2/process/2d/mosaic_001_image_003_processed_biref.nii";
ref_ori   = "/local_mount/space/megaera/1/users/kchai/code/psoct-renew/vol_recon_test2/process/2d/mosaic_001_image_003_processed_orientation.nii";

% Surface not supplied: use 0 (handled internally)
surfaceFile = "";           % or set to path of a surface NIfTI if you have one

% --------- Parameters (match your pipeline) ---------
depth = 100;                % pixels below surface (depth_lin_ret)
zSize_um = 0.0025;   % derive from NIfTI header (mm -> µm)
lambda_um = 0.0013;           % Telesto ~1300 nm; adjust if needed

% --------- Temp output directory ---------
out_dir = fullfile(tempdir, "psoct_verify");
if ~exist(out_dir, "dir"); mkdir(out_dir); end

out_biref = fullfile(out_dir, "mosaic_001_image_003_processed_biref.nii");
out_ori   = fullfile(out_dir, "mosaic_001_image_003_processed_orientation.nii");

% (We only need biref + orientation here; skip others by giving "")
aip=""; mip=""; ret=""; O3D=""; R3D=""; dBI3D="";

fprintf('Recomputing to temp dir: %s\n', out_dir);

% --------- Run your new function ---------
Complex2Processed( ...
    in_complex, surfaceFile, depth, zSize_um, ...
    aip, mip, ret, out_ori, out_biref, O3D, R3D, dBI3D, ...
    "WavelengthUm", lambda_um);

% --------- Verify: biref ---------
fprintf('\n=== Verifying biref ===\n');
verify_nifti_equal(ref_biref, out_biref, 'biref', ...
    'value_tol', 1e-6, 'rmse_tol', 1e-6);

% --------- Verify: orientation (mod 180°) ---------
fprintf('\n=== Verifying orientation ===\n');
verify_orientation_equal(ref_ori, out_ori, ...
    'abs_tol_deg', 1e-3, 'rmse_tol_deg', 1e-3);

fprintf('\nAll checks completed.\n');
end


function verify_nifti_equal(ref_path, test_path, label, varargin)
    % name-value tolerances
    p = inputParser;
    addParameter(p,'value_tol',1e-6);   % max abs diff
    addParameter(p,'rmse_tol',1e-6);
    parse(p,varargin{:});
    vt = p.Results.value_tol;
    rt = p.Results.rmse_tol;

    [A,IA] = loadNii(ref_path);
    [B,IB] = loadNii(test_path);

    assert(isequal(size(A), size(B)), ...
        '%s size mismatch: %s vs %s', label, mat2str(size(A)), mat2str(size(B)));

    % header geometry sanity (pixdim)
    % assert(roughEq(IA.PixelDimensions, IB.PixelDimensions, 1e-9), ...
    %     '%s PixelDimensions differ: %s vs %s', label, mat2str(IA.PixelDimensions), mat2str(IB.PixelDimensions));

    A = single(A); B = single(B);
    diff = B - A;
    maxAbs = max(abs(diff(:)));
    rmse   = sqrt(mean((diff(:)).^2,'omitnan'));

    fprintf('%s: max|Δ|=%.3g, RMSE=%.3g\n', label, maxAbs, rmse);
    assert(maxAbs <= vt || rmse <= rt, ...
        '%s differences exceed tolerances (max %.3g > %.3g; rmse %.3g > %.3g).', ...
        label, maxAbs, vt, rmse, rt);
end

function verify_orientation_equal(ref_path, test_path, varargin)
    % Compare modulo 180° (orientation is 180°-periodic).
    p = inputParser;
    addParameter(p,'abs_tol_deg',1e-3);
    addParameter(p,'rmse_tol_deg',1e-3);
    parse(p,varargin{:});
    at = p.Results.abs_tol_deg;
    rt = p.Results.rmse_tol_deg;

    [A,IA] = loadNii(ref_path);
    [B,IB] = loadNii(test_path);

    assert(isequal(size(A), size(B)), ...
        'orientation size mismatch: %s vs %s', mat2str(size(A)), mat2str(size(B)));
    % assert(roughEq(IA.PixelDimensions, IB.PixelDimensions, 1e-9), ...
    %     'orientation PixelDimensions differ.');

    A = single(A); B = single(B);
    % wrap absolute difference into [0,90], effectively min(|Δ|, 180-|Δ|)
    d = abs(B - A);
    d = mod(d,180);
    d = min(d, 180 - d);

    maxAbs = max(d(:));
    rmse   = sqrt(mean((d(:)).^2,'omitnan'));
    fprintf('orientation: max|Δ|=%.6g° (mod 180), RMSE=%.6g°\n', maxAbs, rmse);

    assert(maxAbs <= at || rmse <= rt, ...
        'orientation differences exceed tolerances (max %.3g° > %.3g°; rmse %.3g° > %.3g°).', ...
        maxAbs, at, rmse, rt);
end

function tf = roughEq(a,b,tol)
    a = double(a); b = double(b);
    tf = isequal(size(a),size(b)) && all(abs(a(:)-b(:)) <= tol);
end

function [img,info] = loadNii(path)
    if ~isfile(path)
        error('Missing file: %s', path);
    end
    info = niftiinfo(path);
    img = niftiread(info);
end