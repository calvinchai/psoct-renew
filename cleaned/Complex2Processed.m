
function Complex2Processed(input, surfaceFile, depth, zSize, aip, mip, ret, ori, biref, O3D, R3D, dBI3D, oriNew, birefNew, varargin)
% Complex2Processed
% Compute 3D metrics (dBI3D, R3D, O3D) and optional enface 2D maps (AIP, MIP, RET, ORI),
% plus optional birefringence from a single complex PS-OCT NIfTI volume.
%
% INPUTS (use "" to skip any output you don't want):
%   input        : string, path to input complex NIfTI. Layout: 4*X × Y × Z
%                  where blocks along dim-1 are [J1_real; J1_imag; J2_real; J2_imag].
%   surfaceFile  : string, path to NIfTI with surface z-index (per X,Y). If "", assume 0.
%                  Indices are assumed 0-based; internally converted to MATLAB 1-based.
%   depth        : scalar nonneg int, #pixels below surface for enface window and biref fit.
%   zSize        : scalar >0, axial voxel size in micrometers; ONLY needed if biref output path is non-empty.
%   aip          : string, path to write AIP (mean dBI over [surf, surf+depth]). NIfTI.
%   mip          : string, path to write MIP (max dBI over [surf, surf+depth]). NIfTI.
%   ret          : string, path to write RET (mean R3D over [surf, surf+depth]) in degrees. NIfTI.
%   ori          : string, path to write ORI (circular mean O3D over [surf, surf+depth]) in degrees [0,180). NIfTI.
%   biref        : string, path to write birefringence Δn (slope of OPD vs depth) over window; NIfTI.
%   O3D          : string, path to write full 3D optic axis orientation (deg), flipped along z to match original script. NIfTI.
%   R3D          : string, path to write full 3D retardance (deg), flipped along z. NIfTI.
%   dBI3D        : string, path to write full 3D backscatter in dB, flipped along z. NIfTI.
%
% NAME-VALUE (optional):
%   "WavelengthUm" : wavelength in micrometers for birefringence fit of retardance (default 1.3).
%                    Only used if 'biref' path is provided.
%
% Example:
% Complex2Processed("complex.nii.gz","surf.nii.gz",100,3.3,"aip.nii.gz","", "", "", "biref.nii.gz", "", "", "");

    arguments
        input string
        surfaceFile string = ""
        depth (1,1) {mustBeInteger, mustBeNonnegative} = 100
        zSize (1,1) double {mustBePositive} = 1.0
        aip string = ""
        mip string = ""
        ret string = ""
        ori string = ""
        biref string = ""
        O3D string = ""
        R3D string = ""
        dBI3D string = ""
        oriNew string = ""
        birefNew string = ""
    end
    arguments (Repeating)
        varargin
    end
    % Parse optional wavelength for birefringence (micrometers)
    p = inputParser;
    addParameter(p,"WavelengthUm",1.3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
    parse(p,varargin{:});
    lambda_um = p.Results.WavelengthUm;

    % -------- Load input complex NIfTI --------
    infoIn = niftiinfo(input);
    V = niftiread(infoIn);               % single or double; dimensions: (4*X) × Y × Z
    V = single(V);                       % enforce single precision downstream

    % Infer dimensions and split Jones components
    X4 = size(V,1);
    if mod(X4,4)~=0
        error('Input first dimension must be a multiple of 4: got %d.', X4);
    end
    nx = X4/4; ny = size(V,2); nz = size(V,3);

    J1r = V(1:nx,         :, :);
    J1i = V(nx+1:2*nx,    :, :);
    J2r = V(2*nx+1:3*nx,  :, :);
    J2i = V(3*nx+1:4*nx,  :, :);

    J1  = complex(J1r, J1i);
    J2  = complex(J2r, J2i);

    clear V J1r J1i J2r J2i

    % -------- Core 3D metrics (flip along z to match your original code) --------
    IJones = abs(J1).^2 + abs(J2).^2;
    dBI3D_vol = flip( 10*log10(max(IJones, eps('single'))), 3 );

    R3D_vol  = flip( atan( abs(J1)./max(abs(J2), eps('single')) )/pi*180, 3 );

    phase1 = angle(J1);
    phase2 = angle(J2);
    phi = (phase1 - phase2);                       % wrap to [-pi,pi]
    phi(phi >  pi) = phi(phi >  pi) - 2*pi;
    phi(phi < -pi) = phi(phi < -pi) + 2*pi;
    O3D_vol = flip( (phi/(2*pi))*180, 3 );         % degrees, nominally [-90,90]

    % -------- Surface handling for enface & biref windows --------
    if strlength(surfaceFile) > 0
        surf = single(niftiread(surfaceFile));     % expected X × Y of 0-based z-indices
        if ~isequal(size(surf), [nx, ny])
            error('Surface size mismatch: expected %dx%d, got %s.', nx, ny, mat2str(size(surf)));
        end
        surf = surf + 1;                           % convert to MATLAB 1-based indexing
    else
        surf = ones(nx,ny,'single');               % "0" → index 1 in MATLAB
    end

    % Clamp surfaces inside [1, nz]
    surf = max(1, min(nz, round(surf)));

    % Compute per-(x,y) window stop index
    stopIdx = min(nz, surf + depth);

    % -------- Write requested 3D outputs --------
    writeIfPath(dBI3D, dBI3D_vol, infoIn);
    writeIfPath(R3D,   R3D_vol,   infoIn);
    writeIfPath(O3D,   O3D_vol,   infoIn);

    % -------- Enface 2D maps over [surf : surf+depth] --------
    
        % Pre-allocate slices per XY on demand
        if strlength(aip)>0 || strlength(mip)>0
            % AIP/MIP use dBI
            [aipMap, mipMap] = enfaceStat(dBI3D_vol, surf, stopIdx, ...
                                          strlength(aip)>0, strlength(mip)>0);
            writeIfPath(aip, aipMap, shrinkHeader(infoIn));
            writeIfPath(mip, mipMap, shrinkHeader(infoIn));
        end

        if strlength(ret)>0
            retMap = enfaceMean(R3D_vol, surf, stopIdx);      % degrees
            writeIfPath(ret, retMap, shrinkHeader(infoIn));
        end

        if strlength(ori)>0
            % oriMap = enfaceOrientation(O3D_vol, surf, stopIdx); % degrees in [0,180)
            oriMap = orien_enface(O3D_vol,5);
            writeIfPath(ori, oriMap, shrinkHeader(infoIn));
        end

        if strlength(biref)>0
            if ~(isscalar(zSize) && zSize>0)
                error('zSize must be provided and > 0 (micrometers) to compute biref.');
            end
            birefMap = fitBirefringence(R3D_vol, surf, stopIdx, zSize, lambda_um);
            writeIfPath(biref, birefMap, shrinkHeader(infoIn));
        end

        [Q,U,V,norms] = stokes_from_ori_ret(O3D_vol, R3D_vol, 3);
        [ori_full, ori_0_60, ori_0_100, ori_0_140, ori_0_180] = ...
            enface_orientation_stokes(O3D_vol, R3D_vol, Q, U, V, dBI3D_vol, surf, ...
                              'StartOffset',10, 'MaxLen',180, ...
                              'CopyBackN',20, 'LowRetGate',pi/8);
        if strlength(oriNew)>0
            writeIfPath(oriNew, single(ori_full), shrinkHeader(infoIn));
        end
        if strlength(birefNew)>0
            birefMap = fitBirefringenceNew(R3D_vol, dBI3D_vol, surf, stopIdx,zSize,  lambda_um);
            writeIfPath(birefNew, birefMap, shrinkHeader(infoIn));
        end
    
end

% ================== Helpers ==================
function writeIfPath(pathStr, data, infoLike)
    if strlength(pathStr) == 0 || isempty(pathStr), return; end

    % ---- Normalize class (NIfTI has no logical) ----
    if islogical(data)
        data = uint8(data);  % 0/1
    end

    % ---- Build output header from input ----
    infoOut = infoLike;

    % Dimension fixups
    sz = size(data);
    nDims = ndims(data);
    infoOut.ImageSize = sz;

    % PixelDimensions length must match dims
    if isfield(infoOut, 'PixelDimensions')
        pd = infoOut.PixelDimensions;
        if numel(pd) < nDims, pd(end+1:nDims) = 1; end
        if numel(pd) > nDims, pd = pd(1:nDims); end
        infoOut.PixelDimensions = pd;
    end

    % ---- Datatype/BitsPerPixel must match class(data) ----
    [dtype, bpp, dtcode] = class2niftiMeta(class(data));
    infoOut.Datatype     = dtype;
    infoOut.BitsPerPixel = bpp;

    % Raw fields
    if ~isfield(infoOut,'Raw'), infoOut.Raw = struct; end
    if ~isfield(infoOut.Raw,'dim'),    infoOut.Raw.dim    = ones(1,8); end
    if ~isfield(infoOut.Raw,'pixdim'), infoOut.Raw.pixdim = ones(1,8); end

    infoOut.Raw.dim(:) = 1;
    infoOut.Raw.dim(1) = nDims;                 % number of dims
    infoOut.Raw.dim(2:1+numel(sz)) = sz;

    infoOut.Raw.pixdim(:) = 1;
    if isfield(infoOut,'PixelDimensions') && ~isempty(infoOut.PixelDimensions)
        infoOut.Raw.pixdim(2:1+numel(infoOut.PixelDimensions)) = infoOut.PixelDimensions;
    end

    infoOut.Raw.datatype = dtcode;              % NIfTI datatype code (e.g., 16 for FLOAT32)
    infoOut.Raw.bitpix   = bpp;

    % Write
    isgz = endsWith(string(pathStr), ".nii.gz");
    niftiwrite(data, pathStr, infoOut, 'Compressed', isgz);
end

% ---- Helper: map MATLAB class -> NIfTI meta ----
function [dtype, bpp, dtcode] = class2niftiMeta(cls)
    switch cls
        case 'uint8',   dtype='uint8';   bpp=8;  dtcode=2;
        case 'int16',   dtype='int16';   bpp=16; dtcode=4;
        case 'int32',   dtype='int32';   bpp=32; dtcode=8;
        case 'single',  dtype='single';  bpp=32; dtcode=16;   % FLOAT32
        case 'double',  dtype='double';  bpp=64; dtcode=64;   % FLOAT64
        case 'int8',    dtype='int8';    bpp=8;  dtcode=256;
        case 'uint16',  dtype='uint16';  bpp=16; dtcode=512;
        case 'uint32',  dtype='uint32';  bpp=32; dtcode=768;
        case 'int64',   dtype='int64';   bpp=64; dtcode=1024;
        case 'uint64',  dtype='uint64';  bpp=64; dtcode=1280;
        otherwise
            error('Unsupported data class for NIfTI: %s', cls);
    end
end

function info2 = shrinkHeader(infoIn)
    % Make a 2D-compatible header derived from the input header.
    info2 = infoIn;
    info2.ImageSize = infoIn.ImageSize(1:2);
    if isfield(info2,'PixelDimensions') && numel(info2.PixelDimensions)>=2
        info2.PixelDimensions = info2.PixelDimensions(1:2);
    end
    if isfield(info2,'Raw') && isfield(info2.Raw,'dim')
        info2.Raw.dim(1) = 2;          % number of dims
        info2.Raw.dim(2) = infoIn.ImageSize(1);
        info2.Raw.dim(3) = infoIn.ImageSize(2);
        info2.Raw.dim(4) = 1;
        info2.Raw.pixdim(1) = 1;
        info2.Raw.pixdim(2) = info2.PixelDimensions(1);
        info2.Raw.pixdim(3) = info2.PixelDimensions(2);
        info2.Raw.pixdim(4) = 1;
    end
end

function [aipMap, mipMap] = enfaceStat(vol, surf, stopIdx, doAIP, doMIP)
    nx = size(vol,1); ny = size(vol,2);
    if doAIP, aipMap = zeros(nx,ny,'single'); else, aipMap = []; end
    if doMIP, mipMap = zeros(nx,ny,'single'); else, mipMap = []; end

    for x = 1:nx
        for y = 1:ny
            z1 = surf(x,y);
            z2 = stopIdx(x,y);
            if z2 < z1, z2 = z1; end
            slice = vol(x,y,z1:z2);
            if doAIP
                aipMap(x,y) = mean(slice,'omitnan');
            end
            if doMIP
                mipMap(x,y) = max(slice,[],'omitnan');
            end
        end
    end
end

function retMap = enfaceMean(volDeg, surf, stopIdx)
    nx = size(volDeg,1); ny = size(volDeg,2);
    retMap = zeros(nx,ny,'single');
    for x = 1:nx
        for y = 1:ny
            z1 = surf(x,y); z2 = stopIdx(x,y);
            if z2 < z1, z2 = z1; end
            retMap(x,y) = mean( volDeg(x,y,z1:z2), 'omitnan' );
        end
    end
end

function oriMap = enfaceOrientation(O3D_deg, surf, stopIdx)
    % Circular mean of 180°-periodic orientations using doubled-angle trick.
    nx = size(O3D_deg,1); ny = size(O3D_deg,2);
    oriMap = zeros(nx,ny,'single');

    for x = 1:nx
        for y = 1:ny
            z1 = surf(x,y); z2 = stopIdx(x,y);
            if z2 < z1, z2 = z1; end
            theta = squeeze(O3D_deg(x,y,z1:z2));          % degrees (can be [-90,90])
            theta = mod(theta,180);                        % map into [0,180)
            z = exp(1i*deg2rad(2*theta));
            m = angle(mean(z,'omitnan'))/2;               % radians
            if isnan(m), m = 0; end
            oriMap(x,y) = mod(rad2deg(m),180);
        end
    end
end

function birefMap = fitBirefringence(R3D_deg, surf, stopIdx, zSize_um, lambda_um)
    % Fit slope of OPD (cycles*lambda) vs depth to estimate Δn.
    % R3D_deg is retardance in degrees; convert to cycles: R/360.
    % OPD (um) = (R/360)*lambda_um. Slope(OPD vs depth) ≈ Δn.
    nx = size(R3D_deg,1); ny = size(R3D_deg,2);
    birefMap = zeros(nx,ny,'single');

    for x = 1:nx
        for y = 1:ny
            z1 = surf(x,y); z2 = stopIdx(x,y);
            if z2 < z1, z2 = z1; end
            rp = single(squeeze(R3D_deg(x,y,z1:z2)));
            % Stop early if zeros (e.g., crop) appear
            zZeros = find(rp==0,1,'first');
            if ~isempty(zZeros)
                rp = rp(1:zZeros);
                z2 = z1 + numel(rp) - 1;
            end
            cycles = rp / 360;
            OPD = cycles * lambda_um;                 % micrometers
            depth_um = zSize_um * (0:numel(OPD)-1)';  % relative to z1

            if numel(OPD) >= 3 && any(OPD) && any(depth_um)
                p = polyfit(double(depth_um), double(OPD), 1);  % slope ~ Δn (unitless)
                birefMap(x,y) = single(p(1));
            else
                birefMap(x,y) = 0;
            end
        end
    end
end

function birefMap = fitBirefringenceNew(R3D_deg, inten, surf, stopIdx,zSize_um, lambda_um)
    % Fit slope of OPD (cycles*lambda) vs depth to estimate Δn.
    % R3D_deg is retardance in degrees; convert to cycles: R/360.
    % OPD (um) = (R/360)*lambda_um. Slope(OPD vs depth) ≈ Δn.
    [nx, ny, nz] = size(R3D_deg);
    birefMap = zeros(nx,ny,'single');
    RET = R3D_deg / 180 * pi;
    depth_per_pixel = zSize_um;
    for i = 1:nx
        for j = 1:ny
            z1 = surf(i,j); z2 = stopIdx(i,j);
            depth = z2-z1;
            if z2 < z1, z2 = z1; end
            remaining = sum(inten(i,j,:) > 0.01) - surf(i,j);
            if remaining > depth
                base_len = 40;
                iter_len = 60;
                endind = base_len + iter_len;
            elseif remaining > 40
                base_len = 40;
                iter_len = remaining - 40;
                endind = base_len + iter_len - 1;
            elseif remaining > 0
                base_len = remaining;
                iter_len = 5;
                endind = base_len - 1;
            else
                biref(i,j) = 0;
                continue;
            end
            start_idx = max(1, z1);
            end_idx = min(nz, z1 + endind);
            ret_line = squeeze(RET(i,j, start_idx+1:end_idx));
            slopes = [];
            for tempi = 0:5:(iter_len-1)
                endcut = base_len + tempi;
                endcut = min(endcut, length(ret_line));
                y = squeeze(ret_line(1:endcut)) * lambda_um / (2 * pi * zSize_um);
                sy=size(y);
                if sy(2)==1;
                    y=y';
                end
                x = 0:(length(y)-1);

                y_mean = mean(y);
                x_mean = mean(x);
                numer = sum( squeeze(x - x_mean) .* squeeze(y - y_mean) );

                denom = sum((x - x_mean).^2);
                slope = abs(numer / denom);
                slopes(end+1) = slope;
            end
    
            if ~isempty(slopes)
                birefMap(i,j) = prctile(slopes, 95);  % Use 95th percentile
            else
                birefMap(i,j) = 0;
            end
        end
    end
end





function [Q,U,V,norms] = stokes_from_ori_ret(O3D_deg, R3D_deg, sigma)
%STOKES_FROM_ORI_RET Build smoothed, normalized Stokes volumes from ORI/RET (deg).
%   [Q,U,V,norms] = stokes_from_ori_ret(O3D_deg, R3D_deg, sigma)
%   sigma can be scalar or 3-element vector for imgaussfilt3.

if nargin < 3, sigma = 3; end
ORI = deg2rad(O3D_deg);         % radians
RET = deg2rad(R3D_deg);

Q = imgaussfilt3( sin(2*ORI) .* sin(2*RET), sigma );
U = imgaussfilt3( cos(2*ORI) .* sin(2*RET), sigma );
V = imgaussfilt3( cos(2*RET),               sigma );

norms = sqrt(Q.^2 + U.^2 + V.^2) + eps('single');
Q = Q ./ norms;
U = U ./ norms;
V = V ./ norms;
end

function [enface_full, enface_0_60, enface_0_100, enface_0_140, enface_0_180] = ...
    enface_orientation_stokes(O3D_deg, R3D_deg, Q, U, V, inten, surf, varargin)
%ENFACE_ORIENTATION_STOKES Unwrap ORI along z using Stokes cues, then circular-mean.
%   Outputs are [ny,nx] to mirror your saved NIfTIs (transpose before writing if needed).
%   Optional args:
%     'StartOffset' (default 10) samples below surf to start if room exists
%     'MaxLen'      (default 180) maximum depth samples considered
%     'CopyBackN'   (default 20)  number of initial samples forced to original ORI
%     'LowRetGate'  (default pi/8) threshold to bail to original ORI when local ret too small

p = inputParser;
addParameter(p,'StartOffset',10);
addParameter(p,'MaxLen',180);
addParameter(p,'CopyBackN',20);
addParameter(p,'LowRetGate',pi/8);
parse(p,varargin{:});
startOff = p.Results.StartOffset;
maxLen   = p.Results.MaxLen;
copyN    = p.Results.CopyBackN;
lowGate  = p.Results.LowRetGate;

[nx,ny,nz] = size(O3D_deg);
ORI = deg2rad(O3D_deg);
RET = deg2rad(R3D_deg);

% We'll accumulate per-(x,y) lines, unwrap, then mean
enface_full  = zeros(ny,nx,'single');
enface_0_60  = zeros(ny,nx,'single');
enface_0_100 = zeros(ny,nx,'single');
enface_0_140 = zeros(ny,nx,'single');
enface_0_180 = zeros(ny,nx,'single');

for y = 1:ny
    % local slices for speed
    ORI_xy = squeeze(ORI(:,y,:));   % [nx,nz]
    RET_xy = squeeze(RET(:,y,:));
    Q_xy   = squeeze(Q(:,y,:));
    U_xy   = squeeze(U(:,y,:));
    V_xy   = squeeze(V(:,y,:));
    inten_xy = squeeze(inten(:,y,:));

    % derive a local retardance from V for the gradient sign test
    localRET_fromV = acos( max(min(V_xy,1),-1) ); % clamp and acos

    for x = 1:nx
        lineI = inten_xy(x,:).';
        valid_len = sum(lineI > 0.01);
        s = max(0, round(surf(x,y))); % 0-based

        % choose start
        if valid_len - s > startOff
            starts_ = s + startOff;
        else
            starts_ = s;
        end
        if starts_ == 10, starts_ = 0; end
        starts_ = max(0, starts_);

        % compute end
        ends_ = min(nz-1, starts_ + min(valid_len - starts_, maxLen) - 1);

        if ends_ < starts_
            % nothing usable
            continue;
        end

        % segments (convert to 1-based for MATLAB indexing)
seg = (starts_+1):(ends_+1);

% Force everything to be column vectors
ori_seg   = squeeze(ORI_xy(x, seg));   ori_seg   = ori_seg(:);
q_seg     = squeeze(Q_xy(x,   seg));   q_seg     = q_seg(:);
u_seg     = squeeze(U_xy(x,   seg));   u_seg     = u_seg(:);
ret_local = squeeze(localRET_fromV(x, seg)); 
ret_local = ret_local(:);              % <-- column

% gradient of local retardance, pad one at end (vertical cat)
grad_ret = [diff(ret_local); 0];       % <-- semicolon, not comma


        % pick q-u or q+u based on dynamic range
        if (max(q_seg - u_seg) - min(q_seg - u_seg)) > (max(q_seg + u_seg) - min(q_seg + u_seg))
            pick = (q_seg - u_seg);
        else
            pick = (q_seg + u_seg);
        end

        % phase flip rule
        ori_adj = ori_seg;
        if mean(pick(grad_ret > 0), 'omitnan') > 0
            idx = (pick < 0);
            ori_adj(idx) = ori_adj(idx) + pi/2;
            ori_adj(ori_adj >  pi/2) = ori_adj(ori_adj >  pi/2) - pi;
        else
            idx = (pick > 0);
            ori_adj(idx) = ori_adj(idx) + pi/2;
            ori_adj(ori_adj >  pi/2) = ori_adj(ori_adj >  pi/2) - pi;
        end

        % copy back the first N samples from original ORI
        cbN = min(copyN, numel(ori_adj));
        ori_adj(1:cbN) = ori_seg(1:cbN);

        % bail out if local retardance is too small
        if mean(ret_local - prctile(ret_local,1), 'omitnan') < lowGate
            ori_adj = ori_seg;
        end

        % pack different windowed circular means
        enface_full(y,x)  = circ180_mean_deg(ori_adj);
        enface_0_60(y,x)  = circ180_mean_deg(ori_adj(1 : min(60,  numel(ori_adj))));
        enface_0_100(y,x) = circ180_mean_deg(ori_adj(1 : min(100, numel(ori_adj))));
        enface_0_140(y,x) = circ180_mean_deg(ori_adj(1 : min(140, numel(ori_adj))));
        enface_0_180(y,x) = circ180_mean_deg(ori_adj(1 : min(180, numel(ori_adj))));
    end
end
end

function mdeg = circ180_mean_deg(theta_rad)
% helper: circular mean with 180° periodicity (theta in radians)
if isempty(theta_rad)
    mdeg = single(0);
    return;
end
ang = 2*theta_rad(:);
c = cos(ang); s = sin(ang);
mean_angle = atan2(sum(s,'omitnan'), sum(c,'omitnan')); % [-pi,pi]
mdeg = single( (mean_angle/2) * (180/pi) );
end

function biref = biref_percentile_slope(V, RET_rad, surf, inten, zSize_um, lambda_mm, varargin)
%BIREFF_PERCENTILE_SLOPE Robust Δn via percentile of linear-fit slopes.
%   biref = biref_percentile_slope(V, RET_rad, surf, inten, zSize_um, lambda_um, ...
%                                  'PreferV',true, 'BaseLen',40, 'Step',5, ...
%                                  'MaxExtra',100, 'Percentile',80)
%
%   Units:
%     zSize_um     : axial pixel size [µm]
%     lambda_um    : wavelength [µm]
%   Δn is unitless (slope of OPD vs depth). OPD from retardance r (radians) is:
%     OPD = (r / (2*pi)) * lambda  .  Slope(OPD vs depth) ≈ Δn

p = inputParser;
addParameter(p,'PreferV',true);
addParameter(p,'BaseLen',40);
addParameter(p,'Step',5);
addParameter(p,'MaxExtra',100);
addParameter(p,'Percentile',80);
addParameter(p,'IntenThresh',0.01);
parse(p,varargin{:});
preferV   = p.Results.PreferV;
baseLen   = p.Results.BaseLen;
step      = p.Results.Step;
maxExtra  = p.Results.MaxExtra;
percP     = p.Results.Percentile;
thI       = p.Results.IntenThresh;

[nx,ny,nz] = size(V);
biref = zeros(nx,ny,'single');
depth_per_px_um = zSize_um;

for i = 1:nx
for j = 1:ny
    remaining = sum(inten(i,j,:) > thI) - surf(i,j);
    s0 = max(0, round(surf(i,j)));           % 0-based
    if remaining <= 0
        biref(i,j) = 0; continue;
    end

    if remaining > (baseLen + maxExtra)
        iter_len = maxExtra;
        endind   = baseLen + iter_len;
    elseif remaining > baseLen
        iter_len = remaining - baseLen;
        endind   = baseLen + iter_len;
    else
        iter_len = min(5, remaining);
        endind   = remaining;
    end

    % clamp indices to [1..nz]
    start_idx = min(nz, s0 + 1);              % convert to 1-based
    end_idx   = min(nz, s0 + endind);

    if end_idx <= start_idx, biref(i,j) = 0; continue; end

    Vseg = squeeze(V(i,j,start_idx:end_idx));
    retV = 0.5 * acos( max(min(Vseg,1),-1) );   % radians in [0,pi/2]

    % choose which retardance line to use
    retR = squeeze(RET_rad(i,j,start_idx:end_idx)); % radians
    if preferV
        % If V-based has *enough* modulation, use it; else fall back to RET
        useV = mean(retV - prctile(retV,1), 'omitnan') >= (pi/8);
        if (useV)

            ret_line = retV(:);
        else
            ret_line = retR(:);
        end
    else
        ret_line = retR(:);
    end

    % Construct a “built” cumulative retardance if needed (optional):
    % Here we simply use ret_line as-is. Keep the placeholders in case you
    % later want to reproduce your accumulative trick exactly.

    slopes = [];
    for add = 0:step:(iter_len-1)
        takeN = min(baseLen + add, numel(ret_line));
        y = (ret_line(1:takeN) / (2*pi)) * lambda_mm;        % OPD [µm]
        x = (0:(numel(y)-1))' * depth_per_px_um;             % depth [µm]
        % robust linear fit (ordinary LS here; replace with robustfit if you want)
        xmean = mean(x); ymean = mean(y);
        denom = sum( (x - xmean).^2 );
        if denom <= 0, continue; end
        slope = sum( (x - xmean).*(y - ymean) ) / denom;     % Δn (unitless)
        slopes(end+1) = abs(slope);
    end

    if ~isempty(slopes)
        biref(i,j) = single( prctile(slopes, percP) );
    else
        biref(i,j) = 0;
    end
end
end
end


