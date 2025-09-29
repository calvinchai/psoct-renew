% Set mosaic number
function recalculate_biref_orientation_I80(mosaicnum)
mosaicnum=str2num(mosaicnum);
% Define file path
file_path = '/vast/fiber/projects/20250920_CCtest4_30degrees_CW/';
outdir = '/autofs/space/megaera_001/users/kchai/project/20250920_CCtest4_30degrees_CW_UpdatedBiref/processed/';
if mod(mosaicnum,2)==1
    all_tissues = [1:13*5];
else
    all_tissues = [1:9*5];
end

% disp('start parpool');
% poolobj = parpool(12);
% disp('parpool started');
disp(mosaicnum/2)
% parpool('threads');

% Loop over all tiles (you should define `all_tissues` in MATLAB beforehand, e.g., all_tissues = 1:588;)
parfor tilenum = all_tissues
    tilenum
    % Construct full filename
    filename = [file_path, sprintf(['mosaic_%03d_image_%03d_processed_cropped.nii'], mosaicnum, tilenum)];
    % disp(file_path)
    % disp(sprintf(['mosaic_%03d_image_%04d_processed_cropped_focus.nii'], mosaicnum, tilenum))
    % Load .nii image
    aa = double(niftiread(filename));  % Convert to double for computation

    % Calculate dimension split
    s = floor(size(aa,1) / 4);

    % Construct complex Jones matrices
    Jones1 = aa(1:s,:,:) + 1i * aa(s+1:2*s,:,:);
    Jones2 = aa(2*s+1:3*s,:,:) + 1i * aa(3*s+1:4*s,:,:);

    % Calculate R3D (retardance angle in degrees), flip along Z
    R3D = flip(atan(abs(Jones1) ./ abs(Jones2)) / pi * 180, 3);

    % Phase calculations
    offset = 0;  % in radians
    phase1 = angle(Jones1);
    phase2 = angle(Jones2);
    phi = (phase1 - phase2) + offset * 2;

    % Wrap phase to [-π, π]
    phi(phi > pi) = phi(phi > pi) - 2*pi;
    phi(phi < -pi) = phi(phi < -pi) + 2*pi;

    % Final orientation map
    O3D = flip(phi / (2 * pi) * 180, 3);

    R3D=shift_non_nan_up(R3D);
    O3D=shift_non_nan_up_fill(O3D);

    % Orientation and Retardance in radians
    ORI = O3D / 180 * pi;
    RET = R3D / 180 * pi;
    
    % Compute intensity
    inten = 10 * log10( sqrt( abs(Jones1).^2 + abs(Jones2).^2 ) );
    
    % Flip intensity along Z-axis (3rd dimension)
    inten = inten(:,:,end:-1:1);
    
    % Shift non-Inf values upward
    inten = shift_non_inf_up(inten);


    
    % Compute Q, U, V directly without allocating ORI2 or RET2
    Q = imgaussfilt3(sin(2 * ORI) .* sin(2 * RET), 3);
    U = imgaussfilt3(cos(2 * ORI) .* sin(2 * RET), 3);
    V = imgaussfilt3(cos(2 * RET), 3);
    
    % Normalize the Stokes vectors
    norms = sqrt(Q.^2 + U.^2 + V.^2) + eps;  % eps to avoid divide-by-zero
    Q = Q ./ norms;
    U = U ./ norms;
    V = V ./ norms;
    

    % Assume inten, ORI, and RET are all defined and 3D volumes (nx, ny, nz)

    [nx, ny, nz] = size(inten);
    surf = zeros(nx, ny);  % Output surface map
    
    % w = 5;
    % w2 = 5;
    % kernel = [-ones(1, w)/w, ones(1, w2)/w2];  % Gradient kernel
    % 
    % % Loop over each (x,y) A-line
    % for i = 1:nx
    %     for j = 1:ny
    %         line = squeeze(inten(i, j, :));
    %         valid_len = sum(line > 0.01);
    %         if valid_len > w + w2
    %             data = imgaussfilt(line(1:valid_len), 5);  % 1D Gaussian smoothing
    %             grad = -conv(data, kernel, 'valid');
    %             positions = (w+1):(valid_len-w+1);
    %             [grad_min, idx_min] = max(grad);
    %             i_min = positions(idx_min);
    %             surf(i,j) = i_min;
    %         else
    %             surf(i,j) = 0;
    %         end
    %     end
    % end
    % 
    % surf = zeros(nx, ny);  % Initialize surface map



    % 
    % % Compute Q, U, V directly without allocating ORI2 or RET2
    % Q = imgaussfilt3(sin(2 * ORI) .* sin(2 * RET), 3);
    % U = imgaussfilt3(cos(2 * ORI) .* sin(2 * RET), 3);
    % V = imgaussfilt3(cos(2 * RET), 3);
    % 
    % % Normalize the Stokes vectors
    % norms = sqrt(Q.^2 + U.^2 + V.^2) + eps;  % eps to avoid divide-by-zero
    % Q = Q ./ norms;
    % U = U ./ norms;
    % V = V ./ norms;

    V3=V(:,:,:);

    w = 5;
    w2 = 5;
    kernel = [-ones(1, w)/w, ones(1, w2)/w2];  % Gradient kernel
    
    for i = 1:nx
        for j = 1:ny
            % line = squeeze(inten(i, j, :));
            % valid_len = sum(line > 0.01);
            % norm_profile = squeeze(norms(i, j, 1:min(10, end)));  % match [:10] behavior
            % 
            % if mean(norm_profile(:)) > 0.5
            %     surf(i, j) = 0;  % override surface
            % elseif valid_len > w + w2
            % 
            % 
            %     data = imgaussfilt(line(1:valid_len), 5);  % 1D Gaussian smoothing
            %     grad = -conv(data, kernel, 'valid');
            %     positions = (w+1):(valid_len-w+1);
            %     [grad_min, idx_min] = max(grad);
            %     i_min = positions(idx_min);
            %     surf(i,j) = i_min;
            % else
            surf(i, j) = 0;
            % end
        end
    end


    % Median filtering of surface map
    surf = medfilt2(surf, [3, 3],'symmetric');


    % volume_ori: cell array of depth profiles
    % volume_ori = cell(size(ORI, 2), size(ORI, 1));  % [Y, X]
    % [nx, ny, nz] = size(ORI);
    % for lays_a = 1:ny  % Y axis
    % inputret = squeeze(RET(:, lays_a, :));  % [X, Z]
    % inputori = squeeze(ORI(:, lays_a, :));
    % for lays = 1:nx  % X axis
    % inten_line = squeeze(inten(lays, lays_a, :));
    % valid_len = sum(inten_line > 0.01);
    % s = surf(lays, lays_a);  % surface index
    % if valid_len - s > 10
    % starts_ = s + 10;
    % else
    % starts_ = s;
    % end
    % ends_ = starts_ + min(valid_len - starts_, 40);
    % starts_ = max(1, round(starts_));
    % ends_   = min(nz, round(ends_));
    % measured_ret_noisy = inputret(lays, starts_:ends_);
    % measured_ori_noisy = inputori(lays, starts_:ends_);
    % volume_ori{lays_a, lays} = measured_ori_noisy;  % [Y, X]
    % end
    % end

    volume_ori = cell(ny, nx);       % [Y, X]
    volume_ori_ref = cell(ny, nx);   % [Y, X]
    
    for lays_a = 1:ny  % Y-axis
        inputret = acos(squeeze(V(:, lays_a, :)));                     % [X, Z]
        inputori = atan2(squeeze(Q(:, lays_a, :)), squeeze(U(:, lays_a, :)));  % [X, Z]
    
        for lays = 1:nx  % X-axis
            inten_line = squeeze(inten(lays, lays_a, :));
            valid_len = sum(inten_line > 0.01);
            s = surf(lays, lays_a);  % surface index
    
            if valid_len - s > 10
                starts_ = round(s + 10);
            else
                starts_ = round(s);
            end
    
            if starts_ == 10
                starts_ = 0;
            end
    
            % Boundaries
            starts_ = max(1, starts_);
            ends_ = starts_ + min(valid_len - starts_, 120);
            ends_ = min(nz, ends_);
    
            % Extract slices
            actual_ori = squeeze(ORI(lays, lays_a, starts_:ends_));
            ref_ori = actual_ori(:);  % keep for reference
    
            q_seg = squeeze(Q(lays, lays_a, starts_:ends_));
            u_seg = squeeze(U(lays, lays_a, starts_:ends_));
    
            % Pick + or - logic
            if (max(q_seg - u_seg) - min(q_seg - u_seg)) > (max(q_seg + u_seg) - min(q_seg + u_seg))
                pick_ = q_seg - u_seg;
            else
                pick_ = q_seg + u_seg;
            end
    
            localret = squeeze(inputret(lays, starts_:ends_));
            % gradients_ret = [diff(localret); 0];  % pad one at end

            localret = localret(:);  % ensure column
            gradients_ret = [diff(localret); 0];  % pad at the end
    
            % Phase flip logic
            if mean(pick_(gradients_ret > 0), 'omitnan') > 0
                idx = pick_ < 0;
                actual_ori(idx) = actual_ori(idx) + pi/2;
                actual_ori(actual_ori > pi/2) = actual_ori(actual_ori > pi/2) - pi;
            else
                idx = pick_ > 0;
                actual_ori(idx) = actual_ori(idx) + pi/2;
                actual_ori(actual_ori > pi/2) = actual_ori(actual_ori > pi/2) - pi;
            end
    
            % Enforce first 20 values from ORI
            copy_len = min(20, ends_ - starts_ + 1);
            actual_ori(1:copy_len) = squeeze(ORI(lays, lays_a, starts_:starts_ + copy_len - 1));
    
            % Fallback if local retardance is low
            if mean(localret(:) - prctile(localret(:),1)  )/2 < pi/6
                actual_ori = squeeze(ORI(lays, lays_a, starts_:ends_));
            end
    
            volume_ori{lays_a, lays} = actual_ori;
            volume_ori_ref{lays_a, lays} = ref_ori;
        end
    end

    % Compute enface image by circular mean of ORI along Z
    enface = zeros(ny, nx);  % [Y, X] for visualization
    enface_ref = zeros(ny, nx);

    for i = 1:ny
    for j = 1:nx
    ori_line = volume_ori{i, j};
    ori_line_ref = volume_ori_ref{i, j};

    if isempty(ori_line)
    enface(i, j) = 0;
    else
    % MATLAB circular mean over [-pi/2, pi/2]
    angles = 2*ori_line(:);
    x = cos(angles);
    y = sin(angles);
    mean_angle = atan2(sum(y), sum(x));  % result in [-pi, pi]
    % Clamp to [-pi/2, pi/2]
    enface(i, j) = mean_angle / 2 / pi * 180;
    end

    if isempty(ori_line_ref)
    enface_ref(i, j) = 0;
    else
    % MATLAB circular mean over [-pi/2, pi/2]
    angles = 2*ori_line_ref(:);
    x = cos(angles);
    y = sin(angles);
    mean_angle = atan2(sum(y), sum(x));  % result in [-pi, pi]
    % Clamp to [-pi/2, pi/2]
    enface_ref(i, j) = mean_angle / 2 / pi * 180;
    end

    end
    end


    % % Constants
    % lambda = 0.013;  % Wavelength in mm
    % depth_per_pixel = 2.5;  % Axial pixel spacing in µm
    % [nx, ny, nz] = size(V);
    % biref = zeros(nx, ny);
    % for i = 1:nx
    % for j = 1:ny
    % remaining = sum(inten(i,j,:) > 0.01) - surf(i,j);
    % if remaining > 140
    % base_len = 40;
    % iter_len = 100;
    % elseif remaining > 40
    % base_len = 40;
    % iter_len = remaining - 40;
    % elseif remaining > 0
    % base_len = remaining;
    % iter_len = 5;
    % else
    % biref(i,j) = 0;
    % continue;
    % end
    % slopes = [];
    % for tempi = 0:5:(iter_len-1)
    % z_start = round(surf(i,j));
    % z_end = min(nz, z_start + base_len + tempi);
    % z_range = z_start:z_end;
    % % Extract V segment
    % vseg = V(i, j, z_range);
    % vseg = squeeze(vseg);
    % if isempty(vseg)
    % continue;
    % end
    % % Compute phase retardation
    % y = acos(vseg) / 2 * lambda / (2 * pi * depth_per_pixel);
    % x = 0:(length(y)-1);
    % % Linear regression
    % y_mean = mean(y);
    % x_mean = mean(x);
    % numer = sum((x - x_mean) .* (y' - y_mean));
    % denom = sum((x - x_mean).^2);
    % slope = abs(numer / denom);
    % slopes(end+1) = slope;
    % end
    % if ~isempty(slopes)
    % biref(i,j) = prctile(slopes, 95);  % 95th percentile
    % else
    % biref(i,j) = 0;
    % end
    % end
    % end



    % Compute Q, U, V directly without allocating ORI2 or RET2
    Q = imgaussfilt3(sin(2 * ORI) .* sin(2 * RET), [3 3 10]);
    U = imgaussfilt3(cos(2 * ORI) .* sin(2 * RET), [3 3 10]);
    V = imgaussfilt3(cos(2 * RET), [3 3 10]);
    
    % Normalize the Stokes vectors
    norms = sqrt(Q.^2 + U.^2 + V.^2) + eps;  % eps to avoid divide-by-zero
    Q = Q ./ norms;
    U = U ./ norms;
    V = V ./ norms;
    
    % Constants
    lambda = 0.013;                % Wavelength in mm
    depth_per_pixel = 2.5;  % Convert µm to mm
    
    [nx, ny, nz] = size(V);
    biref = zeros(nx, ny);
    biref_v = zeros(nx, ny);
    
    for i = 1:nx
        for j = 1:ny
            remaining = sum(inten(i,j,:) > 0.01) - surf(i,j);
            surf_ij = round(surf(i,j));
    
            if remaining > 100
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
                biref_v(i,j) = 0;
                continue;
            end
    
            % Ensure we don’t go out of bounds
            start_idx = max(1, surf_ij);
            end_idx = min(nz, surf_ij + endind);
    
            % Compute V-based retardance
            Vseg = squeeze(V(i, j, start_idx:end_idx));
            ret_raw = acos(Vseg) / 2;
    
            % Use RET directly if mean retardance < π/5
            if mean(ret_raw-prctile(ret_raw,1), 'omitnan') < pi/8
                ret_line = squeeze(RET(i,j, start_idx+1:end_idx));
            else
                Vseg1 = squeeze(V(i,j, start_idx+1:end_idx));
                Vseg0 = squeeze(V(i,j, start_idx:(end_idx-1)));
                diffs = abs((acos(Vseg1) - acos(Vseg0)) / 2);
                rstret = acos(V(i,j,start_idx)) / 2;
                ret_build = zeros(size(diffs)+[1,0]);
                ret_build(1) = rstret;
    
                for k = 1:length(diffs)
                    ret_build(k+1) = ret_build(k) + diffs(k);
                end
    
                ret_line = reshape( ret_build(:) ,[],1) + ...
                    reshape( squeeze(RET(i,j,start_idx:end_idx))',[],1 ) - reshape( (acos(Vseg)/2),[],1);
            end
    
            slopes = [];
            for tempi = 0:5:(iter_len-1)
                endcut = base_len + tempi;
                endcut = min(endcut, length(ret_line));
                y = squeeze(ret_line(1:endcut)) * lambda / (2 * pi * depth_per_pixel);
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
                biref(i,j) = prctile(slopes, 95);  % Use 95th percentile
            else
                biref(i,j) = 0;
            end


            slopes = [];
            for tempi = 0:5:(iter_len-1)
                z_start = max(1,round(surf(i,j)));
                z_end = min(nz, z_start + base_len + tempi);
                z_range = z_start:z_end;
                % Extract V segment
                vseg = V3(i, j, z_range);
                vseg = squeeze(vseg);
                if isempty(vseg)
                continue;
                end
                % Compute phase retardation
                y = acos(vseg) / 2 * lambda / (2 * pi * depth_per_pixel);
                x = 0:(length(y)-1);
                % Linear regression
                y_mean = mean(y);
                x_mean = mean(x);
                numer = sum((x - x_mean) .* (y' - y_mean));
                denom = sum((x - x_mean).^2);
                slope = abs(numer / denom);
                slopes(end+1) = slope;
            end
            if ~isempty(slopes)
            biref_v(i,j) = prctile(slopes, 95);  % 95th percentile
            else
            biref_v(i,j) = 0;
            end




        end
    end



    % Define file paths
    real_name = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_reprocessed_orientation.nii', mosaicnum, tilenum));
    real_name_ref = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_reprocessed_orientation_shallow.nii', mosaicnum, tilenum));
    replace_name = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_processed_orientation.nii', mosaicnum, tilenum));

    % If backup doesn't exist, rename the original

    % if ~isfile(replace_name)
    %     movefile(real_name, replace_name);  % rename
    % end

    % Save the updated orientation map as a NIfTI (transposed to match original orientation)
    enface_transposed = single(enface');  % must be numeric for NIfTI
    niftiwrite(enface_transposed, real_name, 'Compressed', false);  % save new NIfTI

    % Save the reference orientation map as a NIfTI (transposed to match original orientation)
    enface_transposed = single(enface_ref');  % must be numeric for NIfTI
    niftiwrite(enface_transposed, real_name_ref, 'Compressed', false);  % save new NIfTI


    % Load the original orientation NIfTI image (for comparison)

    if isfile(replace_name)
        oriimg = niftiread(replace_name);
    else
        oriimg = zeros(size(enface_transposed));
    end

    oriimg_transposed = double(oriimg);  % transpose to match enface display
    % % Plot side-by-side: before and after unwarping
    % fig = figure('Visible', 'off');  % suppress GUI
    % tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    % % Before
    % nexttile;
    % imagesc(oriimg_transposed);
    % axis image off;
    % colormap(hsv);
    % caxis([-90 90]);
    % title('before');
    % % After
    % nexttile;
    % imagesc(enface_transposed);
    % axis image off;
    % colormap(hsv);
    % caxis([-90 90]);
    % title('after');
    % outdir = [outdir num2str(mosaicnum) '/'];
    % if ~exist(outdir, 'dir')
    % mkdir(outdir);
    % end
    % saveas(fig, fullfile(outdir, sprintf('%d.png', tilenum)));
    % close(fig);




    % Define file paths
    real_name = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_reprocessed_biref.nii', mosaicnum, tilenum));
    real_name_v = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_reprocessed_biref_v.nii', mosaicnum, tilenum));

    replace_name = fullfile(outdir, sprintf('mosaic_%03d_image_%04d_processed_biref.nii', mosaicnum, tilenum));
    % % Backup the old file if not already backed up
    % if ~isfile(replace_name)
    % movefile(real_name, replace_name);  % rename original
    % end
    % Save the new biref map
    biref_single = single(biref);  % ensure single precision for NIfTI
    niftiwrite(biref_single, real_name, 'Compressed', false);
    biref_single_v = single(biref_v);  % ensure single precision for NIfTI
    niftiwrite(biref_single_v, real_name_v, 'Compressed', false);
    % Load the original biref map
    if isfile(replace_name)
        oriimg = niftiread(replace_name);
        img_old = double(oriimg);  % transpose to match

        vals_old = img_old(:);
        low1 = prctile(vals_old, 1);
        high1 = prctile(vals_old, 99);
    else
        oriimg = zeros(size(biref));
        img_old = double(oriimg);
        low1 = 0;
        high1 = 1;
    end
    % Compute display limits using 1st and 99th percentiles

    vals_new = biref(:);
    low2 = prctile(vals_new, 1);
    high2 = prctile(vals_new, 99);
    % % Create figure
    % fig = figure('Visible', 'off');
    % tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    % % Before
    % nexttile;
    % imagesc(img_old);
    % colormap gray;
    % axis image off;
    % caxis([low1, high1]);
    % title('before');
    % % After
    % nexttile;
    % imagesc(biref);
    % colormap gray;
    % axis image off;
    % caxis([low2, high2]);
    % title('after');
    % 
    % 
    % % outdir = [outdir num2str(mosaicnum) '/'];
    % if ~exist(outdir, 'dir')
    % mkdir(outdir);
    % end
    % saveas(fig, fullfile(outdir, sprintf('%d_biref.png', tilenum)));
    % close(fig);





end
end






function out = shift_non_nan_up_fill(inten)
% For each (x,y) A-line along the 3rd axis, remove zeros and shift up,
% then pad the remaining with the last valid value.

[nx, ny, nz] = size(inten);
out = zeros(nx, ny, nz, 'like', inten);

for i = 1:nx
    for j = 1:ny
        line = squeeze(inten(i,j,:));
        valid = line(line ~= 0);
        if isempty(valid)
            out(i,j,:) = NaN;
        else
            new_line = zeros(nz, 1, 'like', inten);
            new_line(1:length(valid)) = valid;
            new_line(length(valid)+1:end) = valid(end);
            out(i,j,:) = new_line;
        end
    end
end
end



function out = shift_non_nan_up(inten, pad_with)
% Shift non-NaN values up along 3rd axis, pad with specified value.

% if nargin < 2
%     pad_with = 0;
% end

pad_with = 0;
[nx, ny, nz] = size(inten);
out = pad_with * ones(nx, ny, nz, 'like', inten);

for i = 1:nx
    for j = 1:ny
        line = squeeze(inten(i,j,:));
        valid = line(~isnan(line));
        if ~isempty(valid)
            out(i,j,1:length(valid)) = valid;
        end
    end
end
end


function out = shift_non_inf_up(inten)
% shift_non_inf_up shifts non-Inf values upward along the 3rd dimension
% and pads with the specified value.
%
% Parameters:
%   inten     : 3D array (nx, ny, nz)
%   pad_with  : scalar, value to pad with (default = 0)
%
% Returns:
%   out       : 3D array, same size, values compacted along axis-3

% if nargin < 2
%     pad_with = 0;
% end
pad_with = 0;
[nx, ny, nz] = size(inten);
flat = reshape(inten, [], nz);  % (nx*ny, nz)

% Initialize output
if isnan(pad_with)
    out_flat = nan(size(flat), 'like', inten);
else
    out_flat = pad_with * ones(size(flat), 'like', inten);
end

% Loop through each line and shift non-Inf values up
for idx = 1:size(flat, 1)
    line = flat(idx, :);
    valid = line(~isinf(line));
    if ~isempty(valid)
        out_flat(idx, 1:numel(valid)) = valid;
        % rest remains pad_with
    end
end

% Reshape back to 3D
out = reshape(out_flat, nx, ny, nz);
end