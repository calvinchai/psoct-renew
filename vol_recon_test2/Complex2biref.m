% new function to calculate biref


focus = niftiread('/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic9plus.nii');
% focus = niftiread('/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic10plus.nii');

%%% load data file
mosaicnum = 15;
% [192,193,194,195,196,229,228,227,226,225,234,235,236,237,238,271,270,269,268,267,279,278,277,276]
% [354,353,360,361,362,363,397,396,395,394,393,403,404,405,406,407,438,437,436,435,434,445,446,447,448,480,479,478,487,488]
for tilenum = [353]%1:1008
    %[354,353,360,361,362,363,397,396,395,394,393,403,404,405,406,407,438,437,436,435,434,445,446,447,448,480,479,478,487,488]
    %[313]%[309,310] %[312,313,314,317,318,319,320,275,276,277,278,279,280,281,266,267]
    disp(tilenum)
    % file_path_ = '/autofs/cluster/connects2/users/Nate/I80premotor_focustest_mineraloil_05092025/';
    file_path_ = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData/';
    filename = [file_path_ ,sprintf('mosaic_%03d_image_%04d_processed_cropped_focus.nii',mosaicnum, tilenum)];
    aa=niftiread(filename);
    s = size(aa,1)/4;
    Jones1 = aa(1:s,:,:) + 1j*aa(s+1:s*2,:,:);
    Jones2 = aa(s*2+1:s*3,:,:) + 1j*aa(s*3+1:s*4,:,:);

    dBI3D = flip(10*log10( abs(Jones1).^2+abs(Jones2).^2 ),3);

    R3D=flip(atan(abs(Jones1)./abs(Jones2))/pi*180,3);

    offset = 0;%100./180*pi;
    phase1 = (angle(Jones1));
    phase2 = (angle(Jones2));
    phi=(phase1-phase2)+offset*2;
    index1=phi>pi;
    phi(index1)=phi(index1)-2*pi;
    index2=phi<-pi;
    phi(index2)=phi(index2)+2*pi;
    O3D = flip( (phi)/2/pi*180 ,3);



    ORI = O3D / 180 * pi;
    RET = R3D / 180 * pi;


    sigma = [1,1,5]; %

    Q = imgaussfilt3( sin( ORI * 2 ) .* sin( RET * 2 ), sigma);
    U = imgaussfilt3( cos( ORI * 2 ) .* sin( RET * 2 ), sigma);
    V = imgaussfilt3( cos( RET * 2 ), sigma);

    norms = sqrt( Q.^2 + U.^2 + V.^2 );

    Q = Q ./ norms;
    U = U ./ norms;
    V = V ./ norms;
    tempR=acos(V)/2.;



%% calculate biref
    wavelength = .0013;

    surf_filename = [file_path_ ,sprintf('mosaic_%03d_image_%04d_processed_surface_finding.nii',mosaicnum, tilenum)]; 
    surf = niftiread(surf_filename); 
    surf = round((1024 - surf) - focus + 150); 
    surf = round( replace_outliers_iterative_matlab(surf, 21) );
    
    % z_px_size = Scan.Resolution(1,3);
    % biref = single(zeros(size(R3D,1),size(R3D,2)));
    % num_row = size(R3D,1);
    % num_col = size(R3D,2);
    % if(strcmpi(Scan.CropMethod,'surface'))
    %     parfor ii = 1:num_row
    %         for jj = 1:num_col
    %             ret_profile = squeeze(R3D(ii,jj,1:1 + depth_lin_ret))/360*wavelength;
    %             depth = z_px_size*(1:1 + depth_lin_ret);
    %             [p,~] = polyfit(depth,ret_profile',1);
    %             biref(ii,jj) = p(1); % Units: degree/um
    %         end
    %     end
    % 
    % 
    % else % CropMethod 'focus' or 'none'
    %     parfor ii = 1:num_row
    %     % for ii = 1:num_row
    %         for jj = 1:num_col
    % 
    %             basethresh = 0.;
    % 
    %             for depths_ = [ min(max(1,surf(ii,jj)+10),450):20:min(max(1,surf(ii,jj)+10)+200,450)]
    % 
    %                 ret_profile = squeeze(R3D(ii,jj,depths_:depths_+50)/(360)*wavelength);
    %                 depth = z_px_size*(depths_:depths_+50);
    % 
    %                 [p,~] = polyfit(depth,ret_profile',1);
    %                 if abs(p(1))> basethresh
    %                     basethresh = abs(p(1));
    %                     biref(ii,jj) = abs(p(1)); % Units: degree/um
    %                 end
    %             end
    %         end
    %     end
    % end 
    % 
    % %% save file
    % % fname = replace(sprintf('mosaic_%03i_image_%04i_processed_[modality].nii',mosaicnum,tilenum),'[modality]','biref');
    % % fpath = sprintf('%s/%s','/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData',fname);
    % % 
    % % fname_old = replace(sprintf('mosaic_%03i_image_%04i_processed_[modality]_old.nii',mosaicnum,tilenum),'[modality]','biref');
    % % fpath_old = sprintf('%s/%s','/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData',fname_old);
    % % 
    % % if ~isfile(fpath_old);
    % %     movefile(fpath,fpath_old);
    % %     niftiwrite(biref,fpath)
    % % else;
    % %     niftiwrite(biref,fpath)
    % % end
    % 
    % fname = replace(sprintf('mosaic_%03i_image_%04i_processed_[modality]_newmethod.nii',mosaicnum,tilenum),'[modality]','biref');
    % fpath = sprintf('%s/%s','/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData',fname);
    % niftiwrite(biref,fpath)

end


function corrected_image = replace_outliers_iterative_matlab(image, window_size)
% Replaces pixel values that are statistical outliers (outside mean +/- 3*std)
% in an image with the average of their local non-outlier neighbors.
%
% Args:
%     image (numeric matrix): The input grayscale image.
%     window_size (integer): The size of the square window (e.g., 3 for 3x3, 5 for 5x5)
%                            to calculate the local average. Must be an odd integer.
%                            Defaults to 5 if not provided.
%
% Returns:
%     double matrix: The image with outlier values corrected.

    % --- Argument Handling & Validation ---
    if nargin < 1
        error('Input image is required.');
    end
    if nargin < 2
        window_size = 5; % Default window_size if not provided
    end

    if ~isnumeric(image) || ~ismatrix(image) % Check if it's a 2D numeric array
        error('Input image must be a 2D numeric matrix.');
    end
    if ~isscalar(window_size) || ~isnumeric(window_size) || window_size < 1 || mod(window_size, 2) == 0 || floor(window_size) ~= window_size
        error('Window size must be a positive odd integer.');
    end

    % --- Calculate Global Statistics ---
    img_as_double = double(image); % Ensure double for calculations and for corrected_image
    
    if isempty(img_as_double)
        corrected_image = img_as_double; % Return empty if input is empty
        return;
    end

    mean_global = mean(img_as_double(:));
    std_global = std(img_as_double(:)); % MATLAB's default std normalizes by N-1

    % Define outlier thresholds
    lower_bound = mean_global - 2 * std_global ;
    upper_bound = mean_global + 2 * std_global;

    % --- Initialization ---
    corrected_image = img_as_double; % Work on a copy of the double-precision image
    [rows, cols] = size(img_as_double);

    % If standard deviation is zero, all pixels are the same (equal to mean_global).
    % No outliers are possible unless the image is empty (handled above).
    if std_global == 0 
        % No action needed, corrected_image is already a copy of the original.
        return;
    end

    % Find indices of outlier pixels in the original image
    % An outlier is a pixel p such that: p < lower_bound OR p > upper_bound
    [r_coords, c_coords] = find(img_as_double < lower_bound | img_as_double > upper_bound);

    pad_width = floor(window_size / 2);

    % --- Iterative Replacement ---
    for k = 1:length(r_coords)
        r = r_coords(k); % Current row of outlier pixel
        c = c_coords(k); % Current column of outlier pixel

        % Define window boundaries for the current pixel (r, c)
        % MATLAB is 1-indexed
        r_start = max(1, r - pad_width);
        r_end = min(rows, r + pad_width);
        c_start = max(1, c - pad_width);
        c_end = min(cols, c + pad_width);

        % Extract the neighborhood from the ORIGINAL image data (img_as_double)
        neighborhood = img_as_double(r_start:r_end, c_start:c_end);

        % Consider only non-outlier values (i.e., within global bounds) 
        % in this neighborhood for the average.
        valid_neighbors_mask = (neighborhood >= lower_bound) & (neighborhood <= upper_bound);
        valid_neighbors = neighborhood(valid_neighbors_mask);

        if ~isempty(valid_neighbors)
            % Calculate the mean of these valid, non-outlier neighbors
            mean_val = mean(valid_neighbors(:)); % Use (:) to ensure it's a vector mean
            corrected_image(r, c) = mean_val;
        else
            % Fallback: If all neighbors in the window are also outliers.
            % Replace with the global mean of the image. This is often a more
            % reasonable fallback than 0.0 when dealing with statistical outliers.
            corrected_image(r, c) = mean_global; 
            % Optional warning:
            % fprintf('Warning: No valid non-outlier neighbors for pixel (%d,%d). Replaced with global mean.\n', r, c);
        end
    end
end
