nii_folder = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/nii2D_transpose';

number_folder = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/StitchingFiji/Stitching_diagram/';

output_folder = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/OverlayOutput';
 
% Create output directory if it doesn't exist

if ~exist(output_folder, 'dir')

    mkdir(output_folder);

end
 
% Loop through all 588 images

for idx = 1:588

    % Format filenames

    nii_filename = sprintf('mosaic_015_image_%04d_processed_mip.nii', idx);

    nii_path = fullfile(nii_folder, nii_filename);

    number_filename = sprintf('number_img_%03d.tiff', idx);

    number_path = fullfile(number_folder, number_filename);
 
    % Skip if files don't exist

    if ~isfile(nii_path)

        warning('Missing NII file: %s', nii_filename);

        continue;

    end

    if ~isfile(number_path)

        warning('Missing number image: %s', number_filename);

        continue;

    end
 
    % Load .nii image (make sure you have the NIfTI toolbox)

    nii_img = niftiread(nii_path);

    % nii_img = squeeze(nii.img);  % Remove singleton dimensions
 
    % Normalize and convert to 8-bit grayscale

    nii_img = double(nii_img)';

    nii_img = uint8(255 * mat2gray(nii_img));
 
    % Read number image and resize

    number_img = imread(number_path);

    number_img = imresize(number_img, size(nii_img));  % Ensure sizes match
 
    % Convert both to RGB if needed

    if size(nii_img,3) == 1

        nii_rgb = repmat(nii_img, 1, 1, 3);

    else

        nii_rgb = nii_img;

    end

    if size(number_img,3) == 1

        number_rgb = repmat(number_img, 1, 1, 3);

    else

        number_rgb = number_img;

    end
 
    % Overlay with alpha blending

    alpha = 0.3;  % Opacity of number overlay

    overlay = uint8((1 - alpha) * double(nii_rgb) + alpha * double(number_rgb));
 
    % Save to .tiff

    output_name = sprintf('mosaic_015_image_%04d_overlay.tiff', idx);

    output_path = fullfile(output_folder, output_name);

    imwrite(overlay, output_path);
 
    fprintf('Saved: %s\n', output_path);

end


fmacro = WriteMacro.SaveMacro;
if exist(fmacro,'file'); delete(fmacro);
    fprintf(' --deleting existing Macro file-- \n%s\n',fmacro);end
fid_Macro = fopen(fmacro, 'a'); %'a');


m1=sprintf(['run("Grid/Collection stitching", ' ...
    '"type=[Grid: snake by rows] order=[Down & Left] ' ...
    'grid_size_x=%i grid_size_y=%i tile_overlap=%i ' ...
    'first_file_index_i=%i directory=[%s] ' ...
    'file_names=%s output_textfile_name=[%s] ' ...
    'fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ' ...
    'compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n'],...
    y_tile, x_tile, Experiment.PercentOverlap, curr_start_tile, tiledir, fname, fmosaic);


m2=['saveAs("Tiff", "' fdir '/' modality mod_str '_Slice' slice_num '.tif");\n'];
m3 = 'close();\n';

fprintf(fid_Macro,m1);
fprintf(fid_Macro,m2);
fprintf(fid_Macro,m3);

