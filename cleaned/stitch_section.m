
function stitch_section(basename, grid1, grid2, mosaic1, mosaic2, fijiPath)
% basename: project folder (e.g. '/.../project')
% grid1:    [grid_size_x grid_size_y] for mosaic1
% grid2:    [grid_size_x grid_size_y] for mosaic2
% mosaic1:  (optional) base name prefix, default 'mosaic_001'
% mosaic2:  (optional) base name prefix, default 'mosaic_002'
% fijiPath: (optional) Fiji executable path

if nargin < 4 || isempty(mosaic1)
    mosaic1 = 'mosaic_001';
end
if nargin < 5 || isempty(mosaic2)
    mosaic2 = 'mosaic_002';
end
if nargin < 6 || isempty(fijiPath)
    fijiPath = '/autofs/cluster/octdata2/users/Hui/Fiji.app/ImageJ-linux64';
end

procDir     = fullfile(basename,'processed');
macroPath   = fullfile(basename,'stitch.ijm');

procDirIJ   = strrep(procDir, filesep, '/');
baseIJ      = strrep(basename, filesep, '/');

% Output file names (Tiff + Jpeg)
m1_biref_tif = [baseIJ '/' 'mosaic_001_biref.tif'];
m1_biref_jpg = [baseIJ '/' 'mosaic_001_biref.jpg'];
m1_ori_tif   = [baseIJ '/' 'mosaic_001_ori.tif'];
m1_ori_jpg   = [baseIJ '/' 'mosaic_001_ori.jpg'];

m2_biref_tif = [baseIJ '/' 'mosaic_002_biref.tif'];
m2_biref_jpg = [baseIJ '/' 'mosaic_002_biref.jpg'];
m2_ori_tif   = [baseIJ '/' 'mosaic_002_ori.tif'];
m2_ori_jpg   = [baseIJ '/' 'mosaic_002_ori.jpg'];

% Build macro text
macro = sprintf([ ...
'run("Grid/Collection stitching", "type=[Grid: snake by columns] order=[Up & Right] grid_size_x=%d grid_size_y=%d tile_overlap=20 first_file_index_i=1 directory=[%s] file_names=%s_image_{iii}_processed_biref.nii output_textfile_name=[Mosaic_depth001_slice001.txt] fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n' ...
'saveAs("Tiff", "%s");\n' ...
'saveAs("Jpeg", "%s");\n' ...
'close();\n' ...
'run("Grid/Collection stitching", "type=[Grid: snake by columns] order=[Up & Right] grid_size_x=%d grid_size_y=%d tile_overlap=20 first_file_index_i=1 directory=[%s] file_names=%s_image_{iii}_processed_ori.nii    output_textfile_name=[Mosaic_depth001_slice001.txt] fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n' ...
'saveAs("Tiff", "%s");\n' ...
'saveAs("Jpeg", "%s");\n' ...
'close();\n' ...
'run("Grid/Collection stitching", "type=[Grid: snake by columns] order=[Up & Right] grid_size_x=%d grid_size_y=%d tile_overlap=20 first_file_index_i=1 directory=[%s] file_names=%s_image_{iii}_processed_biref.nii output_textfile_name=[Mosaic_depth001_slice001.txt] fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n' ...
'saveAs("Tiff", "%s");\n' ...
'saveAs("Jpeg", "%s");\n' ...
'close();\n' ...
'run("Grid/Collection stitching", "type=[Grid: snake by columns] order=[Up & Right] grid_size_x=%d grid_size_y=%d tile_overlap=20 first_file_index_i=1 directory=[%s] file_names=%s_image_{iii}_processed_ori.nii    output_textfile_name=[Mosaic_depth001_slice001.txt] fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n' ...
'saveAs("Tiff", "%s");\n' ...
'saveAs("Jpeg", "%s");\n' ...
'close();\n' ...
], ...
    grid1(1), grid1(2), procDirIJ, mosaic1, m1_biref_tif, m1_biref_jpg, ...
    grid1(1), grid1(2), procDirIJ, mosaic1, m1_ori_tif,   m1_ori_jpg, ...
    grid2(1), grid2(2), procDirIJ, mosaic2, m2_biref_tif, m2_biref_jpg, ...
    grid2(1), grid2(2), procDirIJ, mosaic2, m2_ori_tif,   m2_ori_jpg);

% Write macro
fid = fopen(macroPath,'w');
assert(fid>0, 'Cannot open %s for writing.', macroPath);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', macro);
clear cleanup;

% Run Fiji headless
cmd = sprintf('"%s" --headless --console -macro "%s"', fijiPath, macroPath);
status = system(cmd);
if status ~= 0
    error('Fiji/ImageJ returned non-zero status (%d). Command was:\n%s', status, cmd);
end
end