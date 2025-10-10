function stitch_section_modalities(basename, grids, varargin)
% stitch_section_modalities
% General Fiji grid-stitching launcher for multiple mosaics and modalities.
%
% Required:
%   basename : project folder (e.g. '/.../project')
%   grids    : [N x 2] numeric array of [grid_size_x grid_size_y]
%
% Optional name-value pairs:
%   'Mosaics'       : cellstr (default {'mosaic_001','mosaic_002'})
%   'Modalities'    : cellstr (default {'biref','ori','aip'})
%   'FijiPath'      : Fiji executable path
%   'ProcDir'       : processed directory name (default 'processed')
%   'MacroName'     : output macro file name (default 'stitch.ijm')
%   'OutputFormats' : cellstr of save formats (default {'Tiff','Jpeg'})
%   'IJ'            : struct of ImageJ stitching parameters
%
% Example:
%   grids = [10 8; 12 9];
%   stitch_section_modalities('/data/projectX', grids);
%
%   stitch_section_modalities('/data/projectX', grids, ...
%       'Mosaics', {'mosaic_001','mosaic_002','mosaic_003'}, ...
%       'Modalities', {'biref','ori','aip','ret'});

% -------------------
% Parse inputs
% -------------------
p = inputParser;
p.addRequired('basename', @(s)ischar(s) || isstring(s));
p.addRequired('grids', @(x)isnumeric(x) && size(x,2)==2);

p.addParameter('Mosaics', {'mosaic_001','mosaic_002'}, @(c)iscellstr(c) && ~isempty(c));
p.addParameter('Modalities', {'biref','ori','aip'}, @(c)iscellstr(c) && ~isempty(c));
p.addParameter('FijiPath', '/autofs/cluster/octdata2/users/Hui/Fiji.app/ImageJ-linux64', @(s)ischar(s) || isstring(s));
p.addParameter('ProcDir', 'processed', @(s)ischar(s) || isstring(s));
p.addParameter('MacroName', 'stitch.ijm', @(s)ischar(s) || isstring(s));
p.addParameter('OutputFormats', {'Tiff','Jpeg'}, @(c)iscellstr(c) && ~isempty(c));

IJdefaults = struct( ...
    'type', '[Grid: snake by columns]', ...
    'order', '[Up & Right]', ...
    'tile_overlap', 20, ...
    'first_file_index_i', 1, ...
    'output_textfile_name', 'Mosaic_depth001_slice001.txt', ...
    'fusion_method', '[Linear Blending]', ...
    'regression_threshold', 0.30, ...
    'max_avg_displacement_threshold', 2.50, ...
    'absolute_displacement_threshold', 3.50, ...
    'computation_parameters', '[Save memory (but be slower)]', ...
    'image_output', '[Fuse and display]' ...
);
p.addParameter('IJ', IJdefaults, @(s)isstruct(s));

p.parse(basename, grids, varargin{:});
basename       = char(p.Results.basename);
grids          = p.Results.grids;
mosaics        = p.Results.Mosaics;
modalities     = p.Results.Modalities;
fijiPath       = char(p.Results.FijiPath);
procDirName    = char(p.Results.ProcDir);
macroName      = char(p.Results.MacroName);
outputFormats  = p.Results.OutputFormats;
IJ             = mergeStruct(IJdefaults, p.Results.IJ);

% Ensure grids and mosaics match
if numel(mosaics) ~= size(grids,1)
    if size(grids,1) == 1 && numel(mosaics) > 1
        grids = repmat(grids, numel(mosaics), 1);
    elseif numel(mosaics) == 1 && size(grids,1) > 1
        mosaics = repmat(mosaics, 1, size(grids,1));
    else
        error('Number of mosaics (%d) must match number of grid rows (%d).', numel(mosaics), size(grids,1));
    end
end

% -------------------
% Prepare paths
% -------------------
procDir   = fullfile(basename, procDirName);
macroPath = fullfile(basename, macroName);
procDirIJ = strrep(procDir, filesep, '/');
baseIJ    = strrep(basename, filesep, '/');

% -------------------
% Build macro
% -------------------
macro = buildMacro(procDirIJ, baseIJ, mosaics, grids, modalities, outputFormats, IJ);

% -------------------
% Write macro
% -------------------
fid = fopen(macroPath, 'w');
assert(fid > 0, 'Cannot open %s for writing.', macroPath);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', macro);
clear cleanup;

% -------------------
% Run Fiji headless
% -------------------
cmd = sprintf('"%s" --headless --console -macro "%s"', fijiPath, macroPath);
status = system(cmd);
if status ~= 0
    error('Fiji/ImageJ returned non-zero status (%d). Command was:\n%s', status, cmd);
end
end


% ============================================================
% Build ImageJ macro
% ============================================================
function macro = buildMacro(procDirIJ, baseIJ, mosaics, grids, modalities, outputFormats, IJ)
lines = {};
for m = 1:numel(mosaics)
    mosaic = mosaics{m};
    gx = grids(m,1);
    gy = grids(m,2);
    for k = 1:numel(modalities)
        mod = modalities{k};
        pattern = sprintf('%s_image_{iiii}_processed_%s.nii', mosaic, mod);

        % Stitching command
        stitchCmd = sprintf(['run("Grid/Collection stitching", ' ...
            '"type=%s order=%s grid_size_x=%d grid_size_y=%d tile_overlap=%d ' ...
            'first_file_index_i=%d directory=[%s] file_names=%s ' ...
            'output_textfile_name=[%s] fusion_method=%s ' ...
            'regression_threshold=%.2f max/avg_displacement_threshold=%.2f ' ...
            'absolute_displacement_threshold=%.2f computation_parameters=%s ' ...
            'image_output=%s");'], ...
            IJ.type, IJ.order, gx, gy, IJ.tile_overlap, IJ.first_file_index_i, ...
            procDirIJ, pattern, IJ.output_textfile_name, IJ.fusion_method, ...
            IJ.regression_threshold, IJ.max_avg_displacement_threshold, ...
            IJ.absolute_displacement_threshold, IJ.computation_parameters, IJ.image_output);
        lines{end+1} = stitchCmd; %#ok<AGROW>

        % Save outputs
        outBase = sprintf('%s/%s_%s', baseIJ, mosaic, mod);
        for f = 1:numel(outputFormats)
            fmt = outputFormats{f};
            ext = lower(fmt);
            if strcmp(ext, 'tiff'), ext = 'tif'; end
            if strcmp(ext, 'jpeg'), ext = 'jpg'; end
            lines{end+1} = sprintf('saveAs("%s", "%s.%s");', fmt, outBase, ext);
        end
        lines{end+1} = 'close();';
    end
end
macro = strjoin(lines, sprintf('\n'));
macro = [macro sprintf('\n')];
end


% ============================================================
% Merge structs (b overwrites a)
% ============================================================
function c = mergeStruct(a, b)
c = a;
if isempty(b), return; end
fn = fieldnames(b);
for i = 1:numel(fn)
    c.(fn{i}) = b.(fn{i});
end
end
