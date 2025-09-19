function WriteMacro_Serpentine( ParameterFile )

% WriteMacro_Serpentine( ParameterFile )
%   create Fiji Macro script that initiate serpentine stitching for
%   specified slices and modality
%
% USAGE:
%
%  INPUTS:
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters and WriteMacro:
%
%   ~2024-03-06~
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);

%%%% SETTING PARAMETERS
% origin          = WriteMacro.SpiralOrigin;
modality        = WriteMacro.Modality;
sliceid         = WriteMacro.SliceID;
mosaic_nums     = Parameters.MosaicID;
tile_nums       = Parameters.TileID;
fname_format    = replace(Enface.output_format,'[modality]',WriteMacro.Modality);
ind_format      = regexp(fname_format,'image_%0(\d+)i','tokens');
t_num_format    = regexp(fname_format,'image_%0\d+i','match');
t_num_length    = str2double(ind_format{1,1});
tile_format_str = ['image_{',repmat('i',1,t_num_length),'}'];
fiji_fname      = replace(fname_format,t_num_format,tile_format_str);
is_transposed   = Parameters.transpose;

%%%% SETTING DIRECTORIES
tiledir     = WriteMacro.TileDir;
AgarInfoDir = WriteMacro.AgarInfoDir;
fdir        = WriteMacro.FijiOutDir; if ~exist(fdir,'dir'); mkdir(fdir); end

fprintf(' - Tile directory = \n%s\n',tiledir);
fprintf(' - Fiji Input/Output coordinate directory = \n%s\n',fdir);
fprintf(' - Agar directory = \n%s\n',AgarInfoDir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS
fprintf(' - Loading Experiment file...\n %s\n', Parameters.ExpBasic);
S = whos('-file',Parameters.ExpBasic);
load(Parameters.ExpBasic);
Experiment  = eval(S(find(contains({S(:).name},'Experiment'))).name);

% Get index values of mosaic numbers from mosaic_nums
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    num_macro = 2; % 2 macros to write per slice
    steps = [Experiment.TilesPerSlice,Experiment.TilesPerSlice_tilt];
    idx = 0;
    idx_start = 1;
    i = 0;
    while idx_start <= length(mosaic_nums)
        i = i + 1;
        idx(i) = idx_start;
        idx_start = idx_start + steps(mod(i - 1,2) + 1);
    end
    mosaic_n = mosaic_nums(idx(1:2:end));
    mosaic_n_tilt = mosaic_nums(idx(2:2:end));
else
    num_macro = 1;
    idx = 1:Experiment.TilesPerSlice:length(mosaic_nums);
    mosaic_n = mosaic_nums(idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:num_macro
    switch n
        case 1 % Normal incidence
            fmacro = WriteMacro.SaveMacro;
            if exist(fmacro,'file'); delete(fmacro);
                fprintf(' --deleting existing Macro file-- \n%s\n',fmacro);end
            fid_Macro = fopen(fmacro, 'a'); %'a');
            mod_str = '';
            x_tile = Experiment.X_tile;
            y_tile = Experiment.Y_tile;
            z_tile = Experiment.Z_tile;
        case 2 % Tilted incidence
            fmacro = WriteMacro.SaveMacro_tilt;
            if exist(fmacro,'file'); delete(fmacro);
                fprintf(' --deleting existing Macro file-- \n%s\n',fmacro);end
            fid_Macro = fopen(fmacro, 'a'); %'a');
            mod_str = '_tilt';
            x_tile = Experiment.X_tile_tilt;
            y_tile = Experiment.Y_tile_tilt;
            z_tile = Experiment.Z_tile_tilt;
    end
    
    for ssind=1:length(sliceid)
        curr_slice = sliceid(ssind);
        
        % currstarttile = (currslice-1)*Experiment.TilesPerSlice + Experiment.First_Tile;
        curr_start_tile = 1; % Every mosaic restarts with tile 1
        slice_num = sprintf('%03i',curr_slice);
        
        switch n
            case 1 % Normal incidence 
                fname = sprintf(fiji_fname,mosaic_n(ssind));
            case 2 % Tilted incidence
                fname = sprintf(fiji_fname,mosaic_n_tilt(ssind));
        end

        for zz = 1:z_tile
            fmosaic = ['Mosaic',mod_str,'_depth' sprintf('%03i',zz) '_slice' slice_num '.txt']; % relative path to tiledir
            % fmosaic = [fdir,'/Mosaic',mod_str,'_depth' sprintf('%03i',zz) '_slice' slice_num '.txt']; % relative path to tiledir
            
            % Create the command line for the actual stitching
            if(~is_transposed) %%%%%%%%%%% If tiles have not been transposed %%%%%%%%%%%
                m1=sprintf(['run("Grid/Collection stitching", ' ...
                    '"type=[Grid: snake by columns] order=[Up & Right] ' ...
                    'grid_size_x=%i grid_size_y=%i tile_overlap=%i ' ...
                    'first_file_index_i=%i directory=[%s] ' ...
                    'file_names=%s output_textfile_name=[%s] ' ...
                    'fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ' ...
                    'compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n'],...
                    x_tile, y_tile, Experiment.PercentOverlap, curr_start_tile, tiledir, fname, fmosaic);
            elseif(is_transposed) %%%%%%%%%%% If tiles have been transposed %%%%%%%%%%%
                m1=sprintf(['run("Grid/Collection stitching", ' ...
                    '"type=[Grid: snake by rows] order=[Down & Left] ' ...
                    'grid_size_x=%i grid_size_y=%i tile_overlap=%i ' ...
                    'first_file_index_i=%i directory=[%s] ' ...
                    'file_names=%s output_textfile_name=[%s] ' ...
                    'fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 ' ...
                    'compute_overlap computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");\n'],...
                    y_tile, x_tile, Experiment.PercentOverlap, curr_start_tile, tiledir, fname, fmosaic);
            end

            m2=['saveAs("Tiff", "' fdir '/' modality mod_str '_Slice' slice_num '.tif");\n'];
            m3 = 'close();\n';
            fprintf(fid_Macro,m1);
            fprintf(fid_Macro,m2);
            fprintf(fid_Macro,m3);
        end
    end
    fclose(fid_Macro);
end




end
%% spiral / roll-by-roll / multiple small patches (get step & offset) / y strip first then x direction

