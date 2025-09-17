function [ Experiment_Fiji ] = ReadFijiTileCoordinates_NB( ParameterFile, saveflag )
%
% [ Experiment_Fiji ] = ReadFijiTileCoordinates( ParameterFile, saveflag )
%   Improve coordinates from registered coordinates files of Fiji computed overlap
%
% USAGE:
%
%  INPUTS:
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters and ReadCoord:
%       saveflag          =   save Experiment_Fiji if true
%
%   ~2024-03-06~
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);

%%%% DISPLAY METHOD
Method      = ReadCoord.Method;

if      Method==1; fprintf(' -- Saving median coordinates -- \n');
elseif  Method==2; fprintf(' -- Saving reconstructed coordinates from median stepsize and offset -- \n');
end

%%%% INCLUDE SLICE
sliceidx     = ReadCoord.ReadSliceID;
nslices{1,1} = length(sliceidx);
%%%% SAVE EXPERIMENT
fout        = ReadCoord.SaveExp;
fprintf(' - Save Experiment path = %s\n',fout);
%%%% LOAD TEMPLATE
Experiment  = ReadCoord.ExperimentBasic;

%%%% SET COORDINATE & TILE BASENAME
if FijiInfo.isRegistered==1
    fiji_file_ext = '.registered.txt';else;fiji_file_ext = '.txt';
end

pathname_coord{1,1} = [FijiInfo.OutDir '/' FijiInfo.FileBase '%03i' fiji_file_ext];
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    num_mosaics_per_slice = 2;
    sliceidx_tilt = ReadCoord.ReadSliceID_tilt;
    nslices{2,1} = length(sliceidx_tilt);
    pathname_coord{2,1} = [FijiInfo.OutDir '/' FijiInfo.FileBase_tilt '%03i' fiji_file_ext];
    mosaic_nums{1,1} = 2*sliceidx - 1;
    mosaic_nums{2,1} = 2*sliceidx_tilt;
else
    num_mosaics_per_slice = 1;
    mosaic_nums{1,1} = sliceidx;
end


% pathname_tile  = [WriteMacro.TileDir filesep WriteMacro.Modality '_'];
fname_format    = replace(Enface.output_format,'[modality]',WriteMacro.Modality);
m_num_format    = regexp(fname_format,'mosaic_%0\d+i','match');
m_num_length    = length(m_num_format{1,1});
split_idx       = strfind(fname_format,m_num_format) + m_num_length;
fname_format1   = fname_format(1:split_idx);
fname_format2   = fname_format(split_idx + 1:end);
t_num_format    = regexp(fname_format,'image_%0\d+i','match');
fname_format2   = replace(fname_format2,t_num_format,'image_%d');
AgarInfoDir     = Parameters.Agar.dir;
Agar_file_format = Parameters.Agar.file_format;

% if exist([WriteMacro.TileDir '../StitchingFiji_coord'],"dir")
%     pathname_coord=[WriteMacro.TileDir '/../StitchingFiji_coord/' FijiInfo.FileBase '%03i' fiji_file_ext];
%     fprintf('coordinate path is used : %s\n',pathname_coord)
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Experiment_Fiji = Experiment;

% Initialize variables
all_coords_x = cell(num_mosaics_per_slice,1);
all_coords_y = cell(num_mosaics_per_slice,1);
MapIndex_Tot = cell(num_mosaics_per_slice,1);
for n = 1:num_mosaics_per_slice
    switch n
        case 1
            X_Tot = Experiment_Fiji.X_Tot;
            Y_Tot = Experiment_Fiji.Y_Tot;
            MI_Tot = Experiment_Fiji.MapIndex_Tot;
        case 2
            X_Tot = Experiment_Fiji.X_Tot_tilt;
            Y_Tot = Experiment_Fiji.Y_Tot_tilt;
            MI_Tot = Experiment_Fiji.MapIndex_Tot_tilt;
    end
    all_coords_x{n,1} = zeros([size(X_Tot) nslices{n,1}]);
    all_coords_y{n,1} = zeros([size(Y_Tot) nslices{n,1}]);
    MapIndex_Tot{n,1} = zeros([size(MI_Tot) nslices{n,1}]);
end

% Read coordinates (normal)

for n = 1:num_mosaics_per_slice % Loop between normal and tilted incidence
    for ssind=1:nslices{n,1}
        switch n
            case 1
                currslice = sliceidx(ssind);
                fprintf('Slice %d\n',currslice);
            case 2
                currslice = sliceidx_tilt(ssind);
                fprintf('Slice %d (tilted)\n',currslice);
        end
        % Loop z tiles (should always have Z_tile==1 though)
        for zz = 1:Experiment_Fiji.Z_tile

            %         % Get tile layout (tile numbers)
            %         MapIndex = Experiment_Fiji.MapIndex_Tot(:,:,zz) + (currslice-1)*Experiment_Fiji.TilesPerSlice;
            %
            %         % Get "Basic" tile coordinates (fixed overlap/uniform positions)
            %         XX = squeeze(Experiment_Fiji.X_Tot(:,:,zz));
            %         YY = squeeze(Experiment_Fiji.Y_Tot(:,:,zz));
            %
            %         % Mtxs for Fiji optimal tile coordinates
            %         X1 = XX; %zeros(size(MapIndex)); %XX;
            %         Y1 = YY; %zeros(size(MapIndex)); %YY;
            switch n
                case 1 % Normal incidence
                    % MapIndex = Experiment_Fiji.MapIndex_Tot(:,:,zz)+(currslice-1)*Experiment_Fiji.TilesPerSlice;
                    MapIndex = Experiment_Fiji.MapIndex_Tot(:,:,zz); % New naming format re-cycles through tile numbers for every slice
                    XX = squeeze(Experiment_Fiji.X_Tot(:,:,zz));
                    YY = squeeze(Experiment_Fiji.Y_Tot(:,:,zz));
                    X1 = zeros(size(MapIndex));
                    Y1 = zeros(size(MapIndex));
                case 2 % Tilted incidence
                    % MapIndex = Experiment_Fiji.MapIndex_Tot_tilt(:,:,zz)+(currslice-1)*Experiment_Fiji.TilesPerSlice_tilt;
                    MapIndex = Experiment_Fiji.MapIndex_Tot_tilt(:,:,zz); % New naming format re-cycles through tile numbers for every slice
                    XX = squeeze(Experiment_Fiji.X_Tot_tilt(:,:,zz));
                    YY = squeeze(Experiment_Fiji.Y_Tot_tilt(:,:,zz));
                    X1 = zeros(size(MapIndex));
                    Y1 = zeros(size(MapIndex));
            end
            

            % Read Fiji stitching coordinates file
            %  (Can use 'registered' from compute_overlap, or simply
            %       TileCoordinates.txt from manual stitching in GUI)
            fid = fopen(sprintf(pathname_coord{n,1},currslice),'r');

            % Skip 1st 4 lines
            for ii=1:4
                kk = fgets(fid);
            end

            % Read coords
            ii=0;
            while ~feof(fid)
                ii = ii +1;
                kk = fgets(fid);
                scan_string = [sprintf(fname_format1,mosaic_nums{n,1}(ssind)),fname_format2,'; ; (%f, %f)'];
                img(:,ii) = sscanf(kk,scan_string);
            end
            fclose(fid);

            % Extract coords and save to matrices
            MapIndex2 = -1*ones(size(MapIndex));
            coord = (img');
            warningcounter = 0;
            for ii=1:size(coord,1)
                [row,col] = find(MapIndex==coord(ii,1));
                if isempty(row) && isempty(col)
                    warningcounter = warningcounter+1;
                end
                MapIndex2(row,col) = coord(ii,1);
                % round coords to nearest integer
                Y1(row,col) = round(coord(ii,2));
                X1(row,col) = round(coord(ii,3));
            end

            if warningcounter>0
                warning('Tile number mismatch detected between MapIndex and Fiji coord text file');
                disp(' -- Make sure ReadCoord.ReadSliceID is correct. ');
                disp(' -- Make sure tile nbs in MapIndex match those in text file.');
            end

            % agarose tile
            MapIndex_agar = false(size(MapIndex));
            for ii = 1:length(MapIndex(:))
                load([AgarInfoDir,'/',sprintf(Agar_file_format,mosaic_nums{n,1}(ssind),ii)],'tileinfo');
                if(tileinfo(3) == 0) % if tileinfo(3) is 0 which is agar, then skip
                    MapIndex_agar(ii) = true;
                end
            end

            % Set min. coord value in X,Y = 1; and find nans;
            clear img coord;

            X        = X1;
            Y        = Y1;

            X           = X-min(min(X))+1;
            Y           = Y-min(min(Y))+1;
            indexNaN    = (MapIndex2==-1) | MapIndex_agar;


            X(indexNaN) = NaN;
            Y(indexNaN) = NaN;

            % Uniform the 'Origin' coordinate
            X = X - X(round(size(X,1)/2),round(size(X,2)/2));
            Y = Y - Y(round(size(Y,1)/2),round(size(Y,2)/2));

            % Store coords in 3d matrix
            all_coords_x{n,1}(:,:,ssind) = X;
            all_coords_y{n,1}(:,:,ssind) = Y;
            MapIndex_Tot{n,1}(:,:,ssind) = MapIndex2;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normal incidence coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Experiment_Fiji.X_Tot = all_coords_x{1,1};
Experiment_Fiji.Y_Tot = all_coords_y{1,1};
Experiment_Fiji.MapIndex_Tot = MapIndex_Tot{1,1};

% Method 1: Set XY_Mean coords to the median of all slices (omit NaNs)
Method_1.X_Mean = round(squeeze(median(all_coords_x{1,1},3,'omitnan')));
Method_1.Y_Mean = round(squeeze(median(all_coords_y{1,1},3,'omitnan')));
Method_1.X_std = squeeze(std(all_coords_x{1,1},[],3,'omitnan'));
Method_1.Y_std = squeeze(std(all_coords_y{1,1},[],3,'omitnan'));
Experiment_Fiji.Method_1 = Method_1;

Experiment_Fiji.X_Mean = round(squeeze(median(all_coords_x{1,1},3,'omitnan')));
Experiment_Fiji.Y_Mean = round(squeeze(median(all_coords_y{1,1},3,'omitnan')));


% Method 2: Set XY_Mean coords from median offset and median step size
X_tile = Experiment_Fiji.X_tile;
Y_tile = Experiment_Fiji.Y_tile;

step_x          = median(   diff(all_coords_x{1,1}, 1, 1),           'all','omitnan'); step_x = round(step_x);
matrix_x        = repmat(   [1: step_x : step_x * X_tile].',    [1 Y_tile]);
offset_X        = median(   diff(all_coords_x{1,1}, 1, 2),           'all','omitnan'); offset_X = round(offset_X);
offset_map_X    = repmat(   [1:Y_tile]-1,                       [X_tile,1])         *offset_X;
Method_2.X_Mean = matrix_x + offset_map_X;
Method_2.step_x = step_x;
Method_2.offset_X = offset_X;

step_y          = median(   diff(all_coords_y{1,1}, 1, 2),           'all','omitnan'); step_y = round(step_y);
matrix_y        = repmat(   [1: step_y : step_y * Y_tile],      [X_tile 1]);
offset_Y        = median(   diff(all_coords_y{1,1}, 1, 1),           'all','omitnan'); offset_Y = round(offset_Y);
offset_map_Y    = repmat(   [1:X_tile].'-1,                     [1,Y_tile])         *offset_Y;
Method_2.Y_Mean = matrix_y + offset_map_Y;
Method_2.step_y = step_y;
Method_2.offset_Y = offset_Y;

Experiment_Fiji.Method_2 = Method_2;

switch Method
    case 1
        Experiment_Fiji.X_Mean = Experiment_Fiji.Method_1.X_Mean;
        Experiment_Fiji.Y_Mean = Experiment_Fiji.Method_1.Y_Mean;
        Experiment_Fiji.X_std = Experiment_Fiji.Method_1.X_std;
        Experiment_Fiji.Y_std = Experiment_Fiji.Method_1.Y_std;
    case 2
        Experiment_Fiji.X_Mean = Experiment_Fiji.Method_2.X_Mean;
        Experiment_Fiji.Y_Mean = Experiment_Fiji.Method_2.Y_Mean;
end


% Find tiles that were NaN for all slices - exclude in Mosaic)
Experiment_Fiji.Nantile_Prc = sum(isnan(Experiment_Fiji.X_Tot),3)/size(Experiment_Fiji.X_Tot,3);
nantiles = Experiment_Fiji.Nantile_Prc>0.95; % Moved this up from below
Experiment_Fiji.X_Mean(nantiles) = nan;
Experiment_Fiji.Y_Mean(nantiles) = nan;

% nantiles = Experiment_Fiji.Nantile_Prc>0.95;

nantiles = isnan(Experiment_Fiji.X_Mean);

Experiment_Fiji.MapIndex_Fiji = Experiment_Fiji.MapIndex_Tot_offset;
Experiment_Fiji.MapIndex_Fiji(nantiles) = -1;

fprintf('Y std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.Y_std(:)),median(Experiment_Fiji.Y_std(:),'omitnan'));
fprintf('X std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.X_std(:)),median(Experiment_Fiji.X_std(:),'omitnan'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tilted incidence coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    Experiment_Fiji.X_Tot_tilt = all_coords_x{2,1};
    Experiment_Fiji.Y_Tot_tilt = all_coords_y{2,1};
    Experiment_Fiji.MapIndex_Tot_tilt = MapIndex_Tot{2,1};

    % Method 1: Set XY_Mean coords to the median of all slices (omit NaNs)
    Method_1.X_Mean_tilt = round(squeeze(median(all_coords_x{2,1},3,'omitnan')));
    Method_1.Y_Mean_tilt = round(squeeze(median(all_coords_y{2,1},3,'omitnan')));
    Method_1.X_std_tilt = squeeze(std(all_coords_x{2,1},[],3,'omitnan'));
    Method_1.Y_std_tilt = squeeze(std(all_coords_y{2,1},[],3,'omitnan'));
    Experiment_Fiji.Method_1_tilt = Method_1;

    Experiment_Fiji.X_Mean_tilt = round(squeeze(median(all_coords_x{2,1},3,'omitnan')));
    Experiment_Fiji.Y_Mean_tilt = round(squeeze(median(all_coords_y{2,1},3,'omitnan')));


    % Method 2: Set XY_Mean coords from median offset and median step size
    X_tile = Experiment_Fiji.X_tile_tilt;
    Y_tile = Experiment_Fiji.Y_tile_tilt;

    step_x          = median(   diff(all_coords_x{2,1}, 1, 1),           'all','omitnan'); step_x = round(step_x);
    matrix_x        = repmat(   [1: step_x : step_x * X_tile].',    [1 Y_tile]);
    offset_X        = median(   diff(all_coords_x{2,1}, 1, 2),           'all','omitnan'); offset_X = round(offset_X);
    offset_map_X    = repmat(   [1:Y_tile]-1,                       [X_tile,1])         *offset_X;
    Method_2.X_Mean_tilt = matrix_x + offset_map_X;
    Method_2.step_x_tilt = step_x;
    Method_2.offset_X_tilt = offset_X;

    step_y          = median(   diff(all_coords_y{2,1}, 1, 2),           'all','omitnan'); step_y = round(step_y);
    matrix_y        = repmat(   [1: step_y : step_y * Y_tile],      [X_tile 1]);
    offset_Y        = median(   diff(all_coords_y{2,1}, 1, 1),           'all','omitnan'); offset_Y = round(offset_Y);
    offset_map_Y    = repmat(   [1:X_tile].'-1,                     [1,Y_tile])         *offset_Y;
    Method_2.Y_Mean_tilt = matrix_y + offset_map_Y;
    Method_2.step_y_tilt = step_y;
    Method_2.offset_Y_tilt = offset_Y;

    Experiment_Fiji.Method_2_tilt = Method_2;

    switch Method
        case 1
            Experiment_Fiji.X_Mean_tilt = Experiment_Fiji.Method_1_tilt.X_Mean_tilt;
            Experiment_Fiji.Y_Mean_tilt = Experiment_Fiji.Method_1_tilt.Y_Mean_tilt;
            Experiment_Fiji.X_std_tilt = Experiment_Fiji.Method_1_tilt.X_std_tilt;
            Experiment_Fiji.Y_std_tilt = Experiment_Fiji.Method_1_tilt.Y_std_tilt;
        case 2
            Experiment_Fiji.X_Mean_tilt = Experiment_Fiji.Method_2_tilt.X_Mean_tilt;
            Experiment_Fiji.Y_Mean_tilt = Experiment_Fiji.Method_2_tilt.Y_Mean_tilt;
    end


    % Find tiles that were NaN for all slices - exclude in Mosaic)
    Experiment_Fiji.Nantile_Prc_tilt = sum(isnan(Experiment_Fiji.X_Tot_tilt),3)/size(Experiment_Fiji.X_Tot_tilt,3);
    nantiles = Experiment_Fiji.Nantile_Prc>0.95; % Moved this up from below
    Experiment_Fiji.X_Mean(nantiles) = nan;
    Experiment_Fiji.Y_Mean(nantiles) = nan;

    % nantiles = Experiment_Fiji.Nantile_Prc>0.95;

    nantiles = isnan(Experiment_Fiji.X_Mean_tilt);

    Experiment_Fiji.MapIndex_Fiji_tilt = Experiment_Fiji.MapIndex_Tot_offset_tilt;
    Experiment_Fiji.MapIndex_Fiji_tilt(nantiles) = -1;

    fprintf('Y std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.Y_std_tilt(:)),median(Experiment_Fiji.Y_std_tilt(:),'omitnan'));
    fprintf('X std = %.2f/%.2f (max/median [px]); \n', max(Experiment_Fiji.X_std_tilt(:)),median(Experiment_Fiji.X_std_tilt(:),'omitnan'));
end

% Save Experiment_Fiji struct to .mat file (to load during Mosaic)
if saveflag==1
    save(fout,'Experiment_Fiji','FijiInfo','Experiment');
end








