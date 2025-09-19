function Experiment = ReadTLSSLogFile(ParameterFile, fid, sliceid)

% Read parameters from TLSS log file.
% 
% PSOCT scan on Octopus system.
% TLSS log file is from stage computer.
% 
% fid = importdata([pathname filename_log],'\n');
% sliceid = current slice number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % NB - Edit notes, 04/21/2025 % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added: 
%   Extra read-ins for multiple Z stage configurations
%   Extra read-ins for motorized rotational stage for tilted illumination
% Removed:
%   Offset of 1 on the tile index since mosaic_array(TILE_IND, : ,MOSAIC_IND) = 
%      [IND_X, IND_Y, IND_Z, POS_X, POS_Y, POS_ZL, POS_ZR, TIME_STAMP] was
%      removed
% NOTE: 
%   The two Z-stage coordinates ZL and ZR are added


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % RJJ - Edit notes, 12/22/2020 % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Replaced:
%   row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
% with:
%   row_idx = find(contains(fid,FormatSpec));
% 
% NOTE: must keep str2num(str(a+1:end-2))
%   (Using str2double() does not work, returns NaN)
% 
% Added:
%  - Instantiate Experiment struct
%  - Instantiate coord matrix
%  - Experiment.TilesPerSlice field
%  - Experiment.SliceThickness field
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Mosaic Parameters % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(ParameterFile);

%struct to save parameters to
Experiment = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read tile step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FormatSpec = 'mosaic_step = [';
a = length(FormatSpec);
% row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
row_idx = find(contains(fid,FormatSpec));
str = fid{row_idx(1)};
% str = fid{row_idx};
Tile_dim = str2num(str(a+1:end-2));
Experiment.X_step = Tile_dim(1);
Experiment.Y_step = Tile_dim(2);
Experiment.Z_step = Tile_dim(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read number of tiles in x and y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FormatSpec = 'mosaic_ntiles = [';
a = length(FormatSpec);
% row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
row_idx = find(contains(fid,FormatSpec));
str = fid{row_idx(1)};
Tile_num = str2num(str(a+1:end-2));
Experiment.X_tile = Tile_num(1);
Experiment.Y_tile = Tile_num(2);
Experiment.Z_tile = Tile_num(3);
Total_tile = prod(Tile_num);
Experiment.TilesPerSlice = Total_tile;

if(strcmpi(Scan.TiltedIllumination,'Yes'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read tilted tile step size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FormatSpec = 'mosaic_step_tilt = [';
    a = length(FormatSpec);
    % row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
    row_idx = find(contains(fid,FormatSpec));
    str = fid{row_idx(1)};
    % str = fid{row_idx};
    Tile_dim = str2num(str(a+1:end-2));
    Experiment.X_step_tilt = Tile_dim(1);
    Experiment.Y_step_tilt = Tile_dim(2);
    Experiment.Z_step_tilt = Tile_dim(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read number of tilted tiles in x and y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FormatSpec = 'mosaic_ntiles_tilt = [';
    a = length(FormatSpec);
    % row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
    row_idx = find(contains(fid,FormatSpec));
    str = fid{row_idx(1)};
    Tile_num = str2num(str(a+1:end-2));
    Experiment.X_tile_tilt = Tile_num(1);
    Experiment.Y_tile_tilt = Tile_num(2);
    Experiment.Z_tile_tilt = Tile_num(3);
    Total_tile = prod(Tile_num);
    Experiment.TilesPerSlice_tilt = Total_tile;
    Total_tile = Experiment.TilesPerSlice + Experiment.TilesPerSlice_tilt;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read mosaic Coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FormatSpec = 'mosaic_array';
% row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
row_idx = find(contains(fid,FormatSpec));

% tile_index = ((Total_tile)*(sliceid- 1)+1:(Total_tile)*sliceid);
tile_index = 1 + ((Total_tile)*(sliceid- 1)+1:(Total_tile)*sliceid);

 % Add offset of 1 b/c first instance of 'mosaic_array' is:
    % mosaic_array(TILE_IND, : ,MOSAIC_IND) = 
     % [IND_X, IND_Y, IND_Z, POS_X, POS_Y, POS_ZL, POS_ZR, TIME_STAMP]

% Now, loop through all tiles in the slice
num_coord = 10;
line_format = ['mosaic_array(%f,:,%f) = [%f, %f, %f, %f, %f, %f, %f, %f, %f];'];

coord = zeros(length(tile_index),num_coord); 
for ii=1:length(tile_index)
    % ii
    str = fid{row_idx(tile_index(ii))};
    coord(ii,:) = sscanf(str,line_format);
end

%Save tile coordinates to Experiment
ii_tilt = 0;
for ii = 1:size(coord,1)
    if(strcmpi(Scan.TiltedIllumination,'Yes') && ii > Experiment.TilesPerSlice)
        ii_tilt = ii_tilt + 1;
        Experiment.MapIndex_Tot_tilt(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,1); 
        Experiment.X_Tot_tilt(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,6);
        Experiment.Y_Tot_tilt(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,7);

        % Use normal incidence Z values for the tilted acquisition (an additional Z offset is required to allow for the physical tilting) 
        if(strcmpi(Scan.ZStageConfig,'Single'))
            Experiment.Z_Tot_tilt(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(1,8);
        elseif(strcmpi(Scan.ZStageConfig,'Stacked'))
            Experiment.Z_Tot_tilt(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(1,8) + coord(1,9);
        end
    else
        Experiment.MapIndex_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,1); 
        Experiment.X_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,6);
        Experiment.Y_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,7);
        if(strcmpi(Scan.ZStageConfig,'Single'))
            Experiment.Z_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,8);
        elseif(strcmpi(Scan.ZStageConfig,'Stacked'))
            Experiment.Z_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,8) + coord(ii,9);
        end
    end
    
    
end

%MapIndex is tile number corresponding to each position (x,y)
for ii = 1:size(Experiment.MapIndex_Tot,3)
    Experiment.MapIndex_Tot(:,:,ii)= squeeze(Experiment.MapIndex_Tot(:,:,ii));
end
    
%Normalize so x,y pixel positions start at 1
Experiment.X_Tot = Experiment.X_Tot - min(Experiment.X_Tot(:)) +1;
Experiment.Y_Tot = Experiment.Y_Tot - min(Experiment.Y_Tot(:)) +1;
Experiment.Z_Tot = Experiment.Z_Tot - min(Experiment.Z_Tot(:)) +1;

if(strcmpi(Scan.TiltedIllumination,'Yes'))
    for ii = 1:size(Experiment.MapIndex_Tot_tilt,3)
        Experiment.MapIndex_Tot_tilt(:,:,ii)= squeeze(Experiment.MapIndex_Tot_tilt(:,:,ii));
    end

    %Normalize so x,y pixel positions start at 1
    Experiment.X_Tot_tilt = Experiment.X_Tot_tilt - min(Experiment.X_Tot_tilt(:)) +1;
    Experiment.Y_Tot_tilt = Experiment.Y_Tot_tilt - min(Experiment.Y_Tot_tilt(:)) +1;
    Experiment.Z_Tot_tilt = Experiment.Z_Tot_tilt - min(Experiment.Z_Tot_tilt(:)) +1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read vibratome slice thickness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FormatSpec = 'mosaic_cut_size_z = ';
a = length(FormatSpec);
% row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
row_idx = find(contains(fid,FormatSpec));
str = fid{row_idx(1)};
Slice_thickness = str2num(str(a+1:end));
Experiment.SliceThickness = Slice_thickness;


end
