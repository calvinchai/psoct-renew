function Experiment = ReadTLSSLogFile(fid, sliceid)

% Read parameters from TLSS log file.
% 
% PSOCT scan on Octopus system.
% TLSS log file is from stage computer.
% 
% fid = importdata([pathname filename_log],'\n');
% sliceid = current slice number
 

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read mosaic Coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FormatSpec = 'mosaic_array';
% row_idx = find(~cellfun('isempty',strfind(fid,FormatSpec)));
row_idx = find(contains(fid,FormatSpec));
tile_index = 1+((Total_tile)*(sliceid- 1)+1:(Total_tile)*sliceid);
%  Add offset of 1 b/c first instance of 'mosaic_array' is:
%    % mosaic_array(TILE_IND, : ,MOSAIC_IND) = 
%      [IND_X, IND_Y, IND_Z, POS_X, POS_Y, POS_ZL, POS_ZR, TIME_STAMP]

%Now, loop through all tiles in the slice
coord = zeros(length(tile_index),10);
for ii=1:length(tile_index)
    ii
    str = fid{row_idx(tile_index(ii))};
    coord(ii,:) = sscanf(str,['mosaic_array(%f,:,%f) = [%f, %f, %f, %f, %f, %f, %f, %f];']);
end

%Save tile coordinates to Experiment
for ii = 1:size(coord,1)
    Experiment.MapIndex_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,1);
    Experiment.X_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,6);
    Experiment.Y_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,7);
    Experiment.Z_Tot(coord(ii,4),coord(ii,3),coord(ii,5)) = coord(ii,8);
end

%MapIndex is tile number corresponding to each position (x,y)
for ii = 1:size(Experiment.MapIndex_Tot,3)
    Experiment.MapIndex_Tot(:,:,ii)= squeeze(Experiment.MapIndex_Tot(:,:,ii));
end
    
%Normalize so x,y pixel positions start at 1
Experiment.X_Tot = Experiment.X_Tot - min(Experiment.X_Tot(:)) +1;
Experiment.Y_Tot = Experiment.Y_Tot - min(Experiment.Y_Tot(:)) +1;
Experiment.Z_Tot = Experiment.Z_Tot - min(Experiment.Z_Tot(:)) +1;


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
