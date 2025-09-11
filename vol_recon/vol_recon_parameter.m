clear
addpath('/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/vol_recon_proj_NB_2025');

% t = tic;
% next step sets path for parameter file put it in your processed directory
ParameterFile  = ['/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData_NB/Parameters.mat'];
%P = whos('-file',ParameterFile);
if exist(ParameterFile,'file')
    load(ParameterFile);
    load(Parameters.ExpBasic);
    load(Parameters.ExpFiji);
    warning(' --- ParameterFile already exists. Script is stopped to avoid overwriting. --- ');
    return; 
end
%% Aquisition Parameters 
% This portion is mostly for record keeping
Scan.SampleName          = 'I80 Premotor Slab';% name of your sample
Scan.Date                = '05_13_2025';% date of acquisition
Scan.OCTOperator         = 'Chris';% who acquired the data
Scan.PostProcessOperator = 'Nate';% who is processing the data
Scan.ImagingNotes        = 'I80 human hemi slab with premotor cortex';% put any additional notes
Scan.SaveMethod          = 'complex';  % either "spectral", "complex", or "surfacefinding"
Scan.TiltedIllumination  = 'Yes'; % either "Yes" or "No"

% Detail the input file name format for the raw data
if contains('spectral',Scan.SaveMethod)
    Scan.FilePrefix      = {'mosaic_*_normal_*.nii','mosaic_*_tilted_*.nii'}; 
    Scan.FileNameFormat  = 'mosaic_%03i_image_%04i_x_%06i_y_%06i_z_%06i_r_%06i_normal_spectral_%04i.nii'; % Only used as a template, do not need 'tilted variant'
    Scan.CropMethod      = 'none'; % either "none", "surface", or "focus"
elseif contains('surfacefinding',Scan.SaveMethod)
    % Scan.FilePrefix      = {'test_processed_'};
    % Scan.FileNameFormat  = {'unused'};
    % Scan.CropMethod      = 'surface'; % either "none", "surface", or "focus"
elseif contains('complex',Scan.SaveMethod)
    Scan.FilePrefix      = {'cropped_focus'}; % Strings used for [modality]
    Scan.FileNameFormat  = 'mosaic_%03i_image_%04i_processed_[modality].nii'; % Only used as a template, do not need 'tilted variant'
    Scan.CropMethod      = 'focus'; % either "none", "surface", or "focus"
end
if(strcmpi(Scan.CropMethod,'focus') == 1)
    Scan.Focus_CropStart = 150; % Starting cropping 150 pixels in Z above the 2D focus map
    % Scan.FocusFile = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic1thru7.nii';
    Scan.FocusFile = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic9plus.nii';
    if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
        % Scan.FocusFile_tilt = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic2thru8.nii';
        Scan.FocusFile_tilt = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/Focus/focus_mosaic10plus.nii';
    end
end
Scan.Objective           = '10x'; % what objective was used, this is normally 10x
Scan.System              = 'Telesto'; IsTelesto = strcmpi(Scan.System,'Telesto');% 
Scan.Bath_Solution       = 'Mineraloil'; % what imaging medium was

Scan.RawDataDir          = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData';
Scan.ProcessDir          = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData_NB';
Scan.TLSS_log            = 'TL_serial_shell_v31_20250509_092640_log.txt'; % used for creating ExperimentBasic
Scan.ZStageConfig        = 'Stacked'; % either "Single" or "Stacked"
Scan.First_Tile          = 1; 

Scan.Thickness           = 500;   % um  How thick sections were
Scan.FoV                 = 3500;  % um  Get this from the acquisition software
Scan.NbPixels            = 350;   % pix  Get this from the acquisition software
Scan.StepSize            = 2800;   % um  Get this from the log file
if strcmpi(Scan.TiltedIllumination,'Yes')
    Scan.FoV_tilt             = 2000;  % um  Along shorter (X) axis, confined by physical tilting geometry
    Scan.NbPixels_tilt        = 200;   % pix  Along shorter (X) axis
    Scan.StepSize_tilt        = 1600;   % pix  Along shorter (X) axis
    Scan.First_Tile_tilt      = 1;
end
Scan.CropDepth           = 500;   % pix  size(niftiread([Scan.RawDataDir '/test_processed_001_cropped.nii']));
                                  % if processing spectral this is just how
                                  % much of the dBI volume you saved when
                                  % cropping
                               
in_planeRes              = 0.010;  % mm % Scan.FoV/Scan.NbPixels/1000;
thru_planeRes            = 0.0025; % mm
Scan.Resolution          = [in_planeRes in_planeRes thru_planeRes];

save(ParameterFile,'Scan');% creates parameter file and saves file portion scan

%% Folder & File names
% This sets up the folders where all the following steps will save to
ProcDir     = Scan.ProcessDir;
RawDataDir  = Scan.RawDataDir;

Proc3D          = [ProcDir '/Processed3D'];
% Enface_Fiji     = [ProcDir '/Enface_Tiff'];
% Enface_Fiji     = RawDataDir;
Enface_Fiji     = [ProcDir '/nii2D_transpose'];
Enface_MAT      = [ProcDir '/Enface_MAT'];
StitchingFiji   = [ProcDir '/StitchingFiji'];
StitchingBasic  = [ProcDir '/Stitching'];
StackNii        = [ProcDir '/StackNii'];
MapIndex        = [ProcDir '/MapIndex'];

TLSSLog         = [ProcDir '/' Scan.TLSS_log];
ExpBasic        = [ProcDir '/Experiment.mat'];
ExpFiji         = [ProcDir '/Experiment_Fiji.mat'];

%% Detail the input/output file name formats for the processed3D data
Processed3D.indir         = RawDataDir; % Spectral or Complex files
Processed3D.input_format  = Scan.FileNameFormat; % Spectral or Complex file format
Processed3D.outdir        = Proc3D;
Processed3D.output_format = 'mosaic_%03i_image_%04i_processed_[modality].nii'; % '[modality]_mosaic_%03i_image_%04i.nii';
Processed3D.outputstr     = {'dBI3D','R3D','O3D','biref'};%,'mus','dBI_crop','R3D_crop','O3D_crop'}; % Strings used for [modality]
Processed3D.save          = Processed3D.outputstr([4]); % ([1,2,3,4]);
Processed3D.dispComp      = '/autofs/cluster/octdata2/users/Hui/tools/dg_utils/spectralprocess/dispComp/mineraloil_LSM03/dispersion_compensation_LSM03_mineraloil_20240829/LSM03_mineral_oil_placecorrectionmeanall2.dat';
Processed3D.RollOff       = '/cluster/octdata2/users/TelestoCalibration/RollOff_Spectrum.mat';

% Check enface file name formats
fprintf('Output file naming format for Processed3D data \n\n');
fprintf(' - TEMPLATE: %s \n',Processed3D.output_format);
fprintf(' - EXAMPLES:  \n');
for i = 1:length(Processed3D.save)
    example_filename = replace(sprintf(Processed3D.output_format,1,1),'[modality]',Processed3D.save{i});
    % example_filename = replace(sprintf(Enface.file_format,1),'[modality]',Enface.save{i});
    fprintf(' -           %s \n',example_filename);
end

save(ParameterFile,'Processed3D','-append');


%% Detail the input/output file name formats for the 2D enface images
Enface.indir          = Scan.RawDataDir; % [ProcDir, '/MAT2d']; % [RawDataDir];
Enface.inputstr       = {'aip','mip','retardance','orientation','surface_finding','biref','mus'}; % Strings used for [modality]
Enface.save           = Enface.inputstr([1,2,3,4,5,6]); % ,5,6,7]); % inputstr [modality] to save 
Enface.input_format   = 'mosaic_%03i_image_%04i_processed_[modality].nii';
Enface.outdir         = Enface_Fiji;
Enface.outputstr      = {'aip','mip','retardance','orientation','surface_finding','biref','mus'}; % Strings used for [modality]
Enface.save_tiff      = Enface.outputstr([1,2,3,4,6]); % outputstr [modality] to save 
Enface.output_format  = 'mosaic_%03i_image_%04i_processed_[modality].nii'; % Needs to match file type used for WriteMacro_Serpentine.m

% Check enface file name formats
fprintf('Input file naming format for Enface and Mosaic2D \n\n');
fprintf(' - TEMPLATE: %s \n',Enface.input_format);
fprintf(' - EXAMPLES:  \n');
for i = 1:length(Enface.save)
    example_filename = replace(sprintf(Enface.input_format,1,1),'[modality]',Enface.save{i});
    % example_filename = replace(sprintf(Enface.file_format,1),'[modality]',Enface.save{i});
    fprintf(' -           %s \n',example_filename);
end

save(ParameterFile,'Enface','-append');
% Macro_ijm is defined later.

%% ExperimentBasic
if ~exist(ExpBasic,'file')
    [ ExperimentBasic ] = PrepareBasicExperiment( ParameterFile, TLSSLog, ExpBasic );
    save(ParameterFile,'ExperimentBasic','-append');
end
load(ExpBasic);
% Parameters.SliceID = [4]; % Mosaic 7, 8 (to run birefringence calculation)
% Parameters.SliceID = [5:7]; % Mosaic 9, 10, 11, 12, 13, 14 (to run birefringence calculation)
% Parameters.SliceID = [4:11]; % Mosaic 7 to 22 
% Parameters.SliceID = [13:19]; % Mosaic 25 to 38 (to run stitching) - stopped at mosaic_030_image_400
Parameters.SliceID = [5:11,13:26]; % Slices to try 3D stitching on (Kaidong)

mosaic_nums = [];
tile_nums = [];
if strcmpi(Scan.TiltedIllumination,'Yes')
    for n = 1:length(Parameters.SliceID)
        tmp_mosaic_nums = [repelem(2*Parameters.SliceID(n) - 1,ExperimentBasic.TilesPerSlice),... % Normal
            repelem(2*Parameters.SliceID(n),ExperimentBasic.TilesPerSlice_tilt)]; % Tilted
        tmp_tile_nums = [1:ExperimentBasic.TilesPerSlice,... % Normal
            1:ExperimentBasic.TilesPerSlice_tilt]; % Tilted
        mosaic_nums = [mosaic_nums,tmp_mosaic_nums];
        tile_nums = [tile_nums,tmp_tile_nums];
    end
else
    mosaic_nums = repelem(1:Parameters.SliceID(end),ExperimentBasic.TilesPerSlice);
    tile_nums = repmat(1:ExperimentBasic.TilesPerSlice,1,Parameters.SliceID(end));
end
Parameters.MosaicID = mosaic_nums;
Parameters.TileID  = tile_nums;
clear n mosaic_nums tile_nums tmp_mosaic_nums tmp_tile_nums

Parameters.OrientationSign = 1;
Parameters.OrienOffset     = 0;

Parameters.XPixClip = 18;% how many pixels to clip in x
Parameters.YPixClip = 0; % how many pixels to clip in y 
Parameters.transpose = false; % if artifact is not on the top side from imshow() or tile tiff, then do transpose. Typically, this is false.

Parameters.ExpFiji = ExpFiji;% where experiment Fiji mat file is saved
Parameters.ExpBasic= ExpBasic; % where exp basic file is saved

save(ParameterFile,'Parameters','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 Begin Processing & Stitching Pipeline                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process Raw Data
if(0 == 1)
    if contains ('spectral',Scan.SaveMethod) % Process Spectral Raw data --> dBI3D, R3D, O3D, and 2D mat files
        istart = 165; % Starting Z-index for SurfaceFinding

        S2P_Spectral_to_Processed(ParameterFile,istart);
    elseif contains ('complex',Scan.SaveMethod) % Process Complex Processed data (Jones1, Jones2) --> Processed3D files (dBI3D, R3D, O3D, mus)
        Complex_to_Processed3D(ParameterFile);
    else
    end
end
%% Quick view concatenated slices 
% Used to get grayscale windowing values for generating tiff images (enter below)
if ~exist('Parameters','var');load(ParameterFile);end
view_modality = 'aip'; % aip,mip,ret,biref,mus

mosaic_num = 9;
mosaic_dim = [ExperimentBasic.Y_tile, ExperimentBasic.X_tile];
transpose_flag = 1;
% Only works for Mosaic 1 at the moment
[c_enface,c,~] = get_enface_concatenate(Enface.indir,Enface.input_format,view_modality,mosaic_num,mosaic_dim,1,transpose_flag);
view3D(c_enface)
figure; histogram(c_enface)


%%
Parameters.AipGrayRange =   [47   68];
Parameters.MipGrayRange =   [58   85];
Parameters.RetGrayRange =   [30   50];
Parameters.BirefGrayRange = [0    1e-4];
Parameters.musGrayRange =   [0.7  3];
% Parameters.Surfacerange = [130      210];

save(ParameterFile,'Parameters','-append');

%% Agarose thresholding based on std (removing tiles with no tissue)

% Implement std check for first 1000 tiles. Determine threshold to
% separate agarose tiles from tissue tiles. Use Caroline's WriteMacro that
% operates based on file location

inspect_slice_num = 10;
modality = 'aip';

mosaic_num_normal = inspect_slice_num*2 - 1;
mosaic_num_tilted = inspect_slice_num*2;

std_slice = cell(2,1);
for n_tilt = 1:2
    switch n_tilt
        case 1
            title_string = 'Normal Illumination';
            idx = find(Parameters.MosaicID == mosaic_num_normal);
            idx_normal = idx';
        case 2
            title_string = 'Tilted Illumination';
            idx = find(Parameters.MosaicID == mosaic_num_tilted);
            idx_tilted = idx';
    end
    std_tmp = zeros(numel(idx),1);
    for n = 1:numel(idx)
        filename = replace(sprintf(Enface.input_format,Parameters.MosaicID(1,idx(n)),Parameters.TileID(1,idx(n))),'[modality]',modality);
        tmp_modality = niftiread(fullfile(Enface.indir,filename));
        std_tmp(n,1) = std2(tmp_modality);
    end
    std_slice{n_tilt,1} = std_tmp;
    figure(n_tilt);
    histogram(std_slice{n_tilt,1},'BinWidth',0.05);
    title(title_string);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Look at aip std tiles in intermediate range
low_range = 1.05; % 0.25;
high_range = 1.15; % 0.3; 

idx_normal_range = find(std_slice{1,1} > low_range & std_slice{1,1} < high_range);
idx_tilted_range = find(std_slice{2,1} > low_range & std_slice{2,1} < high_range);

% for i = 1:length(idx_normal_range)
% for i = 2
%     filename = replace(sprintf(Enface.input_format,Parameters.MosaicID(1,idx_normal(idx_normal_range(i,1))),Parameters.TileID(1,idx_normal(idx_normal_range(i,1)))),'[modality]',modality);
%     tmp_modality = niftiread(fullfile(Enface.indir,filename));
%     figure(5);
%     imagesc(tmp_modality);
%     colorbar;
%     title(sprintf('aip %d, std %d',idx_normal(i,1),std_slice{1,1}(idx_normal_range(i,1))));
%     pause(2);
% end
% for i = 1:length(idx_tilted)
for i = 5
    filename = replace(sprintf(Enface.input_format,Parameters.MosaicID(1,idx_tilted(idx_tilted_range(i,1))),Parameters.TileID(1,idx_tilted(idx_tilted_range(i,1)))),'[modality]',modality);
    tmp_modality = niftiread(fullfile(Enface.indir,filename));
    figure(5);
    imagesc(tmp_modality);
    colorbar;
    title(sprintf('aip %d, std %d',idx_tilted(i,1),std_slice{2,1}(idx_tilted_range(i,1))));
    pause(2);
end

% clear n n_tilt title_string std_tmp filename tmp_modality
%% Sample aip tiles to look at thresholding
mask_val = 56;
stitch_path = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/StitchingFiji';

for inspect_slice_num = Parameters.SliceID
    slice_string = sprintf('Slice %d',inspect_slice_num);
    aip_normal = load([stitch_path,sprintf('/AIP_slice%03i.mat',inspect_slice_num)]);
    aip_normal = aip_normal.MosaicFinal;
    aip_tilted = load([stitch_path,sprintf('/AIP_tilt_slice%03i.mat',inspect_slice_num)]);
    aip_tilted = aip_tilted.MosaicFinal;
    aip_normal_masked = aip_normal;
    aip_normal_masked(aip_normal_masked < mask_val) = 0;
    aip_tilted_masked = aip_tilted;
    aip_tilted_masked(aip_tilted_masked < mask_val) = 0;


    fig1 = figure(6);
    subplot(1,2,1);
    imagesc(aip_normal);
    title([slice_string,' AIP normal']);
    subplot(1,2,2);
    imagesc(aip_normal_masked);
    title([slice_string,' AIP normal masked']);
    colorbar;
    fig1.Position = [1,551,1500,400];

    fig2 = figure(7);
    subplot(1,2,1);
    imagesc(aip_tilted);
    title([slice_string,' AIP tilted']);
    subplot(1,2,2);
    imagesc(aip_tilted_masked);
    title([slice_string,' AIP tilted masked']);
    colorbar;
    fig2.Position = [1,1,1500,400];
    pause(5);
end


%% Set Agar info

% Parameters.Agar.threshold_pct = 0.04;
% Parameters.Agar.threshold_std = 1.5; 
% Parameters.Agar.gm_aip = diff(Parameters.AipGrayRange)/4+Parameters.AipGrayRange(1);%(find thresh at mid-range intensity of graymatter which is roughly at lower half of AipGrayRange)

Parameters.Agar.aip_threshold = 56;


Parameters.Agar.dir = MapIndex;
Parameters.Agar.file_format = 'agar_status_mosaic_%03i_image_%04i.mat';
save(ParameterFile,'Parameters','-append');

% % use std on mat2gray(I,aipgrayscale) may yeild a more distinct result, std thresh may have variation then
% [ agar_status ] = view_tile_aip( ParameterFile, [7 5] );


%% WriteMacro & Read Registered Coordinates
WriteMacro.Modality    = 'aip'; 

if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    WriteMacro.SpiralOrigin = round([ExperimentBasic.X_tile_tilt,ExperimentBasic.Y_tile_tilt]/2);
end
WriteMacro.SpiralOrigin = round([ExperimentBasic.X_tile,ExperimentBasic.Y_tile]/2);

WriteMacro.SliceID     = Parameters.SliceID;
WriteMacro.TileDir     = Enface_Fiji;
WriteMacro.FijiOutDir  = StitchingFiji;
WriteMacro.SaveMacro   = [StitchingFiji '/Macro_' WriteMacro.Modality '.ijm'];
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    WriteMacro.SaveMacro_tilt   = [StitchingFiji '/Macro_' WriteMacro.Modality '_tilt.ijm'];
end
WriteMacro.AgarInfoDir = MapIndex;
save(ParameterFile,'WriteMacro','-append'); % FijiStep

ReadCoord.Method       = 2;     %select a method. 1=median; 2=step&offset
% ReadCoord.ReadSliceID  = Parameters.SliceID;
ReadCoord.ReadSliceID  = [6,7,9,10,11]; % Slices to use for normal stitching coordinates
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    ReadCoord.ReadSliceID_tilt  = [5,7]; % Slices to use for tilted stitching coordinates
end
ReadCoord.CoordDir     = StitchingFiji;
ReadCoord.SaveExp      = ExpFiji;
load(ExpBasic, 'ExperimentBasic');
ReadCoord.ExperimentBasic=ExperimentBasic;

FijiInfo.isRegistered  = 1;
FijiInfo.OutDir        = StitchingFiji;
FijiInfo.FileBase      = 'Mosaic_depth001_slice';
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    FijiInfo.FileBase_tilt  = 'Mosaic_tilt_depth001_slice';
end
save(ParameterFile,'FijiInfo', 'ReadCoord', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Save Enface Tiff
SaveTiff_Enface_script(ParameterFile); % Only for creating agarose info 

%% Create Macro.ijm file
WriteMacro_Serpentine(ParameterFile);

%% Run ImageJ (FIJI) without ssh for speed. e.g. tigervnc
load(ParameterFile,'WriteMacro')

t_stitch = tic;
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    mosaics_per_slice = 2;
else
    mosaics_per_slice = 1;
end

operating_system = 'Linux'; % 'Linux' or 'Windows';
for n = 1:mosaics_per_slice
    if(strcmpi(operating_system,'Linux'))
        % Linux commands do not work on Windows
        switch n
            case 1 % Normal incidence
                cmd = ['/autofs/cluster/octdata2/users/Hui/Fiji.app/ImageJ-linux64 --headless --console -macro ' WriteMacro.SaveMacro];
            case 2 % Tilted incidence
                cmd = ['/autofs/cluster/octdata2/users/Hui/Fiji.app/ImageJ-linux64 --headless --console -macro ' WriteMacro.SaveMacro_tilt];
        end
        system(cmd);

    elseif(strcmpi(operating_system,'Windows'))
        % Windows commands might also work on Linux
        fiji_path = '"C:\Program Files\Fiji.app\ImageJ-win64.exe"';
        switch n
            case 1 % Normal incidence
                macro_path = sprintf('"%s"', WriteMacro.SaveMacro);
            case 2 % Tilted incidence
                macro_path = sprintf('"%s"', WriteMacro.SaveMacro_tilt);
        end

        cmd = sprintf('%s --headless --console -macro %s', fiji_path, macro_path);
        system(cmd);
    end
end
for slice_n = Parameters.SliceID
    movefile([Enface_Fiji,'/',FijiInfo.FileBase,sprintf('%03i.txt',slice_n)],...
        [StitchingFiji,'/',FijiInfo.FileBase,sprintf('%03i.txt',slice_n)]);
    movefile([Enface_Fiji,'/',FijiInfo.FileBase,sprintf('%03i.registered.txt',slice_n)],...
        [StitchingFiji,'/',FijiInfo.FileBase,sprintf('%03i.registered.txt',slice_n)]);
    if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
        movefile([Enface_Fiji,'/',FijiInfo.FileBase_tilt,sprintf('%03i.txt',slice_n)],...
            [StitchingFiji,'/',FijiInfo.FileBase_tilt,sprintf('%03i.txt',slice_n)]);
        movefile([Enface_Fiji,'/',FijiInfo.FileBase_tilt,sprintf('%03i.registered.txt',slice_n)],...
            [StitchingFiji,'/',FijiInfo.FileBase_tilt,sprintf('%03i.registered.txt',slice_n)]);
    end
end
elapsed_time = toc(t_stitch);
fprintf('Elapsed time for stitching slices is %i minutes\n',round(elapsed_time/60));



%% Run this to find median Fiji Tile Coordinates
%%%% (stored as Experiment_Fiji.X_Mean & .Y_Mean)

saveflag = 1;
[ Experiment_Fiji ] = ReadFijiTileCoordinates_NB( ParameterFile, saveflag);
% Note: Had to comment out Method 2 saving in ReadFijiTileCoordinates
% because of an error (?) with the Brainstem I58 data
save(ParameterFile,'Experiment_Fiji','-append');
%% Create sliceidx
tmp.in              = [Parameters.SliceID];                                %sliceid from multiple runs
% tmp.out             = 1:length(tmp.in); 
tmp.out             = Parameters.SliceID; 
tmp.run             = [Parameters.SliceID*0+1];
sliceidx            = [tmp.in;tmp.out;tmp.run];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlab Mosaic Parameter
% load(ExpFiji, 'Experiment_Fiji');
% Parameters.Experiment = Experiment_Fiji;
% save(ParameterFile,'Parameters','-append');

Mosaic2D.sliceidx   = sliceidx;
Mosaic2D.imresize   = 0; 

Mosaic2D.InFileType = 'nifti';    
Mosaic2D.file_format= Enface.input_format;                                            % 'nifti' 'mat'
% Mosaic2D.indir      = {Enface.indir};   
Mosaic2D.indir      = {'/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData'};   
Mosaic2D.outdir     = StitchingFiji;
Mosaic2D.Exp        = ExpFiji;
    
save(ParameterFile,'Mosaic2D','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mosaic 2D slices 

Modality_2D = {'AIP','MIP','Retardance','Birefringence','Orientation'};%,'mus'};

for m = 1:length(Modality_2D)
    modality = Modality_2D{m};       
    Mosaic2D_Telesto( ParameterFile, modality ); % Method_1 : median coordinates in z | Method_2 : median stepsize and offset
end
%% StackNii
Stack.indir         = StitchingFiji;
Stack.outdir        = StackNii;
Stack.sliceid       = sliceidx(2,:);
Stack.Resolution    = [Scan.Resolution(1) Scan.Resolution(2) Scan.Thickness/1000];
save(ParameterFile,'Stack','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% StackNii

Modality_2D = {'AIP','MIP','Retardance','Orientation','Birefringence'};%,'mus'};

for m = 1:4
    modality = Modality_2D{m};
    do_stacking(ParameterFile, modality)
end

%delete(findall(0,'type','figure','tag','TMWWaitbar'));
%% Mosaic3D
Mosaic3D.indir       = {RawDataDir}; % {Scan.RawDataDir};% raw data directory for spectral change to {[ProcDir '/MAT']};
Mosaic3D.outdir      = [Scan.ProcessDir]; % where stitching slices will save
Mosaic3D.file_format = Processed3D.output_format;  
Mosaic3D.InFileType  = 'nifti'; % 'nifti' for saved on surface finding 'mat' for spectral raw
Mosaic3D.invert_Z    = false; 
Mosaic3D.sliceidx    = sliceidx; % give the slice index
Mosaic3D.Exp         = ExpFiji; % pulls the coordinates from ExpFiji.mat
Mosaic3D.sxyz        = [1 1 4]; % [4 4 8]; % [1 1 4]; % Downsample by [x,y,z] = [2,2,8], 10/10/2.5 um --> 20/20/20 um (isotropic)

Mosaic3D.MZL         = 500; % how much of the volume to save in z in pixels
Mosaic3D.z_parameters = [0,60,0]; % Skirt for z-stitching (the middle number is thickness of slice in pixels)


save(ParameterFile,'Mosaic3D','-append');

%% do 3D stitching   
t = tic;
modality_3D = {'dBI','R3D','O3D','mus'};
for m = 1%:length(Modality_3D)
    Mosaic3D_Telesto(ParameterFile,modality_3D{m});
end

elapsed_time = toc(t);
fprintf('Elapsed time for stitching 1 slice is %i minutes\n',round(elapsed_time/60));
%% do Z Stitching
% do_Z_stitching_DG(ParameterFile);


%% Mus
% SaveMus.depth = 100;
% SaveMus.z_mean= 5:80;
% SaveMus.save %2d tile 3d tile%
% SaveMus.Method %Vermeer % Optical_fitting_MGH 

%% FijiParameters
% channels_for_registration=[Red, Green and Blue] 
% rgb_order=rgb 
% fusion_method=[Linear Blending] 
% fusion=1.50 
% regression=0.30 
% max/avg=2.50 
% absolute=3.50

%% 

% rawdatadir
%
% ['test_processed_' sprintf('%03i',i) '_' lower(modality) '.nii']
% sprintf('test_processed_%03i_%s.nii',i,lower(modality));
% test_processed_001_aip.nii
% test_processed_001_mip.nii
% test_processed_001_orientation.nii
% test_processed_001_retardance.nii

% test_processed_001_cropped.nii
% test_processed_001_surface_finding.nii


% Enface_Tiff
%
% AIP_001.tiff
% MIP_001.tiff
% Retardance_001.tiff
% Orientation1_001.tiff
% Orientation2_001.tiff


% StitchingFiji
% 
% AIP_slice001.mat
% MIP_slice001.mat
% Retardance_slice001.mat
% Orientation_slice001.mat
% mus_slice001.mat
% 
% AIP_slice001.tiff
% AIP_Slice001.tif                          <- fiji output
% MIP_slice001.tiff
% Retardance_slice001.tiff
% Orientation1_slice001.tiff
% Orientation2_slice001.tiff
% 
% resized-AIP_slice001.jpg
% resized-MIP_slice001.jpg
% resized-Retardance_slice001.jpg
% resized,nearest-Orientation1_slice001.tiff
% 
% Mosaic_depth001_slice001.txt              <- WriteMacro output
% Mosaic_depth001_slice001.txt.registered   <- fiji output


% StackNii
% 
% Stacked_AIP.nii.gz 
% Stacked_MIP.nii.gz  
% Stacked_Retardance.nii.gz 
% Stacked_mus.nii.gz


