clear
ParameterFile = 'D:/OCT_user/stitching_test_2025_04_14_TestI56_blockA_mineraloil_02282025/processed_data/Parameters.mat'; % <-------------
addpath('D:/OCT_user/stitching_test_2025_04_14_TestI56_blockA_mineraloil_02282025/vol_recon_proj');
slice_nums = 1;
depth_n = 1;
for slice_n = slice_nums
    %% Enface
    SaveTiff_Enface_script( ParameterFile )

    %% Create Macro.ijm file

    % WriteMacro_Spiral( ParameterFile );
    WriteMacro_Serpentine( ParameterFile );

    %% Run ImageJ (FIJI) without ssh for speed. e.g. tigervnc

    load(ParameterFile,'WriteMacro')
    fiji_path = '"C:\Program Files\Fiji.app\ImageJ-win64.exe"';
    macro_path = sprintf('"%s"', WriteMacro.SaveMacro);

    cmd = sprintf('%s --headless --console -macro %s', fiji_path, macro_path);
    system(cmd);
    
    movefile([Enface_Tiff,'/',sprintf('Mosaic_depth%03i_slice%03i.txt',depth_n,slice_n)],...
        [StitchingFiji,'/',sprintf('Mosaic_depth%03i_slice%03i.txt',depth_n,slice_n)]);
    movefile([Enface_Tiff,'/',sprintf('Mosaic_depth%03i_slice%03i.registered.txt',depth_n,slice_n)],...
        [StitchingFiji,'/',sprintf('Mosaic_depth%03i_slice%03i.registered.txt',depth_n,slice_n)]);

        %% Run this to find median Fiji Tile Coordinates
        %%%% (stored as Experiment_Fiji.X_Mean & .Y_Mean)

        saveflag = 1;
        [ Experiment_Fiji ] = ReadFijiTileCoordinates( ParameterFile, saveflag);
        save(ParameterFile,'Experiment_Fiji','-append');

        %% Mosaic 2D slices

        Modality_2D = {'AIP','MIP','Retardance','Orientation'};%,'mus'};%Modality_2D = {'mus_minIP1', 'mus_AIP'}

        for m = 1:4
        modality = Modality_2D{m};
        Mosaic2D_Telesto( ParameterFile, modality ); % Method_1 : median coordinates in z | Method_2 : median stepsize and offset
        end

        %% StackNii

        Modality_2D = {'AIP','MIP','Retardance','Orientation'};%,'mus'};

        for m = 1:4
            modality = Modality_2D{m};
            do_stacking(ParameterFile, modality)
        end

        %delete(findall(0,'type','figure','tag','TMWWaitbar'));

        %% Mosaic 3D slices

        modality = 'dBI';
        Mosaic3D_Telesto( ParameterFile, modality);

end