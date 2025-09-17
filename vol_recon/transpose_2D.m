%% Transpose images for compatibility with stitching pipeline
raw_data_dir = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData';
new_save_dir = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/nii2D_transpose';
file_format = 'mosaic_%03i_image_%04i_processed_[modality].nii';
% modality_str = {'aip','mip','retardance','orientation','surface_finding'}; % 'biref'};
modality_str = {'biref'};

TilesPerSlice = 588;
TilesPerSlice_tilt = 1008;
SliceID = [4:11]; % Mosaic 7 to 22
    
mosaic_nums = [];
tile_nums = [];

for n = 1:length(SliceID)
    tmp_mosaic_nums = [repelem(2*SliceID(n) - 1,TilesPerSlice),... % Normal
        repelem(2*SliceID(n),TilesPerSlice_tilt)]; % Tilted
    tmp_tile_nums = [1:TilesPerSlice,... % Normal
        1:TilesPerSlice_tilt]; % Tilted
    mosaic_nums = [mosaic_nums,tmp_mosaic_nums];
    tile_nums = [tile_nums,tmp_tile_nums];
end
%%

transpose_images(raw_data_dir,new_save_dir,mosaic_nums,tile_nums,file_format,modality_str);


%%

function transpose_images(indir,outdir,mosaic_n,tile_n,fname_format,modality)
    if(~exist(outdir,'dir'))
        mkdir(outdir);
    end

    for i = 1:length(mosaic_n)
        for n = 1:numel(modality)
            fprintf('Transposing mosaic_%03i_image_%04i_processed_%s.nii\n',mosaic_n(i),tile_n(i),modality{n});
            
            filename = sprintf(replace(fname_format,'[modality]',modality{n}),mosaic_n(i),tile_n(i));
            fpath = sprintf('%s/%s',indir,filename); 
            output_fpath = sprintf('%s/%s',outdir,filename);
            % fprintf('filename: %s\n',filename);
            % fprintf('input path: %s\n',fpath);
            % fprintf('output path: %s\n',output_fpath);
            img = niftiread(fpath);
            transposed_img = img';
            niftiwrite(transposed_img,output_fpath);
        end
    end
end