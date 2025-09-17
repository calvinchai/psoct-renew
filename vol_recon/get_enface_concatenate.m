function [c_enface,c_,c_enface_window] = get_enface_concatenate(Path,fname_format,target_str,mosaic_num,mosaic_dim,mosaicnum,transpose_flag,window)
% [c_enface, c_, info ] = get_enface_concatenate(Path,target_str,mosaic_dim,slicenum,transpose_flag)
% INPUT FORMAT:
% Path          the directory of tile enface file (mat/nii)
% target_str    unique string for your file of interest e.g. 'aip'
%               or a pattern e.g. "test_processed_" + digitsPattern + "_aip.nii"
% mosaic_dim    the dimension of the mosaic of the slice.
% slicenum      the # th slice that will be concatenated
% transpose_flag if true, tile enface will be transposed before concatenated
% window        the windowing parameter that determine dynamic range in c_enface_window

% OUTPUT FORMAT:
% c_enface      concatenated enface slice
% c_            cell array of the tile enfaces that created enface slice
% info          directary information filtered by target_str 
% c_enface_window windowed, tile# annotated, and concatenated enface slice


    if ~exist('transpose_flag','var')
        transpose_flag = false;
    end

    if ~exist('window','var')
        window_flag = false;
    else
        c__ = cell(1,prod(mosaic_dim));
        window_flag = true;
    end
    c_ = cell(1,prod(mosaic_dim));

    [mosaic_nums,tile_nums] = get_dir_info(Path,target_str);
    %% get the order
    order = get_serpent(mosaic_dim);
    targ_mosaic_nums = ones(1,prod(mosaic_dim)) + [mosaicnum-1]*ones(1,prod(mosaic_dim)); % mosaic number
    targ_tile_nums = 1:prod(mosaic_dim) + (mosaicnum-1)*prod(mosaic_dim); % tile number
    
    I = zeros(350,350);
    %% tile enface
    for i = 1:length(targ_mosaic_nums(:))
        
        idx_m = find(mosaic_nums == targ_mosaic_nums(i),1);
        idx_t = find(tile_nums == targ_tile_nums(i),1);
        
        if(isempty(idx_m) || isempty(idx_t)) 
            I = I*0;
            fprintf('missing tile in current slice at mosaic %i, tile %i\n',targ_mosaic_nums(i),targ_tile_nums(i));
        else
            filename = replace(sprintf(fname_format,targ_mosaic_nums(i),targ_tile_nums(i)),'[modality]',target_str);
            filepath = [Path '/' filename];
            I = load_input(filepath,transpose_flag);
        end
        
        c_{i} = I;

        if window_flag
            I_rgb = insertText(mat2gray(I,window),[20,20],targ_tiles(i),FontSize=30,BoxOpacity=0,TextColor="yellow");
            c__{i} = I_rgb;
        end
    end
    
    %% slice enface
    c_ = c_(order);
    c_enface = cell2mat(c_);
    if window_flag
        c__ = c__(order);
        c_enface_window = cell2mat(c__);
    else
        c_enface_window = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = load_input(VolName,transpose_flag)
    [~,~,ext] = fileparts(VolName);
    flag_check = 0; % check file content
    switch ext
        case '.nii'
            data = readnifti(VolName);
        case '.mat'
            tmp = load(VolName); 
            tmp = struct2cell(tmp); 
            data = tmp{1};
            flag_check = size(tmp,1)~=1;
    end

    if flag_check
        fprintf('Warning: input file (.mat) has multiple variables.\n')
        return
    end

    if transpose_flag
        % if artifact is not on the top side from imshow() or tile tiff; do transpose  
        data = data';
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mosaic_nums,tile_nums] = get_dir_info(Path,target_str)
    fprintf('the path for the experiment is %s\n',Path)
    
    directory_info = dir(Path);
    idx4niifile = contains({directory_info(:).name},target_str);
    info = directory_info(idx4niifile');

    numbers_in_filenames = str2double(extract([info(:).name],digitsPattern));
    mosaic_nums = numbers_in_filenames(1:2:end);
    tile_nums = numbers_in_filenames(2:2:end);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function order = get_serpent(mosaic_dim)

    tile_count = prod(mosaic_dim);

    order = reshape(1:tile_count,mosaic_dim);
    for odd=1:2:mosaic_dim(2);order(:,odd) = order(end:-1:1,odd);end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%