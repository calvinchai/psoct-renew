function S2P_Spectral_to_Processed(ParameterFile,istart)
load(ParameterFile);

% Input_Path = Processed3D.indir; % Scan.RawDataDir;
Input_Path = Scan.RawDataDir;
Input_File_Format = Scan.FileNameFormat;
Output_Path_3D = Processed3D.outdir;
Output_Files_3D = Processed3D.save;
Output_Format_3D = Processed3D.output_format;

Output_Path_2D = Enface.outdir;
Output_Files_2D = Enface.save;
Output_Format_2D = Enface.file_format;

if ~exist(Output_Path_3D,"dir"); mkdir(Output_Path_3D);end
if ~exist(Output_Path_2D,"dir"); mkdir(Output_Path_2D);end
dispCompFile = Processed3D.dispComp;

spectral_filenames = dir([Input_Path,'/*.nii']);

flag_mus = false;
if flag_mus; fprintf('ATTENTION: mus will be calculated.\n If current machine does not have GPU, it will be very slow.\n');end

% fname_format = 'mosaic_%03i_image_%03i_x_%06i_y_%06i_z_%06i_r_%06i_normal_spectral_%04i.nii';
ind_format = regexp(Input_File_Format, 'mosaic_%0(\d+)i','tokens');
m_num_length = str2double(ind_format{1,1});
ind_format = regexp(Input_File_Format, 'image_%0(\d+)i','tokens');
t_num_length = str2double(ind_format{1,1});
ind_format = regexp(Input_File_Format, '_%0(\d+)i.nii','tokens');
ft_num_length = str2double(ind_format{1,1});

for i = 1:length(spectral_filenames)
    fname = spectral_filenames(i,1).name;
    fpath = spectral_filenames(i,1).folder;
    spectral_filepath = [fpath,'/',fname];

    format_str = 'mosaic_';
    format_ind = strfind(fname,format_str);
    mosaic_num = str2double(fname(format_ind + length(format_str):format_ind + length(format_str) + m_num_length - 1));
    if ~ismember(mosaic_num, Parameters.MosaicID); continue; end;
    format_str = 'image_';
    format_ind = strfind(fname,format_str);
    tile_num = str2double(fname(format_ind + length(format_str):format_ind + length(format_str) + t_num_length - 1));
    % if mosaic_num==19 &tile_num<9; continue; end;
    format_str = '.nii';
    format_ind = strfind(fname,format_str);
    full_tile_num = str2double(fname(format_ind - ft_num_length + 1:format_ind - 1)); % Tile num across whole experiment

    output_fname_format_3D = sprintf(Output_Format_3D,mosaic_num,tile_num);
    output_fname_format_2D = sprintf(Output_Format_2D,mosaic_num,tile_num);

    % Spectral to Processed (dBI3D, R3D, O3D)
    tile_info = niftiinfo(spectral_filepath);
    [dBI3D, R3D, O3D] = Save3D_tile_20250429(spectral_filepath,tile_num,dispCompFile,tile_info.ImageSize(1,2),tile_info.ImageSize(1,3));
    dBI3D = single(dBI3D);
    R3D = single(R3D);
    O3D = single(O3D);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write uncropped volumes
    % Write dBI3D.mat
    idx = strcmpi(Output_Files_3D,'dBI3D');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,dBI3D);
    end

    % Write R3D.mat
    idx = strcmpi(Output_Files_3D,'R3D');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,R3D);
    end

    % Write O3D.mat
    idx = strcmpi(Output_Files_3D,'O3D');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,O3D);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find tissue surface and crop 3D volumes
    refine_range = 220;
    surf = SurfaceFinding_3D(dBI3D,istart,refine_range);
    surf = surf + 5;

    % depth = Scan.CropDepth; % Number of pixels below the surface to save in z
    depth = 120;
    sz = size(dBI3D);

    % [X,Y] = meshgrid(1:sz(1),1:sz(2));
    [X,Y] = meshgrid(1:sz(2),1:sz(1));
    lnridx = sub2ind(sz,Y(:),X(:),surf(:));
    lnridx = reshape(lnridx+([0:depth-1]*prod(sz(1:2))),sz(1),sz(2),depth);

    dBI3D_crop = dBI3D(lnridx);
    R3D_crop = R3D(lnridx);
    O3D_crop = O3D(lnridx);

    % Write dBI3D_crop.mat
    idx = strcmpi(Output_Files_3D,'dBI3D_crop');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,dBI3D_crop)
    end
    % Write R3D_crop.mat
    idx = strcmpi(Output_Files_3D,'R3D_crop');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,R3D_crop)
    end
    % Write O3D_crop.mat
    idx = strcmpi(Output_Files_3D,'O3D_crop');
    if(~isempty(find(idx,1)))
        output_fname = replace(output_fname_format_3D,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_3D,output_fname);
        parsave(fpath,O3D_crop)
    end

    % Write surf.mat
    idx = strcmpi(Output_Files_2D,'surf');
    if(~isempty(find(idx,1)))
        fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        parsave(fpath,surf)
    end

    % Write aip.mat
    idx = strcmpi(Output_Files_2D,'aip');
    if(~isempty(find(idx,1)))
        aip = mean(dBI3D_crop,3,"omitnan");
        fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        parsave(fpath,aip)
    end

    % Write mip.mat
    idx = strcmpi(Output_Files_2D,'mip');
    if(~isempty(find(idx,1)))
        mip = max(dBI3D_crop,[],3,"omitnan");
        fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        parsave(fpath,mip)
    end

    % Write ret.mat
    idx = strcmpi(Output_Files_2D,'ret');
    if(~isempty(find(idx,1)))
        ret = mean(R3D_crop,3,"omitnan");
        fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        parsave(fpath,ret)
    end

    % Write ori.mat
    idx = strcmpi(Output_Files_2D,'ori');
    if(~isempty(find(idx,1)))
        ori = orien_enface(O3D_crop,5);
        fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        parsave(fpath,ori)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit retardance for birefringence and write biref.mat
    idx = strcmpi(Output_Files_2D,'biref');
    if(~isempty(find(idx,1)))
        if(strcmpi(Scan.System,'Telesto'))
            wavelength = .0013; % 1300 nm
        else
            disp('Need wavelength to calculate birefringence (degree/um)');
            return;
        end

        depth_lin_ret = 119; % Depth over which retardance is approximately linear
        biref = zeros(size(R3D,1),size(R3D,2));
        % for ii = 1:size(R3D,1)
        %     for jj = 1:size(R3D,2)
        %         % if(surf(ii,jj) + depth_lin_ret > 500) % If it exceeds crop depth of 500 pixels, only go until end
        %         %     ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):size(R3D,3)))/360*wavelength;
        %         %     depth = px_size_z*[surf(ii,jj):size(R3D,3)];
        %         % else
        %         ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):surf(ii,jj) + depth_lin_ret))/360*wavelength;
        %         depth = Scan.Resolution(1,3)*(surf(ii,jj):surf(ii,jj) + depth_lin_ret);
        %         % end
        %         [p,~] = polyfit(depth,ret_profile',1);
        %         biref(ii,jj) = p(1); % Units: degree/um
        %     end
        % end

        for ii = 1:size(R3D_crop,1)
            for jj = 1:size(R3D_crop,2)
                ret_profile = squeeze(R3D_crop(ii,jj,1:1 + depth_lin_ret))/360*wavelength;
                depth = Scan.Resolution(1,3)*(1:1 + depth_lin_ret);
                [p,~] = polyfit(depth,ret_profile',1);
                biref(ii,jj) = p(1); % Units: degree/um
            end
        end
    end
    fname = replace(output_fname_format_2D,'[modality]',Output_Files_2D{idx});
    fpath = sprintf('%s/%s',Output_Path_2D,fname);
    parsave(fpath,biref)
end

% % mus
% if flag_mus ==1
%     depth = 600;
%     lnridx = sub2ind(sz,Y(:),X(:),surf(:));
%     lnridx = reshape(lnridx+([0:depth-1]*prod(sz(1:2))),sz(1),sz(2),depth);
%
%     dBI_crop = dBI3D(lnridx);
%     mus = dBI2mus2D(dBI_crop); % faster with gpu machine
%
%     fname = sprintf('%s/%s_%03i.mat',Output_Path_2D,"mus",FileNum);
%     parsave(fname,mus)
% end
end



% function I = dBI2mus2D(I3)
%
%
% if ~isa(I3, 'single'); I3 = single(I3);end
% if canUseGPU(); I3 =  gpuArray(I3);end
% % if strcmpi(Scan.System, 'Octopus'); I3 = permute(I3,[2,1,3]);end
% % if transpose_flag;                  I3 = permute(I3,[2,1,3]);end
%
% % I3 = I3(XPixClip+1:end,YPixClip+1:end,:);
% sz = size(I3);
% sz(3) = 120; % average depth
% % I3 = I3(:,:,end:-1:1); % uncomment if z direction is flipped
% I3 = 10.^(I3/10);
%
% data =  zeros(sz,'like',I3);
% for z=1:sz(3)-1
%     data(:,:,z) = I3(:,:,z)./(sum(I3(:,:,z+1:end),3))/2/0.0025;
% end
% if canUseGPU(); data =  gather(data);end
% I = squeeze(mean(data,3)); %figure;plot(squeeze(mean(data,[1 2])))
%
% end
