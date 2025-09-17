function Processed3D_to_Enface(ParameterFile,istart,is_cropped)

load(ParameterFile);

Input_Path = Scan.RawDataDir;
% Input_Path = Processed3D.outdir;
Output_Path = Enface.indir;
Output_File_Format = Enface.input_format;
Output_Files = Enface.save;

if ~exist(Output_Path,"dir"); mkdir(Output_Path);end

mat_files = find(contains(Scan.FileNameFormat,'.mat'));
nii_files = find(contains(Scan.FileNameFormat,'.nii'));
file_prefix = cell(size(Scan.FileNameFormat));
for n = 1:size(Scan.FileNameFormat,1)
    file_prefix{n,1} = Scan.FileNameFormat{n,1}(1:strfind(Scan.FileNameFormat{n,1},'_mosaic_%') - 1);
    % file_prefix{n,1} = Scan.FileNameFormat{n,1}(1:strfind(Scan.FileNameFormat{n,1},'_%') - 1);
end
dBI3D_idx = find(contains(file_prefix,'dBI') == 1);
R3D_idx = find(contains(file_prefix,'R3D') == 1);
O3D_idx = find(contains(file_prefix,'O3D') == 1);

idx = strcmpi(Output_Files,'mus');
if(~isempty(find(idx,1)))
    fprintf('ATTENTION: mus will be calculated.\nIf current machine does not have GPU, it will be very slow.\n');
end

mosaic_nums = Parameters.MosaicID;
tile_nums = Parameters.TileID;
for i = 1:length(mosaic_nums)
    fprintf('Mosaic #%i, Tile #%i\n',mosaic_nums(i),tile_nums(i));

    % Load Processed3D Data
    for n = 1:size(Scan.FileNameFormat,1)
        filename = sprintf(Scan.FileNameFormat{n,1},mosaic_nums(i),tile_nums(i));
        % filename = sprintf(Scan.FileNameFormat{n,1},tile_nums(i));
        switch n
            case dBI3D_idx
                if(~isempty(mat_files))
                    dBI3D = load(sprintf('%s/%s',Input_Path,filename)).(file_prefix{dBI3D_idx});
                elseif(~isempty(nii_files))
                    dBI3D = niftiread(sprintf('%s/%s',Input_Path,filename));
                end
            case R3D_idx
                if(~isempty(mat_files))
                    R3D = load(sprintf('%s/%s',Input_Path,filename)).(file_prefix{R3D_idx});
                elseif(~isempty(nii_files))
                    R3D = niftiread(sprintf('%s/%s',Input_Path,filename));
                end
            case O3D_idx
                if(~isempty(mat_files))
                    O3D = load(sprintf('%s/%s',Input_Path,filename)).(file_prefix{O3D_idx});
                elseif(~isempty(nii_files))
                    O3D = niftiread(sprintf('%s/%s',Input_Path,filename));
                end    
        end
    end

    % Output file name format
    fname_format = sprintf(Output_File_Format,mosaic_nums(i),tile_nums(i));
    % fname_format = sprintf(Output_File_Format,tile_nums(i));

    if(is_cropped == false) % For uncropped 3D data, need to do surface finding
        refine_range = 220;
        surf = SurfaceFinding_3D(dBI3D,istart,refine_range);
        surf = surf + 5;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% Added to make cropping work for test %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        surf(surf > 500) = 500;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        depth = Scan.CropDepth; % Number of pixels below the surface to save in z
        sz = size(dBI3D);

        % [X,Y] = meshgrid(1:sz(1),1:sz(2));
        [X,Y] = meshgrid(1:sz(2),1:sz(1));
        lnridx = sub2ind(sz,Y(:),X(:),surf(:));
        lnridx = reshape(lnridx+([0:depth-1]*prod(sz(1:2))),sz(1),sz(2),depth);

        dBI3D_crop = dBI3D(lnridx);
        R3D_crop = R3D(lnridx);
        O3D_crop = O3D(lnridx);

        idx = strcmpi(Output_Files,'surf');
        if(~isempty(find(idx,1)))
            fname = replace(fname_format,'[modality]',Output_Files{idx});
            fpath = sprintf('%s/%s',Output_Path,fname);
            parsave(fpath,surf)
        end

        save_cropped_3D = false;
        if(save_cropped_3D == true)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%%%%%%%%%%%%%%%%%%%%%% Added to simulate having cropped .nii files to start %%%%%%%%%%%%%%%%%%%%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filename = sprintf('dBI_crop_mosaic_%03i_image_%04i.nii',mosaic_nums(i),tile_nums(i));
            % niftiwrite(dBI_crop,['/autofs/cluster/connects2/users/Nate/test_Case_158_Brainstem_1_31_25_Mosaic_1/nifti3D/',filename]);
            % filename = sprintf('R3D_crop_mosaic_%03i_image_%04i.nii',mosaic_nums(i),tile_nums(i));
            % niftiwrite(R3D_crop,['/autofs/cluster/connects2/users/Nate/test_Case_158_Brainstem_1_31_25_Mosaic_1/nifti3D/',filename]);
            % filename = sprintf('O3D_crop_mosaic_%03i_image_%04i.nii',mosaic_nums(i),tile_nums(i));
            % niftiwrite(O3D_crop,['/autofs/cluster/connects2/users/Nate/test_Case_158_Brainstem_1_31_25_Mosaic_1/nifti3D/',filename]);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Write dBI3D_crop
            fname = replace(fname_format,'[modality]','dBI_crop');
            fpath = sprintf('%s/%s',Input_Path,fname); % Save to 3D input folder
            parsave(fpath,dBI3D_crop)

            % Write R3D_crop
            fname = replace(fname_format,'[modality]','dBI_crop');
            fpath = sprintf('%s/%s',Input_Path,fname); % Save to 3D input folder
            parsave(fpath,dBI3D_crop)

            % Write O3D_crop
            fname = replace(fname_format,'[modality]','dBI_crop');
            fpath = sprintf('%s/%s',Input_Path,fname); % Save to 3D input folder
            parsave(fpath,dBI3D_crop)
        end

    else % For cropped 3D data, directly determine Enface maps
        dBI3D_crop = dBI3D;
        R3D_crop = R3D;
        O3D_crop = O3D;

    end

    % Write aip.mat
    idx = strcmpi(Output_Files,'aip');
    if(~isempty(find(idx,1)))
        aip = mean(dBI3D_crop,3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files{idx});
        fpath = sprintf('%s/%s',Output_Path,fname);
        parsave(fpath,aip)
    end

    % Write mip.mat
    idx = strcmpi(Output_Files,'mip');
    if(~isempty(find(idx,1)))
        mip = max(dBI3D_crop,[],3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files{idx});
        fpath = sprintf('%s/%s',Output_Path,fname);
        parsave(fpath,mip)
    end

    % Write ret.mat
    idx = strcmpi(Output_Files,'ret');
    if(~isempty(find(idx,1)))
        ret = mean(R3D_crop,3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files{idx});
        fpath = sprintf('%s/%s',Output_Path,fname);
        parsave(fpath,ret)
    end

    % Write ori.mat
    idx = strcmpi(Output_Files,'ori');
    if(~isempty(find(idx,1)))
        ori = orien_enface(O3D_crop,5);
        fname = replace(fname_format,'[modality]',Output_Files{idx});
        fpath = sprintf('%s/%s',Output_Path,fname);
        parsave(fpath,ori)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit retardance for birefringence and write biref.mat
    idx = strcmpi(Output_Files,'biref');
    if(~isempty(find(idx,1)))
        if(strcmpi(Scan.System,'Telesto'))
            wavelength = .0013; % 1300 nm
        else
            disp('Need wavelength to calculate birefringence (degree/um)');
            return;
        end

        depth_lin_ret = 200; % Depth over which retardance is approximately linear
        biref = zeros(size(R3D,1),size(R3D,2));
        if(is_cropped == false)
            for ii = 1:size(R3D,1)
                for jj = 1:size(R3D,2)
                    % if(surf(ii,jj) + depth_lin_ret > 500) % If it exceeds crop depth of 500 pixels, only go until end
                    %     ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):size(R3D,3)))/360*wavelength;
                    %     depth = px_size_z*[surf(ii,jj):size(R3D,3)];
                    % else
                    ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):surf(ii,jj) + depth_lin_ret))/360*wavelength;
                    depth = Scan.Resolution(1,3)*(surf(ii,jj):surf(ii,jj) + depth_lin_ret);
                    % end
                    [p,~] = polyfit(depth,ret_profile',1);
                    biref(ii,jj) = p(1); % Units: degree/um
                end
            end

        else
            for ii = 1:size(R3D,1)
                for jj = 1:size(R3D,2)
                    ret_profile = squeeze(R3D(ii,jj,1:1 + depth_lin_ret))/360*wavelength;
                    depth = Scan.Resolution(1,3)*(1:1 + depth_lin_ret);
                    [p,~] = polyfit(depth,ret_profile',1);
                    biref(ii,jj) = p(1); % Units: degree/um
                end
            end
        end
        fname = replace(fname_format,'[modality]',Output_Files{idx});
        fpath = sprintf('%s/%s',Output_Path,fname);
        parsave(fpath,biref)
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine mus and write mus.mat
idx = strcmpi(Output_Files,'mus');
if(~isempty(find(idx,1)))

    disp('mus code not up to date'); % Focus for tilted illumination?
    return;
    % Get Focus and Rolloff data
    if(strcmpi(Scan.System,'Telesto'))
        % Load in previously generated data for Focus?
        path_focus = '/autofs/space/darwin_002/users/Agarose_Focus_tot.mat'; %add this parameter in Parameters.mat
        
        path_rolloff = '/autofs/cluster/octdata2/users/TelestoCalibration/RollOff_Spectrum.mat';
        RollOff = load(path_rolloff).('RollOff');

        RollOff = repmat(RollOff,[1,size(dBI3D,1),size(dBI3D,1)]);
        RollOff = RollOff(25:end,:,:);
        RollOff_reshape = reshape(RollOff(Mask),[size(tab,2) 700 700]);
        RollOff_reshape = permute(RollOff_reshape,[2 3 1]);
    else
        disp('Need focus and rolloff data to determine mus');
        return;
    end

    I = 10.^(dBI3D/10).*RollOff_reshape;
    newdBI3D_reshape = 10 .* log10(I);

    I = smooth3(I,'gaussian');
    runningsum = sum(I(:,:,1:end),3);
    for zz = 1:500
        runningsum = runningsum - sum(I(:,:,zz),3);
        mus(:,:,zz) = I(:,:,zz)./runningsum/2/Scan.Resolution(1,3);
    end
    fname = replace(fname_format,'[modality]',Output_Files{idx});
    fpath = sprintf('%s/%s',Output_Path,fname);
    parsave(fpath,mus)
end



end




