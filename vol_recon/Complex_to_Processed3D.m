function Complex_to_Processed3D(ParameterFile)

load(ParameterFile);

if isempty(gcp('nocreate'))
    parpool(8);
end

Input_Path = Scan.RawDataDir;
Input_File_Format = Scan.FileNameFormat;
Input_Files = Scan.FilePrefix;
Output_Path_3D = Processed3D.outdir;
Output_Files_3D = Processed3D.save;
Output_Format_3D = Processed3D.output_format;

Output_Path_2D = Enface.indir;
Output_Files_2D = Enface.save;
% Output_Format_2D = Enface.input_format;
if(strcmpi(Scan.CropMethod,'focus') == 1)
    surf_idx = find(contains(Output_Files_2D,'surf') == 1);
    focus_normal = niftiread(Scan.FocusFile);
    if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
        focus_tilt = niftiread(Scan.FocusFile_tilt);
    end
end

if ~exist(Output_Path_3D,"dir"); mkdir(Output_Path_3D);end
if ~exist(Output_Path_2D,"dir"); mkdir(Output_Path_2D);end
dispCompFile = Processed3D.dispComp;

mat_files = find(contains(Input_File_Format,'.mat'),1);
nii_files = find(contains(Input_File_Format,'.nii'),1);

complex_idx = find(contains(Input_Files,'cropped') == 1);

idx = strcmpi(Output_Files_3D,'mus');
if(~isempty(find(idx,1)))
    fprintf('ATTENTION: mus will be calculated.\nIf current machine does not have GPU, it will be very slow.\n');
end

mosaic_nums = Parameters.MosaicID;
tile_nums = Parameters.TileID;
for i = 1:length(mosaic_nums)
% for i = 980
    fprintf('Mosaic #%i, Tile #%i\n',mosaic_nums(i),tile_nums(i));

    if (tile_nums(i)<1000) 
        continue; 
    end
    if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
        if(mod(mosaic_nums(i),2) == 1) % Odd mosaic nums --> Normal incidence
            num_x_px = Scan.NbPixels;
            determine_mus = true;
            if(strcmpi(Scan.CropMethod,'focus') == 1)
                focus = focus_normal;
            end
        elseif(mod(mosaic_nums(i),2) == 0) % Even mosaic nums --> Tilted incidence
            num_x_px = Scan.NbPixels_tilt;
            determine_mus = false;
            if(strcmpi(Scan.CropMethod,'focus') == 1)
                focus = focus_tilt;
            end
        end
    else
        num_x_px = Scan.NbPixels;
        determine_mus = true;
    end
    if(strcmpi(Scan.CropMethod,'focus') == 1)
        % Load surface for focus-cropped 3D volumes
        surf_filename = replace(sprintf(Input_File_Format,mosaic_nums(i),tile_nums(i)),'[modality]',Output_Files_2D{1,surf_idx});
        surf = niftiread(sprintf('%s/%s',Output_Path_2D,surf_filename));
        surf = round((1024 - surf) - focus + Scan.Focus_CropStart);
        surf = surf + 15;
        surf(surf <= 0) = 1;
    end

    % Load Complex Processed Data
    for n = 1:length(Input_Files)
        filename = replace(sprintf(Input_File_Format,mosaic_nums(i),tile_nums(i)),'[modality]',Input_Files{1,n});

        switch n
            case complex_idx
                if(~isempty(mat_files))
                    complex3D = load(sprintf('%s/%s',Input_Path,filename)).(Input_Files{complex_idx});
                elseif(~isempty(nii_files))
                    complex3D = niftiread(sprintf('%s/%s',Input_Path,filename));
                end
                Jones1_real = complex3D(0*num_x_px + 1:1*num_x_px,:,:);
                Jones1_imag = complex3D(1*num_x_px + 1:2*num_x_px,:,:);
                Jones2_real = complex3D(2*num_x_px + 1:3*num_x_px,:,:);
                Jones2_imag = complex3D(3*num_x_px + 1:4*num_x_px,:,:);
                Jones1 = Jones1_real + sqrt(-1)*Jones1_imag;
                Jones2 = Jones2_real + sqrt(-1)*Jones2_imag;
        end
    end

    % Output file name format
    fname_format = sprintf(Output_Format_3D,mosaic_nums(i),tile_nums(i));

    %%%% dBI3D %%%%
    IJones = abs(Jones1).^2+abs(Jones2).^2;
    dBI3D = flip(10*log10(IJones),3);

    %%%% R3D %%%%
    R3D = flip(atan(abs(Jones1)./abs(Jones2))/pi*180,3);

    %%%% O3D %%%%
    offset = 0; %100/180*pi; 
    phase1 = angle(Jones1);
    phase2 = angle(Jones2);
    phi = (phase1 - phase2) + offset*2;
    index1 = (phi > pi);
    phi(index1) = phi(index1) - 2*pi;
    index2 = (phi < -pi);
    phi(index2) = phi(index2) + 2*pi;
    O3D = flip(phi/2/pi*180,3);
    % 
    % if(mosaic_nums(i) == 7)
    %     for x = 201:350
    %         for y = 1:350
    %             Aline_dBI3D = dBI3D(y,x,150 + (round(focus(y,x)) - 150):end);
    %             dBI3D(y,x,:) = zeros(1,1,size(dBI3D,3));
    %             dBI3D(y,x,1:numel(Aline_dBI3D)) = Aline_dBI3D; 
    % 
    %             Aline_R3D = R3D(y,x,150 + (round(focus(y,x)) - 150):end);
    %             R3D(y,x,:) = zeros(1,1,size(R3D,3));
    %             R3D(y,x,1:numel(Aline_R3D)) = Aline_R3D; 
    % 
    %             Aline_O3D = O3D(y,x,150 + (round(focus(y,x)) - 150):end);
    %             O3D(y,x,:) = zeros(1,1,size(O3D,3));
    %             O3D(y,x,1:numel(Aline_O3D)) = Aline_O3D; 
    %         end
    %     end
    % end
    % if(mod(mosaic_nums(i),2) == 0 && mosaic_nums(i) > 8)
    %     [y_corr,x_corr] = find(focus <= Scan.Focus_CropStart); % Indices that need correction
    %     for n_corr = 1:length(x_corr)
    %         inf_idx = find(dBI3D(y_corr(n_corr),x_corr(n_corr),:) == -inf);
    % 
    %         Aline_dBI3D = dBI3D(y_corr(n_corr),x_corr(n_corr),max(inf_idx):end);
    %         dBI3D(y_corr(n_corr),x_corr(n_corr),:) = zeros(1,1,size(dBI3D,3));
    %         dBI3D(y_corr(n_corr),x_corr(n_corr),1:numel(Aline_dBI3D)) = Aline_dBI3D; 
    % 
    %         Aline_R3D = R3D(y_corr(n_corr),x_corr(n_corr),max(inf_idx):end);
    %         R3D(y_corr(n_corr),x_corr(n_corr),:) = zeros(1,1,size(R3D,3));
    %         R3D(y_corr(n_corr),x_corr(n_corr),1:numel(Aline_R3D)) = Aline_R3D; 
    % 
    %         Aline_O3D = O3D(y_corr(n_corr),x_corr(n_corr),max(inf_idx):end);
    %         O3D(y_corr(n_corr),x_corr(n_corr),:) = zeros(1,1,size(O3D,3));
    %         O3D(y_corr(n_corr),x_corr(n_corr),1:numel(Aline_O3D)) = Aline_O3D; 
    %     end
    % end

    % Write dBI3D
    idx = strcmpi(Output_Files_3D,'dBI3D');
    if(~isempty(find(idx,1)))
        fname = replace(fname_format,'[modality]','dBI3D');
        fpath = sprintf('%s/%s',Output_Path_3D,fname); % Save to 3D input folder
        if(~isempty(mat_files))
            parsave(fpath,dBI3D)
        elseif(~isempty(nii_files))
            niftiwrite(dBI3D,fpath);
        end
    end

    % Write R3D
    idx = strcmpi(Output_Files_3D,'R3D');
    if(~isempty(find(idx,1)))
        fname = replace(fname_format,'[modality]','R3D');
        fpath = sprintf('%s/%s',Output_Path_3D,fname); % Save to 3D input folder
        if(~isempty(mat_files))
            parsave(fpath,R3D);
        elseif(~isempty(nii_files))
            niftiwrite(R3D,fpath);
        end
    end

    % Write O3D
    idx = strcmpi(Output_Files_3D,'O3D');
    if(~isempty(find(idx,1)))
        fname = replace(fname_format,'[modality]','O3D');
        fpath = sprintf('%s/%s',Output_Path_3D,fname); % Save to 3D input folder
        if(~isempty(mat_files))
            parsave(fpath,O3D);
        elseif(~isempty(nii_files))
            niftiwrite(O3D,fpath);
        end
    end
    % 
    % % Write mus
    % idx = strcmpi(Output_Files_3D,'mus');
    % if(determine_mus == true)
    %     if(~isempty(find(idx,1)))
    %         fprintf('mus code not updated\n');
    %         return;
    %         Focus_pathname = '/autofs/space/darwin_002/users/Focus_tot.mat';%add this parameter in Parameters.mat
    %         load(Focus_pathname);
    %         Focus = Focus_tot;%Focus_median; Surf_tot;
    % 
    %         [X,Y] = meshgrid(1:Bline);
    %         tab = sub2ind([DepthL,Bline,Aline_dBI3D],round(Focus(:)),Y(:),X(:)) + (-150:349);
    %         Mask = false([DepthL,Bline,Aline_dBI3D]);
    %         Mask(tab(:))=true;
    % 
    %         tab_I = sub2ind([DepthL,Bline,Aline_dBI3D],round(Focus(:)),Y(:),X(:)) + (-150:549);
    %         Mask_I = false([DepthL,Bline,Aline_dBI3D]);
    %         Mask_I(tab_I(:))=true;
    % 
    %         load(Processed3D.RollOff);
    %         RollOff = repmat(RollOff,[1,Scan.NbPixels,Scan.NbPixels]);
    %         RollOff = RollOff(25:end,:,:);
    %         RollOff_reshape = reshape(RollOff(Mask_I),[size(tab_I,2) 700 700]);
    %         RollOff_reshape = permute(RollOff_reshape,[2 3 1]);
    % 
    %         I = 10.^(dBI3D/10).*RollOff_reshape;
    %         newdBI3D_reshape = 10.*log10(I);
    % 
    %         I = smooth3(I,'gaussian');
    %         runningsum = sum(I(:,:,1:end),3);
    %         for zz = 1:500
    %             runningsum = runningsum - I(:,:,zz);
    %             % mus(:,:,zz) = I(:,:,zz)./sum(I(:,:,zz+1:end),3)/2/z_pixsize;
    %             mus(:,:,zz) = I(:,:,zz)./runningsum/2/z_pixsize;
    %         end
    % 
    %         fname = replace(fname_format,'[modality]','mus');
    %         fpath = sprintf('%s/%s',Output_Path_3D,fname); % Save to 3D input folder
    %         if(~isempty(mat_files))
    %             parsave(fpath,mus);
    %         elseif(~isempty(nii_files))
    %             niftiwrite(mus,fpath);
    %         end
    %     end
    % end

    % Write aip.mat
    idx = strcmpi(Output_Files_2D,'aip');
    if(~isempty(find(idx,1)))
        aip = mean(dBI3D,3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        % parsave(fpath,aip);
        niftiwrite(aip,fpath)
    end

    % Write mip.mat
    idx = strcmpi(Output_Files_2D,'mip');
    if(~isempty(find(idx,1)))
        mip = max(dBI3D,[],3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        % parsave(fpath,mip);
        niftiwrite(mip,fpath)
    end

    % Write ret.mat
    idx = strcmpi(Output_Files_2D,'ret');
    if(~isempty(find(idx,1)))
        ret = mean(R3D,3,"omitnan");
        fname = replace(fname_format,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        % parsave(fpath,ret);
        niftiwrite(ret,fpath)
    end

    % Write ori.mat
    idx = strcmpi(Output_Files_2D,'orientation');
    if(~isempty(find(idx,1)))
        ori = orien_enface(O3D,5);
        fname = replace(fname_format,'[modality]',Output_Files_2D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        % parsave(fpath,ori);
        niftiwrite(ori,fpath)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fit retardance for birefringence and write biref.mat
    idx = strcmpi(Output_Files_3D,'biref');
    if(~isempty(find(idx,1)))
        if(strcmpi(Scan.System,'Telesto'))
            wavelength = .0013; % 1300 nm
        else
            disp('Need wavelength to calculate birefringence (degree/um)');
            return;
        end
        depth_lin_ret = 100
        % if(mosaic_nums(i) == 7)
        %     depth_lin_ret = 100;
        % else
        %     depth_lin_ret = 170; % Depth over which retardance is approximately linear
        % end
        
        % % Use sliding window to analyze peaks of retardance depth profiles
        % window = 5;
        % w_s = (window - 1)/2;
        % y_window = 3:window:size(R3D,1);
        % x_window = 3:window:size(R3D,2);
        % for y = y_window
        %     for x = x_window
        %         surf_window = surf((y - w_s):(y + w_s),(x - w_s):(x + w_s));
        %         surf_window_mean((y - w_s):(y + w_s),(x - w_s):(x + w_s)) = ones(size(surf_window))*mean2(surf_window);
        %         R3D_window = R3D((y - w_s):(y + w_s),(x - w_s):(x + w_s),:);
        %         R3D_window_mean((y - w_s):(y + w_s),(x - w_s):(x + w_s),:) = repmat(mean(mean(R3D_window,1),2),size(R3D_window,1),size(R3D_window,2));
        %     end
        % end

        % View retardance depth profiles
        if(1 == 0)

            figure(1);hold on;
            for n = 1:2
                switch n
                    case 1
                        % y_px = 43:52;
                        y_px = 43:47;
                        % x_px = 43:52;
                        x_px = 43:47;
                        color = 'red';
                        format_spec = 'r--';
                    case 2
                        y_px = size(R3D,1)-9:size(R3D,1);
                        x_px = size(R3D,2)-9:size(R3D,2);
                        color = 'blue';
                        format_spec = 'b--';
                end
                ret_profile = squeeze(mean(mean(R3D(y_px,x_px,:),1),2));
                avg_surf = mean2(surf(y_px,x_px));
                plot(ret_profile,color);
                xline(avg_surf,format_spec,'LineWidth',2);
                xline(avg_surf + depth_lin_ret,format_spec,'LineWidth',2);
            end

        end


        z_px_size = Scan.Resolution(1,3);
        biref = single(zeros(size(R3D,1),size(R3D,2)));
        num_row = size(R3D,1);
        num_col = size(R3D,2);
        if(strcmpi(Scan.CropMethod,'surface'))
            parfor ii = 1:num_row
                for jj = 1:num_col
                    ret_profile = squeeze(R3D(ii,jj,1:1 + depth_lin_ret))/360*wavelength;
                    depth = z_px_size*(1:1 + depth_lin_ret);
                    [p,~] = polyfit(depth,ret_profile',1);
                    biref(ii,jj) = p(1); % Units: degree/um
                end
            end
        else % CropMethod 'focus' or 'none'
            parfor ii = 1:num_row
            % for ii = 1:num_row
                for jj = 1:num_col
                    if(surf(ii,jj) + depth_lin_ret > size(R3D,3)) % If it exceeds crop depth of 500 pixels, only go until end
                        ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):size(R3D,3)))/360*wavelength;
                        depth = z_px_size*(surf(ii,jj):size(R3D,3));
                    else
                        ret_profile = squeeze(R3D(ii,jj,surf(ii,jj):surf(ii,jj) + depth_lin_ret))/360*wavelength;
                        depth = z_px_size*(surf(ii,jj):surf(ii,jj) + depth_lin_ret);
                        if(any(ret_profile == 0))
                            zero_idx = find(ret_profile == 0);
                            ret_profile = ret_profile(1:min(zero_idx));
                            depth = depth(1:min(zero_idx));
                        end
                    end
                    [p,~] = polyfit(depth,ret_profile',1);
                    biref(ii,jj) = p(1); % Units: degree/um
                end
            end
        end
        
        fname = replace(fname_format,'[modality]',Output_Files_3D{idx});
        fpath = sprintf('%s/%s',Output_Path_2D,fname);
        if(~isempty(mat_files))
            parsave(fpath,biref);
        elseif(~isempty(nii_files))
            niftiwrite(biref,fpath);
        end
    end

end
end




