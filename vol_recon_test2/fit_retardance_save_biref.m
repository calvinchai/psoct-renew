function fit_retardance_save_biref(ParameterFile)
load(ParameterFile);

if isempty(gcp('nocreate'))
    parpool(8);
end

Input_File_Format = Scan.FileNameFormat;
Output_Path_3D = Processed3D.outdir;
Output_Format_3D = Processed3D.output_format;

Output_Path_2D = Enface.indir;
Output_Files_2D = Enface.save;
if(strcmpi(Scan.CropMethod,'focus') == 1)
    surf_idx = find(contains(Output_Files_2D,'surf') == 1);
    focus_normal = niftiread(Scan.FocusFile);
    if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
        focus_tilt = niftiread(Scan.FocusFile_tilt);
    end
end

mosaic_nums = Parameters.MosaicID;
tile_nums = Parameters.TileID;
% for i = 1:length(mosaic_nums)
for i = 17:300
    fprintf('Mosaic #%i, Tile #%i\n',mosaic_nums(i),tile_nums(i));

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
        surf(surf <= 0) = 1;
    end

    filename = replace(Output_Format_3D,'[modality]','R3D');
    filepath = sprintf('%s/%s',Output_Path_3D,sprintf(filename,mosaic_nums(i),tile_nums(i)));
    R3D = niftiread(filepath);
    


    if(strcmpi(Scan.System,'Telesto'))
        wavelength = .0013; % 1300 nm
    else
        disp('Need wavelength to calculate birefringence (degree/um)');
        return;
    end

    if(mosaic_nums(i) == 7)
        depth_lin_ret = 100;
    else
        depth_lin_ret = 170; % Depth over which retardance is approximately linear
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
    else % CropMethod focus or none
        parfor ii = 1:num_row
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

    fname = sprintf(replace(Output_Format_3D,'[modality]','biref'),mosaic_nums(i),tile_nums(i));
    fpath = sprintf('%s/%s',Output_Path_2D,fname);
    niftiwrite(biref,fpath);



end
end