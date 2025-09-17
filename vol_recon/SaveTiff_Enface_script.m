function SaveTiff_Enface_script( ParameterFile, nWorker )
% 
% SaveTiff_Enface_script( ParameterFile )
%   Save .nii or .mat tiles as .tiff files for specified 2D/en-face contrasts
% 
% USAGE:
% 
%  INPUTS:          
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters and Enface: 
% 
%   ~2024-03-06~  
% 


if isdeployed 
    disp('--App is deployed--');
else
%     disp('  --Local--  ');
    addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils/OCTBasic'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(ParameterFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Hard-coded vars
doSaveTif = false;

% display numbers of workers 
% fprintf('--nWorker is %i--\n', nWorker);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MODALITY TO SAVE 
modality        = Enface.save_tiff;

saveflag_aip = strcmpi(modality, 'aip');
saveflag_mip = strcmpi(modality, 'mip');
saveflag_ret = strcmpi(modality, 'Retardance') | strcmpi(modality, 'ret');
saveflag_ori = strcmpi(modality, 'Orientation') | strcmpi(modality, 'ori');
saveflag_biref = strcmpi(modality, 'biref');
saveflag_mus = strcmpi(modality, 'mus');
saveflag_surf = strcmpi(modality, 'surf');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING PARAMETERS 
XPixClip        = Parameters.XPixClip;
YPixClip        = Parameters.YPixClip;
transpose_flag  = Parameters.transpose;
% parameters for tiff (255) image intensity scale
if any(saveflag_aip);    aipGrayRange    = Parameters.AipGrayRange;end
if any(saveflag_mip);    mipGrayRange    = Parameters.MipGrayRange;end
if any(saveflag_ret);    retGrayRange    = Parameters.RetGrayRange;end
if any(saveflag_biref);    birefGrayRange    = Parameters.BirefGrayRange;end
if any(saveflag_mus);    musGrayRange    = Parameters.musGrayRange;end
if any(saveflag_surf);   surfGrayRange   = Parameters.Surfacerange;end
if any(saveflag_ori);    orientationSign = Parameters.OrientationSign;
                    orientationOffset = Parameters.OrienOffset;  end

% parameters to identify agar tile 
threshold_pct   = Parameters.Agar.threshold_pct;
threshold_std   = Parameters.Agar.threshold_std;
tissue_thresh   = Parameters.Agar.gm_aip;

% Vectors of tile and mosaic IDs to organize saving
tileid = Parameters.TileID;
mosaicid = Parameters.MosaicID;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir               = Enface.indir;
input_file_format   = [indir '/' Enface.input_format];
outdir              = Enface.outdir;
output_file_format  = [outdir '/' Enface.output_format];
agardir             = Parameters.Agar.dir;
agar_file_format    = [agardir '/' Parameters.Agar.file_format];

if(~exist(outdir,'dir') && doSaveTif)   
    mkdir(outdir);    
end
if(~exist(agardir,'dir'))
    mkdir(agardir);   
end


fprintf(' - Input directory = \n%s\n',indir);
fprintf(' - Output directory = \n%s\n',outdir);
fprintf(' - Agar directory = \n%s\n',agardir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poolobj = parpool(nWorker);
for i = 1:length(tileid)   
    fprintf('mosaic %d, tile %d\n',mosaicid(i),tileid(i));
    
    %AIP tile
    if any(saveflag_aip)
        VolName1 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_aip});
        if isfile(VolName1) 
            data1 = load_input(VolName1,transpose_flag);
            data1 = (data1(1+XPixClip:end,1+YPixClip:end));
            
            SLO_dBI = (data1);

            tissue_pct = nnz(SLO_dBI>tissue_thresh)/numel(SLO_dBI);
            agar_status = tissue_pct<threshold_pct & std(SLO_dBI,[],'all','omitnan')< threshold_std;
            tileinfo = [mosaicid(i),tileid(i),~agar_status];
            agar_filename = sprintf(agar_file_format,mosaicid(i),tileid(i));
            parsave(agar_filename,tileinfo);

        %     figure,imagesc(data1,AipGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_aip});
                imwrite(mat2gray(data1,aipGrayRange),...
                    output_filepath,...
                    'compression','none');
            end
        else
            fprintf('Missing %s\n',VolName1);
        end
    end
    
    %MIP tile
    if any(saveflag_mip)
        VolName2 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_mip});
        if isfile(VolName2) 
            data2 = load_input(VolName2,transpose_flag);
            data2 = (data2(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data2,MipGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_mip});
                imwrite(mat2gray(data2,mipGrayRange),...
                    output_filepath,...
                    'compression','none');  
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName2);
            end
        end
    end
    
    %Retardance tile
    if any(saveflag_ret)
        VolName3 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_ret});
        if isfile(VolName3) 
            data3 = load_input(VolName3,transpose_flag);
            data3 = (data3(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data3,RetGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_ret});
                imwrite(mat2gray(data3,retGrayRange),...
                    output_filepath,...
                    'compression','none');
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName3);
            end
        end  
    end

    % Orientation tile
    if any(saveflag_ori)
        VolName4 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_ori});
        if isfile(VolName4) 
            data4 = load_input(VolName4,transpose_flag);
            data4 = (data4(1+XPixClip:end,1+YPixClip:end));
            
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',[modality{saveflag_ori},'_gray']);
                imwrite((90+data4)/180,...
                    output_filepath,...
                    'compression','none');
            end
            
            map3D=ones([size(data4) 3]);
            map3D(:,:,1)=(90+data4)/180;
            maprgb=hsv2rgb(map3D);
        %     figure,imshow(data4);
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_ori});
                imwrite(maprgb,...
                    output_filepath,...
                    'compression','none');
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName4);
            end
        end
    end

    % Birefringence tile
    if any(saveflag_biref)
        VolName5 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_biref});
        if isfile(VolName5) 
            data5 = load_input(VolName5,transpose_flag);
            data5 = (data5(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data5,birefGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_biref});
                imwrite(mat2gray(data5,birefGrayRange),...
                    output_filepath,...
                    'compression','none');
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName5);
            end
        end  
    end

    % Mus tile
    if any(saveflag_mus)
        VolName6 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_mus});
        if isfile(VolName6) 
            data6 = load_input(VolName6,transpose_flag);
            data6 = (data6(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data6,musGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_mus});
                imwrite(mat2gray(data6,musGrayRange),...
                    output_filepath,...
                    'compression','none');
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName6);
            end
        end  
    end

    %Surf tile
    if any(saveflag_surf)
        VolName7 = get_volname(input_file_format,mosaicid(i),tileid(i),modality{saveflag_surf});
        if isfile(VolName7) 
            data7 = load_input(VolName7,transpose_flag);
            data7 = (data7(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data7,surfGrayRange);colormap gray;
            if doSaveTif
                output_filepath = replace(sprintf(output_file_format,mosaicid(i),tileid(i)),'[modality]',modality{saveflag_surf});
                imwrite(mat2gray(data7,surfGrayRange)*255,jet,...
                    output_filepath,...
                    'compression','none');  
            end
        else
            if doSaveTif
                fprintf('Missing %s\n',VolName7);
            end
        end
    end
end
delete(gcp('nocreate'));
end


function VolName = get_volname(BaseFileName,mosaic_num,tile_num,modality)
% this function is for the flexibility
% during atypical processing, input file naming was different from the
% standard. e.g. location of num and modality swapped; different modality
% naming; 
% this function provide flexibility when these situation happens
% might still need improvement
    VolName = sprintf(BaseFileName,mosaic_num,tile_num);
    VolName = replace(VolName,'[modality]',modality);
end

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





