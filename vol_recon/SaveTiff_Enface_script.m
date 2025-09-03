function SaveTiff_Enface_script( ParameterFile, nWorker )
% 
% SaveTiff_Enface_Basic(ProcDir, startTile, nbTiles)
%   Save .nii tiles as .tiff files for all 2D/en-face contrasts
% 
% USAGE:
% 
%  INPUTS:          
%       ProcDir     =   Path to processing directory
%                           - expects to find these files here: 
%                               OctScanInfoLog.txt 
%                               OctGrayscaleWindows.txt 
%                           - will create /Enface_Tiff + /Enface_MAT
%                             subdirectories in this directory
%       startTile   =   First tile number to process (1-based)
%       nbTiles     =   Total number of tiles to process
%                        (will process from startTile:startTile+nbTiles-1)
% 
% 
% Original script by HW, modified by DG
%   ~2021-12-22~  
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

% display numbers of workers 
fprintf('--nWorker is %i--\n', nWorker);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MODALITY TO SAVE 
modality        = Enface.save;

saveAIP = any(strcmpi(modality, 'aip'));
saveMIP = any(strcmpi(modality, 'mip'));
saveRet = any(strcmpi(modality, 'Retardance'));
saveOri = any(strcmpi(modality, 'Orientation'));
savemus = any(strcmpi(modality, 'mus'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING PARAMETERS 
XPixClip        = Parameters.XPixClip;
YPixClip        = Parameters.YPixClip;

% parameters for tiff (255) image intensity scale
AipGrayRange    = Parameters.AipGrayRange;
MipGrayRange    = Parameters.MipGrayRange;
RetGrayRange    = Parameters.RetGrayRange;
orientationSign = Parameters.OrientationSign;
orientationOffset = Parameters.OrienOffset;

% parameters to identify agar tile 
threshold_pct   = Parameters.Agar.threshold_pct;
threshold_std   = Parameters.Agar.threshold_std;
tissue_thresh   = Parameters.Agar.gm_aip;

% a vector of tile id to save
tileid = Parameters.TileID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
filetype        = Enface.filetype;
indir           = Enface.indir;
%outdir          =[Scan.ProcessDir filesep 'Enface_Tiff'];
outdir          = Enface.outdir;
agardir         = Parameters.Agar.dir;
if ~exist(outdir,'dir'),   mkdir(outdir);    end
if ~exist(agardir,'dir'),  mkdir(agardir);   end


fprintf(' - Input directory = \n%s\n',indir);
fprintf(' - Output directory = \n%s\n',outdir);
fprintf(' - Agar directory = \n%s\n',agardir);

BaseFileName = [indir filesep Scan.FilePrefix];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Hard-coded vars
doSaveTif = true;
rotflag = 1;
if rotflag == 1; fprintf('Tile will be rotated.\n');end
fprintf('Input format is %s.\n', filetype);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% poolobj = parpool(nWorker);
for i = 1:length(tileid)   
    fprintf('%d',tileid(i));
    
    %AIP tile
    if saveAIP
        if strcmp(filetype, 'nifti'); VolName1 = [BaseFileName sprintf('%03i', tileid(i)) '_aip.nii'];
        elseif strcmp(filetype, 'mat'); VolName1 = [indir filesep 'AIP_' sprintf('%03i', tileid(i)) '.mat'];
        end
        if isfile(VolName1) 
            if strcmp(filetype, 'nifti'); data1 = readnifti(VolName1);
            elseif strcmp(filetype, 'mat'); load(VolName1,'SLO_dBI');data1 = SLO_dBI;
            end
            if rotflag == 1;data1 = rot90(data1,-1);end
            data1 = (data1(1+XPixClip:end,1+YPixClip:end));
            SLO_dBI = (data1);


            tissue_pct = nnz(SLO_dBI>tissue_thresh)/numel(SLO_dBI);
            agar_status = tissue_pct<threshold_pct & std(SLO_dBI,[],'all','omitnan')< threshold_std;
            tileinfo = [tileid(i),~agar_status];
            if agar_status==1; fprintf('-agarose');end
            parsave([agardir filesep sprintf('%03i', tileid(i))], tileinfo);

        %     figure,imagesc(data1,AipGrayRange);colormap gray;
            if doSaveTif
                imwrite(mat2gray(data1,AipGrayRange),...
                    [outdir filesep 'AIP_' sprintf('%03i', tileid(i)) '.tiff'],...
                    'compression','none');
            end
        else
            fprintf('Missing AIP %d\n',tileid(i));
        end
    end
    
    %MIP tile
    if saveMIP
        if strcmp(filetype, 'nifti'); VolName2 = [BaseFileName sprintf('%03i', tileid(i)) '_mip.nii'];
        elseif strcmp(filetype, 'mat'); VolName2 = [indir filesep 'MIP_' sprintf('%03i', tileid(i)) '.mat'];
        end
        if isfile(VolName2) 
            if strcmp(filetype, 'nifti'); data2 = readnifti(VolName2);
            elseif strcmp(filetype, 'mat'); load(VolName2,'SLO_dBI2');data2 = SLO_dBI2;
            end
            if rotflag == 1;data1 = rot90(data1,-1);end
            data2 = (data2(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data2,MipGrayRange);colormap gray;
            if doSaveTif
                imwrite(mat2gray(data2,MipGrayRange),...
                    [outdir filesep 'MIP_' sprintf('%03i', tileid(i)) '.tiff'],...
                    'compression','none');  
            end
        else
            fprintf('Missing MIP %d\n',tileid(i));
        end
    end
    
    %Retardance tile
    if saveRet
        if strcmp(filetype, 'nifti'); VolName3 = [BaseFileName sprintf('%03i', tileid(i)) '_retardance.nii'];
        elseif strcmp(filetype, 'mat'); VolName3 = [indir filesep 'Ret_' sprintf('%03i', tileid(i)) '.mat'];
        end
        
        if isfile(VolName3) 
            if strcmp(filetype, 'nifti'); data3 = readnifti(VolName3);
            elseif strcmp(filetype, 'mat'); load(VolName3,'SLO_R');data3 = SLO_R;
            end
            if rotflag == 1;data1 = rot90(data1,-1);end
            data3 = (data3(1+XPixClip:end,1+YPixClip:end));
        %     figure,imagesc(data3,RetGrayRange);colormap gray;
            if doSaveTif
                imwrite(mat2gray(data3,RetGrayRange),...
                    [outdir filesep 'Retardance_' sprintf('%03i', tileid(i)) '.tiff'],...
                    'compression','none');
            end
        else
            fprintf('Missing Retardance %d\n',tileid(i));
        end  
    end
    fprintf(' \n');
end

delete(gcp('nocreate'));
end









