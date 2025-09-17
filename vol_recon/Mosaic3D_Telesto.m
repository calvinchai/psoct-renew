function Mosaic3D_Telesto(ParameterFile,modality)

% Mosaic3D_Telesto(ParameterFile, modality)
%   stitch .nii or .mat tiles into 3D slice volumes
%
% USAGE:
%
%  INPUTS:
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters, Scan, and Mosaic3D:
%       modality          =  OCT modality to process
%
%   ~2024-06-13~
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(ParameterFile); %Parameters Mosaic3D Scan
ds_flag = false; % Flag to determine if we downsample
% modality = Mosaic3D.modality;
%%% display modality
fprintf('--Modality is %s--\n', modality);

%%%
UseGPUflag = false;
if UseGPUflag; fprintf('--GPU option is considered--\n');
else; fprintf('--GPU option is not considered--\n');end

if isdeployed
    disp('--App is deployed--');
else
    disp('--App is NOT deployed--');
    addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils/OCTBasic'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING SLICE INPUT
sliceidx    = Mosaic3D.sliceidx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS

fprintf(' - Loading Experiment file...\n %s\n', Mosaic3D.Exp);
S = whos('-file',Mosaic3D.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));    isFiji = true;
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));         isFiji = false;
end
load(Mosaic3D.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES
indir       = Mosaic3D.indir;
outdir      = [Mosaic3D.outdir '/' modality]; if ~exist(outdir,'dir'); mkdir(outdir); end
filetype    = Mosaic3D.InFileType; % 'nifti';

fprintf(' - Input directory = \n%s\n',indir{:});
fprintf(' - Output directory = \n%s\n',outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ADDING CLIPPING
switch Scan.System
    case 'Octopus'
        XPixClip    = Parameters.YPixClip; %18;
        YPixClip    = Parameters.XPixClip; %0;
    case 'Telesto'
        XPixClip    = Parameters.XPixClip; %18;
        YPixClip    = Parameters.YPixClip; %0;
end
fprintf('XPixClip = %i\nYPixClip = %i\nSystem is %s\n\n', XPixClip, YPixClip,Scan.System);

if ~isfield(Mosaic3D,'invert_Z')
    flag_flip = true;
else
    flag_flip = Mosaic3D.invert_Z;

end
if flag_flip; str_flip = '       '; else; str_flip = '  NOT  '; end
fprintf('crop volume will %s be Z flipped\n\n', str_flip);

% Set mosaic params
% MapIndex = Experiment.MapIndex_Tot_offset+ Experiment.First_Tile - 1;

% Set unused tiles to -1 in MapIndex
% MapIndex(Experiment_Fiji.MapIndex_Fiji==-1)=-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MODALITY string when inconsistency happens at file name
modality_base        = Processed3D.save;
str_idx = contains(modality_base, modality(1:3),'IgnoreCase',true);
if any(str_idx)
    modality_str = modality_base{str_idx};
else
    modality_str = modality(1:3);
    warning(' %s (current) modality is not included in Enface struct. Mosaic2D might fail. \n', modality);
end

% Set other parameters
MZL     = Mosaic3D.MZL;
fprintf('save %ipx of orignal %ipx\n', MZL, Scan.CropDepth);

sxyz    = Mosaic3D.sxyz;

nslices = size(sliceidx,2);
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    mosaics_per_slice = 2;
    mosaic_nums{1,1} = sliceidx(1,1)*2 - 1:2:sliceidx(1,end)*2 - 1;
    mosaic_nums{2,1} = sliceidx(1,1)*2:2:sliceidx(1,end)*2;
else
    mosaics_per_slice = 1;
    mosaic_nums{1,1} = sliceidx(1,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
for n = 1:mosaics_per_slice
    switch n
        case 1
            mod_str{n,1} = '';
            NbPix{n,1} = [Experiment.NbPix,Experiment.NbPix]; % [y,x]

            X{n,1}       = Experiment.X_Mean;           Y{n,1}  = Experiment.Y_Mean;         % Y = fliplr(Y);
            X{n,1}       = X{n,1}-min(X{n,1}(:))+1;     Y{n,1}  = Y{n,1}-min(Y{n,1}(:))+1;
            sizerow{n,1} = NbPix{n,1}(1,2)-XPixClip;              sizecol{n,1} = NbPix{n,1}(1,1)-YPixClip;
            MXL{n,1}     = max(X{n,1}(:))+sizerow{n,1}-1;         MYL{n,1}     = max(Y{n,1}(:))+sizecol{n,1}-1;
            %%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% GENERATE RAMP

            ramp_x{n,1} = sizerow{n,1}-round(median(diff(Experiment.X_Mean,[],1),'all','omitnan'));
            ramp_y{n,1} = sizecol{n,1}-round(median(diff(Experiment.Y_Mean,[],2),'all','omitnan'));

            ramp_xv{n,1}       = 0:ramp_x{n,1}-1;            ramp_yv{n,1}       = 0:ramp_y{n,1}-1;
            x{n,1}             = ones(1,sizerow{n,1});       y{n,1}             = ones(1,sizecol{n,1});
            x{n,1}(1+ramp_xv{n,1})  = mat2gray(ramp_xv{n,1});     y{n,1}(1+ramp_yv{n,1})  = mat2gray(ramp_yv{n,1});
            x{n,1}(end-ramp_xv{n,1})= mat2gray(ramp_xv{n,1});     y{n,1}(end-ramp_yv{n,1})= mat2gray(ramp_yv{n,1});

            RampOrig{n,1}=x{n,1}.'*y{n,1};
            RampOrig3{n,1} = repmat(RampOrig{n,1}, [1,1,MZL]);
            % %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            tabrow{n,1}=size(Experiment.MapIndex_Tot,1);
            tabcol{n,1}=size(Experiment.MapIndex_Tot,2);

            MapIndex{n,1} = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
        case 2
            mod_str{n,1} = '_tilt';
            NbPix{n,1} = [Experiment.NbPix,Experiment.NbPix_tilt]; % [y,x]

            X{n,1}       = Experiment.X_Mean_tilt;           Y{n,1}  = Experiment.Y_Mean_tilt;         % Y = fliplr(Y);
            X{n,1}       = X{n,1}-min(X{n,1}(:))+1;     Y{n,1}  = Y{n,1}-min(Y{n,1}(:))+1;
            sizerow{n,1} = NbPix{n,1}(1,2)-XPixClip;              sizecol{n,1} = NbPix{n,1}(1,1)-YPixClip;
            MXL{n,1}     = max(X{n,1}(:))+sizerow{n,1}-1;         MYL{n,1}     = max(Y{n,1}(:))+sizecol{n,1}-1;
            %%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% GENERATE RAMP

            ramp_x{n,1} = sizerow{n,1}-round(median(diff(Experiment.X_Mean_tilt,[],1),'all','omitnan'));
            ramp_y{n,1} = sizecol{n,1}-round(median(diff(Experiment.Y_Mean_tilt,[],2),'all','omitnan'));

            ramp_xv{n,1}       = 0:ramp_x{n,1}-1;            ramp_yv{n,1}       = 0:ramp_y{n,1}-1;
            x{n,1}             = ones(1,sizerow{n,1});       y{n,1}             = ones(1,sizecol{n,1});
            x{n,1}(1+ramp_xv{n,1})  = mat2gray(ramp_xv{n,1});     y{n,1}(1+ramp_yv{n,1})  = mat2gray(ramp_yv{n,1});
            x{n,1}(end-ramp_xv{n,1})= mat2gray(ramp_xv{n,1});     y{n,1}(end-ramp_yv{n,1})= mat2gray(ramp_yv{n,1});

            RampOrig{n,1}=x{n,1}.'*y{n,1};
            RampOrig3{n,1} = repmat(RampOrig{n,1}, [1,1,MZL]);
            % %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            tabrow{n,1}=size(Experiment.MapIndex_Tot_tilt,1);
            tabcol{n,1}=size(Experiment.MapIndex_Tot_tilt,2);

            MapIndex{n,1} = Experiment.MapIndex_Tot_offset_tilt+Experiment.First_Tile_tilt -1;
    end
end

% Loop slices
%poolobj = parpool(nWorker);
%parfor s = 1:nslices
for s = 1:nslices
    sliceid_in  = sliceidx(1,s);
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s);
    indir_curr = indir{sliceid_run};
    BaseFileName = [indir_curr '/' Mosaic3D.file_format];

    % Start mosaic
    tic0=tic;
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',modality,sliceid_out, sliceid_run, sliceid_in);

    for n = 1:mosaics_per_slice
        M  = zeros(MXL{n,1},MYL{n,1},MZL,'single'); %flipped xy for telesto
        Ma = zeros(MXL{n,1},MYL{n,1},'single');
        % Mosaic  = zeros(MXL{n,1},MYL{n,1},MZL,'single');
        % Masque = zeros(MXL{n,1},MYL{n,1},'single');

        for ii=1:tabcol{n,1}
            fprintf('col %d / %d\n',ii,tabcol{n,1})
            for jj=1:tabrow{n,1}
                %             fprintf('row %d / %d\n',jj,tabrow)

                if MapIndex{n,1}(jj,ii)>0 && ~isnan(X{n,1}(jj,ii))
                    columns =Y{n,1}(jj,ii):Y{n,1}(jj,ii)+sizecol{n,1}-1;
                    rows     =X{n,1}(jj,ii):X{n,1}(jj,ii)+sizerow{n,1}-1;

                    % currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                    currtile = MapIndex{n,1}(jj,ii); % New naming format re-cycles through tile numbers

                    switch filetype
                        case 'mat'
                            fname = get_volname(BaseFileName,mosaic_nums{n,1}(s),currtile,modality_str);
                            S = whos('-file',fname);
                            load(fname);
                            Imag = eval(S.name);
                        case 'nifti'
                            fname = get_volname(BaseFileName,mosaic_nums{n,1}(s),currtile,modality_str);
                            Imag =  niftiread(fname);
                    end
                    Imag = single(Imag);

                    %if [canUseGPU() & UseGPUflag]; Imag =  gpuArray(Imag);end

                    if strcmpi(Scan.System, 'Octopus');Imag = permute(Imag,[2,1,3]);end

                    Imag =  Imag(XPixClip+1:end,YPixClip+1:end,:); % clip
                    if flag_flip; Imag =  Imag(:,:,end:-1:1);end
                    sz =    size(Imag);

                    switch modality
                        case 'mus'
                            fprintf('mus 3D stitching not updated.');
                            return;
                            % switch filetype
                            %     case 'nifti'
                            %         I =     10.^(Imag/10);
                            % 
                            %         data =  zeros(sz(1),sz(2),MZL,'like',Imag);
                            %         for z=1:MZL
                            %             data(:,:,z) = I(:,:,z)./(sum(I(:,:,z+1:end),3))/2/0.0025;
                            %         end
                            %     case 'mat'
                            %         data = Imag(:,:,1:MZL);
                            % end

                        case 'dBI'
                            data = Imag(:,:,1:MZL);
                        case 'R3D'
                            data = Imag(:,:,1:MZL);
                        case 'O3D'
                            fprintf('Ori 3D stitching not updated.');
                            return;
                    end

                    M(rows,columns,:) = M(rows,columns,:)+data.*RampOrig3{n,1};
                    Ma(rows,columns) = Ma(rows,columns)+RampOrig{n,1};
                end
                % M=M+Mosaic;
                % Ma=Ma+Masque;
            end
        end

        Ma = repmat(Ma, [1,1,MZL]);

        % MosaicFinal=M./(Ma);
        % MosaicFinal(isnan(MosaicFinal)) = 0;
        % MosaicFinal(isinf(MosaicFinal)) = 0;
        M=M./(Ma);
        M(isnan(M)) = 0;
        M(isinf(M)) = 0;
        clear Ma
        
        if 1==0 % approve this flag when 'MosaicFinal' is smaller than 10 G
            disp(' - Saving mosaic...');
            fout = [outdir,'/',modality,mod_str{n,1},sprintf('_slice%03i.mat',sliceid_out)];
            fprintf(' %s \n',fout);
            aatic=tic;
            % save(fout, 'MosaicFinal', 'modality', '-v7.3');
            save(fout, 'M', 'modality', '-v7.3');
            %parsave(fout, MosaicFinal);
            aatoc=toc(aatic);
            fprintf(' - %.1f s to save\n',aatoc);
        end

        %     if strcmpi(modality, 'mus')
        %         fout = [outdir filesep '../StitchingFiji' filesep modality '_slice' sprintf('%03i',sliceid_out(ss)) '.mat'];
        %         mus_2d = mean(MosaicFinal(:,:,5:end),3);
        %         mus_2d = rot90(mus_2d,-1);
        %         save(fout, 'mus_2d');
        %     end


        %reset(gpuDevice())
        %%
        %%%% FOR DOWNSAMPLING

        sx = sxyz(1); % If 2: 10um -> 20um
        sy = sxyz(2); % If 2: 10um -> 20um
        sz = sxyz(3); % If 8: 2.5um -> 20um
        ds_xyz = [sx, sy, sz];

        output_size_xy_px = Scan.Resolution(1,1)*sx*1000; % micron
        output_size_z_px = Scan.Resolution(1,3)*sz*1000; % micron
        if(output_size_xy_px == output_size_z_px)
            ds_string = 'xyz';
        else
            ds_string = 'xy';
        end
        fprintf(' - Downsampling to %d micron iso (%s)...',output_size_xy_px,ds_string);

        % [Id] = do_downsample(MosaicFinal,ds_xyz);
        [Id] = do_downsample(M,ds_xyz);

        %         view3D(Id);
        fout = [outdir,'/',modality,mod_str{n,1},sprintf('_slice%03i_mean_%s_%ium_iso.mat',sliceid_out,ds_string,round(output_size_xy_px))];
        %save(fout, 'Id', 'modality', '-v7.3');
        %parsave(fout, Id);
        save(fout,'modality','-v7.3');
        save(fout,'Id','-append');



    end
    toc0=toc(tic0);
    disp(['Elapsed time to stitch ' num2str(length(sliceid_out)) ' slices: ']);
    disp(['      ' num2str(toc0/60) ' minutes']);
end
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