function Mosaic2D_Telesto(ParameterFile, modality, method)

% Mosaic2D_Telesto(ParameterFile, modality, method)
%   stitch .nii or .mat tiles into slices for specified 2D/en-face modalities
% 
% USAGE:
% 
%  INPUTS:          
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters, Scan, and Mosaic2D: 
%       modality          =  OCT modality to process
%       method            =  Select coordinate that is calculated from the method.
%                            Output will be written to subfolder of the directory
%                           - Method_1 : median coordinates in z
%                           - Method_2 : median stepsize and offset
% 
%   ~2024-03-06~  
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIAL SETTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(ParameterFile); %Parameters Mosaic2D Scan


%%% if using resize_img 
% resize_img  = Mosaic2D.imresize;

%%% if method specified 
method_select = 0;
if exist('method','var'); method_select = 1; end 
%%% display numbers of workers 
% fprintf('--nWorker is %i--\n', nWorker);

%%% display modality 
fprintf('--Modality is %s--\n', modality);

if isdeployed
    disp('--App is deployed--');
else
    disp('--App is NOT deployed--');
    % addpath(genpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils/OCTBasic'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING SLICE INPUT 
sliceidx    = Mosaic2D.sliceidx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING MOSAIC PARAMETERS 

fprintf(' - Loading Experiment file...\n %s\n', Mosaic2D.Exp);
S = whos('-file',Mosaic2D.Exp);
if any(contains({S(:).name},'Experiment_Fiji'))
    idx = find(contains({S(:).name},'Experiment_Fiji'));    isFiji = true;
elseif any(contains({S(:).name},'Experiment'))
    idx = find(contains({S(:).name},'Experiment'));         isFiji = false;
end
load(Mosaic2D.Exp,S(idx).name);
Experiment  = eval(S(idx).name);
fprintf(' - %s is loaded ...\n %s\n', S(idx).name);

if method_select==1
    Experiment.X_Mean = Experiment.(method).X_Mean;
    Experiment.Y_Mean = Experiment.(method).Y_Mean;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir       = Mosaic2D.indir;
outdir      = Mosaic2D.outdir; 

if method_select==1; outdir = [outdir '/' method]; end
if ~exist(outdir,'dir'); mkdir(outdir); end

filetype    = Mosaic2D.InFileType; % 'nifti';
% switch filetype
%     case 'nifti'; % 2d imag file format -> Imag=[indir_curr filesep 'test_processed_' sprintf('%03i',currtile) '_' lower(modality) '.nii'];
%     case 'mat';  % load 2d file like during the Enface step 
% end

transpose_flag = Parameters.transpose;

%fprintf(' - Input directory = \n%s\n',indir{:});
fprintf(' - Output directory = \n%s\n',outdir);
if method_select==1; warning(' - > %s < Method specified \n',method);end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING IMAGE GRAYSCALE RANGE
GrayRange = [];
savetiff  = 0;
savenii   = 1;
switch lower(modality)
    case 'aip';         if isfield(Parameters,'AipGrayRange'); GrayRange = Parameters.AipGrayRange; else; noscale = 0;end
    case 'mip';         if isfield(Parameters,'MipGrayRange'); GrayRange = Parameters.MipGrayRange; else; noscale = 0;end
    case 'retardance';  if isfield(Parameters,'RetGrayRange'); GrayRange = Parameters.RetGrayRange; else; noscale = 0;end
    case 'birefringence';if isfield(Parameters,'BirefGrayRange'); GrayRange = Parameters.BirefGrayRange; else; noscale = 0;end
    case 'mus';         if isfield(Parameters,'musGrayRange'); GrayRange = Parameters.musGrayRange; else; noscale = 0;end
    case 'surf';        if isfield(Parameters,'surfGrayRange'); GrayRange = Parameters.Surfacerange; else; noscale = 0;end    
    otherwise;          if isfield(Parameters,[modality 'GrayRange']); GrayRange = Parameters.([modality 'GrayRange']); else; noscale = 0;end %disp(' - unrecognized modality!');
end
if exist('noscale', 'var')
    savetiff = 0; 
    warning(' %s Grayscale is not found. Will only output MAT file. \n', modality);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MODALITY string when inconsistency happens at file name
modality_base        = Enface.save;
str_idx = contains(modality_base, modality(1:3),'IgnoreCase',true);
if any(str_idx)
    modality_str = modality_base{str_idx};
else
    modality_str = modality(1:3);
    warning(' %s (current) modality is not included in Enface struct. Mosaic2D might fail. \n', modality);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIMENSION OF STITCHED IMAGES
nslices = size(sliceidx,2);
if(strcmpi(Scan.TiltedIllumination,'Yes') == 1)
    if method_select==1
        Experiment.X_Mean_tilt = Experiment.(method).X_Mean_tilt;
        Experiment.Y_Mean_tilt = Experiment.(method).Y_Mean_tilt;
    end
    mosaics_per_slice = 2;
    mosaic_nums{1,1} = sliceidx(1,1)*2 - 1:2:sliceidx(1,end)*2 - 1;
    mosaic_nums{2,1} = sliceidx(1,1)*2:2:sliceidx(1,end)*2;
else
    mosaics_per_slice = 1;
    mosaic_nums{1,1} = sliceidx(1,:);
end

if strcmpi(modality,'orientation');MZL=4;else;MZL=1;end

for n = 1:mosaics_per_slice
    switch n
        case 1
            mod_str{n,1} = '';
            NbPix{n,1} = [Experiment.NbPix,Experiment.NbPix]; % [y,x]

            X{n,1}       = Experiment.X_Mean;           Y{n,1}       = Experiment.Y_Mean;         % Y = fliplr(Y);
            X{n,1}       = X{n,1}-min(X{n,1}(:))+1;               Y{n,1}      = Y{n,1}-min(Y{n,1}(:))+1;
            sizerow{n,1} = NbPix{n,1}(1,2)-XPixClip;              sizecol{n,1} = NbPix{n,1}(1,1)-YPixClip;
            MXL{n,1}     = max(X{n,1}(:))+sizerow{n,1}-1;         MYL{n,1}     = max(Y{n,1}(:))+sizecol{n,1}-1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% GENERATE RAMP

            % ramp_x{n,1} = sizerow{n,1}-round(median(diff(Experiment.X_Mean,[],1),'all','omitnan'));
            ramp_x{n,1} = sizerow{n,1};
            % ramp_y{n,1} = sizecol{n,1}-round(median(diff(Experiment.Y_Mean,[],2),'all','omitnan'));
            ramp_y{n,1} = sizecol{n,1};

            ramp_xv{n,1}       = 0:ramp_x{n,1}-1;            ramp_yv{n,1}       = 0:ramp_y{n,1}-1;
            x{n,1}             = ones(1,sizerow{n,1});       y{n,1}             = ones(1,sizecol{n,1});
            x{n,1}(1+ramp_xv{n,1})  = mat2gray(ramp_xv{n,1});     y{n,1}(1+ramp_yv{n,1})  = mat2gray(ramp_yv{n,1});
            x{n,1}(end-ramp_xv{n,1})= mat2gray(ramp_xv{n,1});     y{n,1}(end-ramp_yv{n,1})= mat2gray(ramp_yv{n,1});

            RampOrig{n,1}=x{n,1}.'*y{n,1};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set mosaic params
            %MapIndex = Experiment.MapIndex_Tot_offset+ Experiment.First_Tile - 1;


            tabrow{n,1} = size(Experiment.MapIndex_Tot,1);
            tabcol{n,1} = size(Experiment.MapIndex_Tot,2);

            MapIndex{n,1} = Experiment.MapIndex_Tot_offset+Experiment.First_Tile -1;
        case 2
            mod_str{n,1} = '_tilt';
            NbPix{n,1} = [Experiment.NbPix,Experiment.NbPix_tilt]; % [y,x]

            X{n,1}       = Experiment.X_Mean_tilt;           Y{n,1}       = Experiment.Y_Mean_tilt;         % Y = fliplr(Y);
            X{n,1}       = X{n,1}-min(X{n,1}(:))+1;               Y{n,1}      = Y{n,1}-min(Y{n,1}(:))+1;
            sizerow{n,1} = NbPix{n,1}(1,2)-XPixClip;              sizecol{n,1} = NbPix{n,1}(1,1)-YPixClip;
            MXL{n,1}     = max(X{n,1}(:))+sizerow{n,1}-1;         MYL{n,1}     = max(Y{n,1}(:))+sizecol{n,1}-1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% GENERATE RAMP

            % ramp_x{n,1} = sizerow{n,1}-round(median(diff(Experiment.X_Mean_tilt,[],1),'all','omitnan'));
            ramp_x{n,1} = sizerow{n,1};
            % ramp_y{n,1} = sizecol{n,1}-round(median(diff(Experiment.Y_Mean_tilt,[],2),'all','omitnan'));
            ramp_y{n,1} = sizecol{n,1};

            ramp_xv{n,1}       = 0:ramp_x{n,1}-1;            ramp_yv{n,1}       = 0:ramp_y{n,1}-1;
            x{n,1}             = ones(1,sizerow{n,1});       y{n,1}             = ones(1,sizecol{n,1});
            x{n,1}(1+ramp_xv{n,1})  = mat2gray(ramp_xv{n,1});     y{n,1}(1+ramp_yv{n,1})  = mat2gray(ramp_yv{n,1});
            x{n,1}(end-ramp_xv{n,1})= mat2gray(ramp_xv{n,1});     y{n,1}(end-ramp_yv{n,1})= mat2gray(ramp_yv{n,1});

            RampOrig{n,1}=x{n,1}.'*y{n,1};
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set mosaic params
            %MapIndex = Experiment.MapIndex_Tot_offset+ Experiment.First_Tile - 1;


            tabrow{n,1} = size(Experiment.MapIndex_Tot_tilt,1);
            tabcol{n,1} = size(Experiment.MapIndex_Tot_tilt,2);

            MapIndex{n,1} = Experiment.MapIndex_Tot_offset_tilt+Experiment.First_Tile_tilt -1;
    end

end

sample = single(1); % double(1); % single(1);
sample_whos = whos('sample');
fprintf(' - Processing in %s\n',sample_whos.class);
fprintf(' - Saving in %s\n','single');

%poolobj = parpool(nWorker);
%parfor s = 1:nslices
for s = 1:nslices
    sliceid_in  = sliceidx(1,s);
    sliceid_out = sliceidx(2,s);
    sliceid_run = sliceidx(3,s);
    indir_curr = indir{sliceid_run}; % this code can accommdate datas from multiple runs
    BaseFileName = [indir_curr '/' Mosaic2D.file_format];

    % Start mosaic
    tic
    fprintf('\nStarting %s mosaic slice %d from run %d slice %d ...\n',modality,sliceid_out, sliceid_run, sliceid_in);

    for n = 1:mosaics_per_slice
        M   = zeros(MXL{n,1},MYL{n,1},MZL,'like',sample);
        Ma  = zeros(MXL{n,1},MYL{n,1},'like',sample);

        for ii=1:tabcol{n,1}
            %     fprintf('Column %d\n',ii);
            for jj=1:tabrow{n,1}
                %             Mosaic  =zeros(MXL,MYL,MZL,'like',sample);
                %             Masque  =zeros(MXL,MYL,'like',sample);

                %         fprintf('row %d\n',jj);
                if MapIndex{n,1}(jj,ii)>0 && ~isnan(X{n,1}(jj,ii))
                    columns  = Y{n,1}(jj,ii):Y{n,1}(jj,ii)+sizecol{n,1}-1;
                    rows      = X{n,1}(jj,ii):X{n,1}(jj,ii)+sizerow{n,1}-1;

                    % currtile = (sliceid_in-1)*Experiment.TilesPerSlice+MapIndex(jj,ii);
                    currtile = MapIndex{n,1}(jj,ii); % New naming format re-cycles through tile numbers

                    % load/generate I (2D image)
                    switch filetype
                        case 'nifti'
                            if strcmpi(modality,'mus')
                                fprintf('2D mus code not updated');
                                return;
                                % Imag = get_volname(BaseFileName,mosaic_nums{n,1}(s),currtile,'cropped');
                                % I3 = niftiread(Imag);
                                % if ~isa(I3, 'single'); I3 = single(I3);end
                                % if canUseGPU(); I3 =  gpuArray(I3);end
                                % if strcmpi(Scan.System, 'Octopus'); I3 = permute(I3,[2,1,3]);end
                                % if transpose_flag;                  I3 = permute(I3,[2,1,3]);end
                                % 
                                % I3 = I3(XPixClip+1:end,YPixClip+1:end,:);
                                % sz = size(I3);
                                % sz(3) = 100; % average depth
                                % I3 = I3(:,:,end:-1:1);
                                % I3 = 10.^(I3/10);
                                % 
                                % data =  zeros(sz,'like',I3);
                                % for z=1:sz(3)-1
                                %     data(:,:,z) = I3(:,:,z)./(sum(I3(:,:,z+1:end),3))/2/0.0025;
                                % end
                                % if canUseGPU(); data =  gather(data);end
                                % I = squeeze(mean(data,3)); %figure;plot(squeeze(mean(data,[1 2])))
                            else % orientation, aip, mip, retardance, biref
                                Imag = get_volname(BaseFileName,mosaic_nums{n,1}(s),currtile,modality_str);
                                I = niftiread(Imag);
                                %if ~isa(I, 'single'); I = single(I);end
                                if strcmpi(Scan.System, 'Octopus'); I = I.';end
                                % if transpose_flag;  I = I.';end % Commented out because I transposed the RawData separately
                                I = I(XPixClip+1:end,YPixClip+1:end);
                            end
                        case 'mat'
                            Imag = get_volname(BaseFileName,mosaic_nums{n,1}(s),currtile,modality_str);
                            % Imag=[indir_curr filesep modality '_' sprintf('%03i',currtile) '.mat'];
                            S = whos('-file',Imag);
                            load(Imag);
                            I = eval(S.name);
                            if strcmpi(Scan.System, 'Octopus'); I = I.';end
                            if transpose_flag;                  I = I.';end
                            I = I(XPixClip+1:end,YPixClip+1:end);
                        otherwise
                            warning('filetype is not recognized: [ %s ]\n', filetype)
                            return
                    end

                    % blend I into Mosaic
                    if strcmpi(modality,'orientation')
                        %I = orientationSign * I + orientationOffset;
                        SLO_OC=zeros([size(I) 4],'double');
                        SLO_OC(:,:,1)=cos(I/180*pi).*cos(I/180*pi);
                        SLO_OC(:,:,2)=cos(I/180*pi).*sin(I/180*pi);
                        SLO_OC(:,:,3)=SLO_OC(:,:,2);
                        SLO_OC(:,:,4)=sin(I/180*pi).*sin(I/180*pi);

                        M(rows,columns,:)=M(rows,columns,:)+SLO_OC.*repmat(RampOrig{n,1},1,1,4);
                        Ma(rows,columns) = Ma(rows,columns)+RampOrig{n,1};
                    else % mus, aip, mip, retardance, birefringence
                        M(rows,columns) = M(rows,columns)+I.*RampOrig{n,1};
                        Ma(rows,columns) = Ma(rows,columns)+RampOrig{n,1};
                    end
                end
                %             M=M+Mosaic;
                %             Ma=Ma+Masque;
            end
        end

        switch lower(modality)
            case 'mip';         modalstr = 'MIP';
            case 'aip';         modalstr = 'AIP';
            case 'retardance';  modalstr = 'Retardance';
            case 'biref';       modalstr = 'Birefringence';
            case 'orientation'; modalstr = 'Orientation';
            case 'mus';         modalstr = 'mus';
            otherwise;          modalstr = modality;%error('foo:bar','Unknown image modality!');
        end


        if strcmpi(modality, 'orientation')
            fprintf('Starting orientation angles eigen decomp...\n');

            %%% orientation angles eigen decomposition
            M=M./(Ma);
            M = double(M);
            a_x = reshape(M,[size(M,1)*size(M,2),2,2]);
            [V,~]= eig2(a_x);
            Data_ = atan(V(:,2,2)./V(:,1,2));
            Data_ = reshape(Data_, size(M,[1 2]));
            Data  = Data_/pi*180;
            O     = Data;

            index1=O<-90;
            index2=O>90;
            O(index1)=O(index1)+180;
            O(index2)=O(index2)-180;
            O=rot90(O,-1);
            MosaicFinal=O;



            %%% save MosaicFinal
            if ~isa(MosaicFinal, 'single'); MosaicFinal = single(MosaicFinal);end
            fprintf('Saving Orientation mosaic .mat...\n');
            parsave([outdir,'/',modalstr,mod_str{n,1},sprintf('_slice%03i.mat',sliceid_out)],MosaicFinal); % ,'-v7.3');

            %%% masking orientation
            if(savetiff == 1)
                RetSliceTiff = [outdir,'/','Retardance',mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)];
                AipSliceTiff = [outdir,'/','AIP',mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)];
                data1 = imread(RetSliceTiff);
                data2 = imread(AipSliceTiff);
            elseif(savenii == 1)
                RetSliceNii = [outdir,'/','Retardance',mod_str{n,1},sprintf('_slice%03i.nii',sliceid_out)];
                AipSliceNii = [outdir,'/','AIP',mod_str{n,1},sprintf('_slice%03i.nii',sliceid_out)];
                data1 = mat2gray(niftiread(RetSliceNii),Parameters.AipGrayRange);
                data2 = mat2gray(niftiread(AipSliceNii),Parameters.AipGrayRange);
            end
            data4 = wiener2(data2,[5 5]);

            map3D=ones([size(O,[1,2]) 3]);
            map3D2=ones([size(O,[1,2]) 3]);
            map3D3=ones([size(O,[1,2]) 3]);
            O_normalized = (O+90)/180;

            %%% Orientation1
            I=mat2gray(double(data1))/(1-0.4);
            %    I=(double(data1)-80)/60; % -lowerthreshold) / DynamicRange
            I(I<0)=0;
            I(I>1)=1;
            I(data4<=20)=0;
            map3D(:,:,1)=O_normalized;
            map3D(:,:,3)=I;
            maprgb=hsv2rgb(map3D);
            %         figure; imshow(maprgb); title('Orientation 1');
            fprintf('Saving Orientation1 mosaic .tiff...\n');
            imwrite(maprgb,[outdir,'/',modalstr,'1',mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)],'compression','none');

            %%% Orientation2
            I=(-mat2gray(double(data2))+1-0.2)/((1-0.2)*0.5); % -() +1 invert; -0.2 clip graymatter; (1-0.2) then dynamic range; 0.5 bring shadow to midtone
            %    I=(double(data2)+20)/220;
            I(I<0)=0;
            I(I>1)=1;
            I(data4<=20)=0;
            map3D2(:,:,1)=O_normalized;
            map3D2(:,:,3)=I;
            maprgb2=hsv2rgb(map3D2);
            %     figure; imshow(maprgb2); title('Orientation 2');
            fprintf('Saving Orientation2 mosaic .tiff...\n');
            imwrite(maprgb2,[outdir,'/',modalstr,'2',mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)],'compression','none');
            
            map3D3(:,:,1)=O_normalized;
            maprgb=hsv2rgb(map3D3);
            imwrite(maprgb,[outdir,'/',modalstr,'3',mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)],'compression','none');
        else
            % If not Orientation, save tiff and mat for slice
            MosaicFinal=rot90(M./(Ma),-1);
            MosaicFinal(isnan(MosaicFinal))=0;

            if(savetiff == 1)
                fprintf('Saving .tiff mosaic...\n');
                foutimg = [outdir,'/',modalstr,mod_str{n,1},sprintf('_slice%03i.tiff',sliceid_out)];
                imwrite(mat2gray(MosaicFinal,GrayRange),foutimg,'compression','none');
            elseif(savenii == 1)
                fprintf('Saving .nii mosaic...\n');
                foutimg = [outdir,'/',modalstr,mod_str{n,1},sprintf('_slice%03i.nii',sliceid_out)];
                niftiwrite(MosaicFinal,foutimg);
            end

            if ~isa(MosaicFinal, 'single'); MosaicFinal = single(MosaicFinal);end
            fprintf('Saving .mat mosaic...\n');
            foutmat = [outdir,'/',modalstr,mod_str{n,1},sprintf('_slice%03i.mat',sliceid_out)];
            parsave(foutmat, MosaicFinal);

        end
    end

    toc
    fprintf('\n   ---   Finished with %s slice %d   ---   \n',modalstr,sliceid_out);

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
