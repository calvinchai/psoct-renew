function [dBI3D, R3D, O3D] = Save3D_tile_20250429(filename, FileNum, dispCompFile, Aline_length, Bline_length)

% Save 3D tile from spectral data
% 
% INPUT ARGUMENTS:
%  Input_Path       = dir with spectral tiles
%  FileNum          = tile number
%  opts             = optional; fields:
%      Output_Path  = path to save (def. Input_Path)
%      modalities   = modalities to save (def. dBI, ret, orien)
% 

% if nargin<3
%     disp('USAGE: Save3D_tile( Input_Path, FileNum, Output_Path );');
%     return
% end

if isdeployed
else
    addpath('/autofs/cluster/octdata2/users/Hui/PSCalibration/code');
    addpath('/autofs/cluster/octdata2/users/Hui/tools/rob_utils');
    addpath('/autofs/cluster/octdata2/users/Chao/code/telesto');
    addpath('/autofs/cluster/octdata2/users/Chao/code/tools/freesurfer')
end

tic
warning off MATLAB:polyfit:RepeatedPointsOrRescale

% clear all;
% % close all;
% clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%% Input_path is where tiles are located 
%%%%%% Output_path is where to write to
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other params to set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%% px = number of pixels (assumes same number of pixels in X and Y
%%%%%%   (if different px in X and Y, then change to appropriate numbers in
%%%%%%   the call to Aline=px; Bline=px; )
%
%%%%%% FileNum = the number of the tile to process (1 = 001, etc.)
% 
%%%%%% phaseDelays = the phase delays to process data with. Can just set to
%%%%%%                  one value to only use one phase delay
%%%%%%       (You can set this to the value you used for calibration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Aline = Aline_length;
Bline = Bline_length;
fprintf('aline = %i; bline = %i \n',Aline,Bline)

MakePlots = 0;
saveenface = 0 ;
saveflag = 0;

SurfaceZOffset = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BaseFileName = [Input_Path filesep 'test_spectral_'];%300:350 % %loop blines
% FileNum= [1]; %[1];
% depthstag = [sprintf('%d',FileNum(1)) ',' sprintf('%d',FileNum(end))];
% NumberFiles=length(FileNum);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% HERE YOU NEED TO PUT THE CALIBRATION FILES TO USE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Spacing = 1e-3/600;
FocalLength = 100*1e-3;
FL=FocalLength;
FocalLengthInterval = 0.2*1e-3;

CenterWavelength = 1300*1e-9;
lambda0=CenterWavelength;
CenterWavelengthInterval = 10*1e-9;

xo = 0;
xoInterval = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get interpolated wavelengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\\
OriginalLineLength1 = 2048;
OriginalLineLength2 = 2048;
Start1=1;
Start2=1;
PaddingFactor=1;
PaddingLength = 2048*PaddingFactor;
InterpolationParameters = [PaddingFactor,PaddingLength,OriginalLineLength1,Start1,OriginalLineLength2,Start2];
[Wavelengths_l, Wavelengths_r,InterpolatedWavelengths2, Ks] = interpolationwave_101620 (InterpolationParameters);
% [Wavelengths_l, Wavelengths_r,InterpolatedWavelengths2, Ks] = interpolationwave_111621 (InterpolationParameters);

InterpolatedWavelengths = (InterpolatedWavelengths2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of values for FocalLength, CenterWavelength, IncAngle and xo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numF = 1;
numWave = 1;
numInc = 1;
numxo = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AutoCorrPeakCut = 24; % multiple of 8
AlineLength = 2048;
DepthL=  1024- AutoCorrPeakCut; %%%
offset = 100/180*pi;

dBOffset = 0 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dispCompFile1 = dispCompFile;
dispCompFile2 = dispCompFile;

fiddc1 = fopen(dispCompFile1, 'rb');
phaseDispersion1 = fread(fiddc1, inf, 'real*8');
fclose(fiddc1);
phaseCorrection1 = exp(-1i .* reshape(phaseDispersion1,AlineLength*PaddingFactor,[]));
phaseCorrection1 = repmat(phaseCorrection1, 1, Aline);
% size(phaseCorrection1)
fiddc2 = fopen(dispCompFile2, 'rb');
phaseDispersion2 = fread(fiddc2, inf, 'real*8');
fclose(fiddc2);
phaseCorrection2 = exp(-1i .* reshape(phaseDispersion2,AlineLength*PaddingFactor,[]));
phaseCorrection2 = repmat(phaseCorrection2, 1, Aline);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read required number of A-lines from each file and average the same
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for groupIndex =    1 :   length(FileNum) %loop file numbers
    fprintf('\n -- -- TileNb / FileNum %g -- \n\n',FileNum(groupIndex));
    
    dBI3D=zeros(Bline,Aline,DepthL,'single');
%     dBIm3D=zeros(Bline,Aline,DepthL);
%     dBIx3D=zeros(Bline,Aline,DepthL);
    R3D=zeros(Bline,Aline,DepthL,'single');
    O3D=zeros(Bline,Aline,DepthL,'single');

%     surface = zeros (Bline,Aline) ;
    
    for FileInd =   1:  Bline 
        if mod(FileInd,50)==0, disp([' bline ' num2str(FileInd)]); end
        
        % fid=fopen([BaseFileName sprintf('%03i', FileNum(groupIndex)) '.nii'], 'rb');
        fid=fopen(filename, 'rb');
        fseek(fid,(FileInd-1)*AlineLength*Aline*2*2+352,'bof');
        data1=fread(fid, AlineLength*Aline, 'uint16');
        fclose(fid);
        WavelengthBuffer1=reshape(data1,2048,[]);
        
        % fid=fopen([BaseFileName sprintf('%03i', FileNum(groupIndex)) '.nii'], 'rb');
        fid=fopen(filename, 'rb');
        fseek(fid,(FileInd-1)*AlineLength*Aline*2*2+AlineLength*Aline*2+352,'bof');
        data2=fread(fid, AlineLength*Aline, 'uint16');
        fclose(fid);
        WavelengthBuffer2=flipud(reshape(data2,2048,[]));
        BufferNum = 1;    % Buffer to read from each file. Any num from 1-10

        refdata1=mean(WavelengthBuffer1,2);
        refdata2=mean(WavelengthBuffer2,2);
        MeanScan1=WavelengthBuffer1-repmat(refdata1,1,Aline);
        MeanScan2=WavelengthBuffer2-repmat(refdata2,1,Aline);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read the validated part of the two polorized spectrums
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        OriginalBuffer1 = MeanScan1(Start1:OriginalLineLength1-1+Start1,:);
        OriginalBuffer2 = MeanScan2(Start2:OriginalLineLength2-1+Start2,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Iterate through all the values of varying interpolation parameters and
        % fit the phase curves and find the errors for each case
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % zeropad of the two spectrums to the same size, and denser in points to
        % PaddingLength

        ZeroPaddedBuffer1 = ZeroPadBuffer(OriginalBuffer1, PaddingFactor);
        ZeroPaddedBuffer2 = ZeroPadBuffer(OriginalBuffer2, PaddingFactor);

        AlineLength1 = size(OriginalBuffer1, 1);
        AlineLength2 = size(OriginalBuffer2, 1);
        NumberAlines = size(OriginalBuffer1, 2);

        % figure; plot(OriginalBuffer);title('Original spectrums');


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        InterpolatedBuffer1 = interp1(Wavelengths_l, ZeroPaddedBuffer1, InterpolatedWavelengths2,'linear','extrap');
        InterpolatedBuffer2 = interp1(Wavelengths_r, ZeroPaddedBuffer2, InterpolatedWavelengths2,'linear','extrap');
        InterpolatedBuffer1 = InterpolatedBuffer1 - repmat(median(InterpolatedBuffer1,2),1,Aline);
        InterpolatedBuffer2 = InterpolatedBuffer2 - repmat(median(InterpolatedBuffer2,2),1,Aline);
        % size(InterpolatedBuffer1)
        % size(phaseCorrection1)
        InterpolatedBuffer1=InterpolatedBuffer1.*phaseCorrection1;
        InterpolatedBuffer2=InterpolatedBuffer2.*phaseCorrection2;
        
        Jones1 = Buffer2Jones(InterpolatedBuffer1, PaddingFactor, AutoCorrPeakCut);
        Jones2 = Buffer2Jones(InterpolatedBuffer2, PaddingFactor, AutoCorrPeakCut);
        IJones = abs(Jones1).^2+abs(Jones2).^2;
        IJones1 = abs(Jones1).^2;
        IJones2 = abs(Jones2).^2;
        dBI=10*log10(IJones);

        dBI1=10*log10(IJones1); % Blue channel
        dBI2=10*log10(IJones2); % Green channel

        retardance=atan(abs(Jones1)./abs(Jones2))/pi*180; % tissue top should be dark
%         retardance=atan(medfilt2(abs(Jones2),[3 3])./medfilt2(abs(Jones1),[3 3]))/pi*180;
        phase1 = (angle(Jones1));
        phase2 = (angle(Jones2));
        phi=(phase1-phase2)+offset*2;Jones1;
        index1=phi>pi;
        phi(index1)=phi(index1)-2*pi;
        index2=phi<-pi;
        phi(index2)=phi(index2)+2*pi;
        orientation = (phi)/2/pi*180;

        dBImage=10*log10(IJones)-dBOffset;     

        if saveenface ==1  
            
            dBImage1=dBImage.*mask2;
            %             figure,imagesc(dBImage1)
            Iedge = edge(medfilt2(dBImage1,[2 2]), 'canny',[0.10 0.75],10); %%% detection of top surface
            %             figure,imagesc(Iedge)
            [edgeMax,edgeIndex]=max(Iedge);
            mask3=zeros(size(dBImage,1),size(dBImage,2));

                for i=1:size(dBImage,2)

                        if (Istart(i)~=0 && Iend(i)-Istart(i)>15)
                            if edgeIndex(i)~=1
                                Istart(i)=edgeIndex(i);
                                %                     mask3(edgeIndex(i):Iend(i),i)=1;
                            else
                                %                     mask3(Istart(i):Iend(i),i)=1;
                            end
           Istart(Istart<10)=100; 
           Istart(Istart>700)=100; 
                                
                                if Istart(i)+DepthL-1<=size(mask2,1)
                                     dBI3D(FileInd,i,:) = dBImage(Istart(i):Istart(i)+DepthL-1,i);
%                                    dBI3D(FileInd,i,:) = dBImage(Istart(i):Istart(i)+DepthL-1,i);
            %                         dBIm3D(FileInd,i,:) = dBI1(Istart(i):Istart(i)+DepthL-1,i);
            %                         dBIx3D(FileInd,i,:) = dBI2(Istart(i):Istart(i)+DepthL-1,i);
                                    R3D(FileInd,i,:) = retardance(Istart(i):Istart(i)+DepthL-1,i);
                                    O3D(FileInd,i,:) = orientation(Istart(i):Istart(i)+DepthL-1,i);
                                else
                                    dBI3D(FileInd,i,:) = dBImage(Istart(i):size(dBImage,1),i);
            %                         dBIm3D(FileInd,i,:) = dBI1(Istart(i):size(dBImage,1),i);
            %                         dBIx3D(FileInd,i,:) = dBI2(Istart(i):size(dBImage,1),i);
                                    R3D(FileInd,i,:) = retardance(Istart(i):size(dBImage,1),i);
                                    O3D(FileInd,i,:) = orientation(Istart(i):size(dBImage,1),i);
                                end

                        end
                end

        end 
    

        dBI3D(FileInd,:,:) = dBImage(:,:)';
        R3D(FileInd,:,:) = retardance(:,:)';
        O3D(FileInd,:,:) = orientation(:,:)';
    end        % for FileInd
    

    if MakePlots
        
        ff=figure('position',[672 466 1164 486]);
        subplot(131);
        imagesc(dBI(1:500,:),[35 75]); set(gca,'colormap',gray); title(['dBI - bline ' num2str(FileInd)]);
         hold on; plot(Istart,'LineWidth', 2); 
         hold on; plot(Istart+DepthL,'LineWidth', 2)

        subplot(132);
        imagesc(retardance(1:500,:),[0 90]); set(gca,'colormap',jet); title('ret');
         hold on; plot(Istart,'LineWidth', 2); 
         hold on; plot(Istart+DepthL,'LineWidth', 2)
        subplot(133);
        imagesc(orientation(1:500,:),[-90 90]); set(gca,'colormap',hsv); title('orien');
         hold on; plot(Istart,'LineWidth', 2); 
         hold on; plot(Istart+DepthL,'LineWidth', 2)
        
    end
    
    if saveflag == 1
%             indd = FileNum(groupIndex);
%                 if indd >270 
%                     indd = indd-1;
%                 end
            
            
            save([Output_Path2 filesep 'AIP_'     sprintf('%03i', FileNum(groupIndex)) '.mat'],'SLO_dBI');
            save([Output_Path2 filesep 'MIP_'     sprintf('%03i', FileNum(groupIndex)) '.mat'],'SLO_dBI2');
            save([Output_Path2 filesep 'Ret_'     sprintf('%03i', FileNum(groupIndex)) '.mat'],'SLO_R');
            save([Output_Path2 filesep 'Orien_'     sprintf('%03i', FileNum(groupIndex)) '.mat'],'SLO_O3');
            save([Output_Path2 filesep 'Bi_'     sprintf('%03i', FileNum(groupIndex)) '.mat'],'SLO_bi');
          

            save([Output_Path filesep 'surface_'   sprintf('%03i', FileNum(groupIndex)) '.mat'],'surface');
             
            save([Output_Path3 filesep 'dBI3D_'   sprintf('%03i', FileNum(groupIndex)) '.mat'],'dBI3D');
            save([Output_Path3 filesep 'Ret3D_'   sprintf('%03i', FileNum(groupIndex)) '.mat'],'R3D');
            save([Output_Path3 filesep 'Orien3D_' sprintf('%03i', FileNum(groupIndex)) '.mat'],'O3D'); 

                imgR1=(SLO_R-15)/40;
                imgR1(imgR1<0)=0;imgR1(imgR1>1)=1;
            %     figure,imagesc(flipud(rot90(imgR1)))
%                     imgR1=flipud((imgR1));

                phi=SLO_O3;
                index1=phi>180;
                phi(index1)=phi(index1)-360;
                index2=phi<-180;
                phi(index2)=phi(index2)+360;
            %     figure,imagesc(flipud(rot90(phi)));colormap hsv

                minO=-90;maxO=90;
                imgOA1=(phi-minO)/(maxO-minO);
                imgOA1(imgOA1>1)=1;imgOA1(imgOA1<0)=0; 

                HMap=ones([size(imgOA1),3]);
                HMap(:,:,1)=imgOA1;
                HMap(:,:,3)=imgR1;
                RGBMap1=hsv2rgb(HMap);
                
            
                imgR1=(SLO_dBI-5)/10;
                imgR1(imgR1<0)=0;imgR1(imgR1>1)=1;
                HMap=ones([size(imgOA1),3]);
                
                HMap(:,:,1)=imgOA1;
                HMap(:,:,3)=imgR1;
                RGBMap2=hsv2rgb(HMap);
                                
               imwrite(mat2gray((SLO_dBI),[10 35]),[Output_Path 'AIP_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               imwrite(mat2gray(SLO_dBI2,[15 40]),[Output_Path 'MIP_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               imwrite(mat2gray((SLO_R),[20 50]),[Output_Path 'Retardance_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               imwrite(RGBMap1,[Output_Path 'Orientation1_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               imwrite(RGBMap2,[Output_Path 'Orientation2_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               
               imwrite(imgOA1,[Output_Path 'Orientation_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
               imwrite(mat2gray((SLO_bi),[0 1]),[Output_Path 'Bi_' sprintf('%03i', FileNum(groupIndex)) '.tiff'],'compression','none');
        

%                    imwrite(mat2gray((SLO_dBI),[10 30]),[Output_Path 'AIP_' sprintf('%03i', indd) '.tiff'],'compression','none');
%                    imwrite(mat2gray(SLO_dBI2,[20 50]),[Output_Path 'MIP_' sprintf('%03i', indd) '.tiff'],'compression','none');
%                    imwrite(mat2gray((SLO_R),[10 50]),[Output_Path 'Retardance_' sprintf('%03i', indd) '.tiff'],'compression','none');
%                    imwrite(RGBMap,[Output_Path 'Orientation1_' sprintf('%03i', indd) '.tiff'],'compression','none');
%                    imwrite(RGBMap,[Output_Path 'Orientation2_' sprintf('%03i', indd) '.tiff'],'compression','none');

    end

end
    
   

toc
disp(' ');

end

