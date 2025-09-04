clear all
addpath ('/autofs/cluster/octdata2/users/Chao/code/demon_registration_version_8f')
addpath('/autofs/cluster/octdata2/users/Chao/code/telesto')
addpath ('/autofs/space/omega_001/users/3d_axis/PAPER/scripts')
addpath ('/autofs/cluster/octdata2/users/Chao/code/tools/cl_utils')

basename1 = '/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/run2_normal/StitchingFiji/';
basename2 = '/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/run2_tilt/StitchingFiji/';

for slice = 1
    slice
    % % read first set
    fixed1 = load([basename1 'Orientation_slice' sprintf('%03i', slice) '.mat']); % orientation fixed
    moving1 = load([basename2 'Orientation_slice' sprintf('%03i', slice) '.mat']);  % orientation moving

    fixed2 = load([basename1 'Birefringence_slice' sprintf('%03i', slice) '.mat']);
    moving2 = load([basename2 'Birefringence_slice' sprintf('%03i', slice) '.mat']);

    %   figure, subplot(1,2,1);imagesc(fixed_bi1);subplot(1,2,2);imagesc(moving_bi1)

%     f1_xrange1 = 11;             f1_xrange2 = 1600;%;1250; 
%     f1_yrange1 = 11;                  f1_yrange2 = 1400;%1200/2;
% 
%     m1_xrange1 = 421;          m1_xrange2 = 1900;%; 1330;
%     m1_yrange1 = 11;        m1_yrange2 = 1350;%1100/2;

    fixed_o1 =  rot90(rot90(fixed1.MosaicFinal));%(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
    fixed_bi1 =  rot90(rot90(fixed2.MosaicFinal));%(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
    moving_o1 = rot90(rot90(moving1.MosaicFinal));%(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);
    moving_bi1 = rot90(rot90(moving2.MosaicFinal));%(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);
    %   figure, subplot(1,2,1);imagesc(fixed_bi1);subplot(1,2,2);imagesc(moving_bi1)

    % fixed_o1 = imresize (fixed_o1, [size(fixed_o1,1) round(size(fixed_o1,2)/2)],'nearest');
    % fixed_bi1 = imresize (fixed_bi1, [size(fixed_bi1,1), round(size(fixed_bi1,2)/2)],'nearest');

    thruplane_reg_optiz_tensor_XY_final (fixed_bi1,moving_bi1, fixed_o1, moving_o1, 15, ['slice_test_' num2str(slice)])

clear fixed_bi1 moving_bi1 fixed_o1 moving_o1

end


%%%%%%%
   
num_iter = 20;option = 2;
delta_t = 1/100;kappa = 1;

for slice = 1:19
    
    load (['slice' num2str(slice) '_data.mat']);
   
    mask = medfilt2(biref_ObsLSQ, [5 5]);
    mask (mask>0.0015)=nan; mask(mask<0.00005)=nan; mask(~isnan(mask))=1;
    alpha = pi/2-Psi_ObsLSQ;
    alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
    alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;
    
    OC = angle2tensor (alpha/pi*180); % input in degrees
    OC1 = anisodiff2D(OC(:,:,1),num_iter,delta_t,kappa,option);
    OC2 = anisodiff2D(OC(:,:,2),num_iter,delta_t,kappa,option);
    OC3 = anisodiff2D(OC(:,:,3),num_iter,delta_t,kappa,option);
    OC4 = anisodiff2D(OC(:,:,4),num_iter,delta_t,kappa,option);
    alpha1 = tensor2angle (cat (3,OC1,OC2,OC3,OC4)); 
%     x = alpha1 (1400:1500,1200:1350); figure,hist(x(:)/pi*180,[-90:3:90])
    
    HMap=ones([size(alpha),3]);
    imgR1=(medfilt2( biref_ObsLSQ, [2 2])-0.0001)/0.0002;
    imgR1(imgR1<0)=0;imgR1(imgR1>1)=1;
    imgOA1=(alpha1+pi/2)/pi;
    imgOA1(imgOA1>1)=1;imgOA1(imgOA1<0)=0;
    HMap(:,:,1)=imgOA1;
    HMap(:,:,3)=imgR1;
    RGBMap2=hsv2rgb(HMap);
    imwrite(RGBMap2,[ 'slice' num2str(slice) '_' 'alpha_cal.tiff'],'compression','none');

    [X, Y, Z] = sph2cart(Theta_ObsLSQ, alpha1, mask);
    oct_vec_3d = cat(3, X, Y, Z);
    outimgname = ['slice' num2str(slice) '_3d'];
    imwrite(oct_vec_3d,[ outimgname '.tiff'],'compression','none');
%     vec3D(:,:,i,:)=oct_vec_3d;
end

resolution = [0.01,0.01,0.15];
MakeNii('Stacked_axis3D_RGB_anisodiff_filtered.nii',vec3D,resolution);
% %%
biref = MRIread('Stacked_biref.nii');
mask (mask>0.0015)=nan; mask(mask<0.00005)=nan; mask(~isnan(mask))=1;


alpha = pi/2-Psi_ObsLSQ;
alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;

size(mri.vol)
resolution = [0.01,0.01,0.15];
vec3D = permute (mri.vol,[2,1,3,4]);
MakeNii('Stacked_axis3D_RGB_despeckle2.nii',vec3D,resolution);

% for i = 1:19
% imwrite (squeeze(mri.vol(:,:,i,:)/255),['3d' num2str(i) '.tif'],'compression','none')
% end
% 
% for i = 1:19
% a = imread (['Stacked_axis3D_RGB_despeckle' sprintf('%02i', i) '.tif']);
% v3 (:,:,i,:)=a;
% end
% name = 'Stacked_axis3D_RGB_despeckle.nii';
% resolution = [0.01,0.01,0.15];
% MakeNii(name,v3,resolution);