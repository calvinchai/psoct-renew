%% Stacked_axis3D_RGB.nii
filebase = '/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/StackNii_normal';

Psi = MRIread([filebase filesep 'Stacked_Psi.nii']);
Theta = MRIread([filebase filesep 'Stacked_Theta.nii']);
% Ret = MRIread([filebase filesep 'Stacked_Retardance.nii']);
% RetRange = [15 30];
% mask = mat2gray(mask, [0.00025 0.0015 ])
%     mask (mask>0.0015)=nan; mask(mask<0.00025)=nan; mask(~isnan(mask))=1;
alpha = pi/2-Psi.vol;
alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;
[X, Y, Z] = sph2cart(Theta.vol, alpha, Theta.vol*0+1);
vec3D = cat(4, X, Y, Z);
RGB = vec3D*255;
figure; imshow(squeeze(vec3D(:,:,1,:)));

mri = Psi;
mri.vol = RGB;
MRIwrite(mri,[filebase filesep 'Stacked_axis3D_RGB.nii'])

%% Stacked_axis3D_RGB.Ret_mask.nii
filebase = '/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/StackNii_normal';

Psi = MRIread([filebase filesep 'Stacked_Psi.nii']);
Theta = MRIread([filebase filesep 'Stacked_Theta.nii']);
Ret = MRIread([filebase filesep 'Stacked_Retardance.nii']);
RetRange = [15 25];
% mask = mat2gray(mask, [0.00025 0.0015 ])
%     mask (mask>0.0015)=nan; mask(mask<0.00025)=nan; mask(~isnan(mask))=1;
alpha = pi/2-Psi.vol;
alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;
[X, Y, Z] = sph2cart(Theta.vol, alpha, mat2gray(Ret.vol,RetRange));
vec3D = cat(4, X, Y, Z);
RGB = vec3D*255;
figure; imshow(squeeze(vec3D(:,:,1,:)));

mri = Psi;
mri.vol = RGB;
MRIwrite(mri,[filebase filesep 'Stacked_axis3D_RGB.Ret_mask.nii'])

% [X, Y, Z] = sph2cart(Theta_ObsLSQ,Psi_ObsLSQ,mat2gray(MosaicFinal,[15 40]));
% vec3D = cat(3,X,Y,Z);
% 
% figure;histogram(X(:))
% figure; imshow(vec3D);

%% Stacked_axis3D_RGB.biref_mask.nii
filebase = '/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/StackNii_normal';

Psi = MRIread([filebase filesep 'Stacked_Psi.nii']);
Theta = MRIread([filebase filesep 'Stacked_Theta.nii']);
biref = MRIread([filebase filesep 'Stacked_biref.nii']);
% RetRange = [15 30];
mask = medfilt3(biref.vol, [5 5 1]);
mask (mask>0.0015)=nan; mask(mask<0.00025)=nan; mask(~isnan(mask))=1;

alpha = pi/2-Psi.vol;
alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;
num_iter = 50;option = 2;
delta_t = 1/100;kappa = 1;
for ii = 1:19
    OC = angle2tensor (alpha(:,:,ii)/pi*180); % input in degrees
    OC1 = anisodiff2D(OC(:,:,1),num_iter,delta_t,kappa,option);
    OC2 = anisodiff2D(OC(:,:,2),num_iter,delta_t,kappa,option);
    OC3 = anisodiff2D(OC(:,:,3),num_iter,delta_t,kappa,option);
    OC4 = anisodiff2D(OC(:,:,4),num_iter,delta_t,kappa,option);
    alpha1(:,:,ii) = tensor2angle (cat (3,OC1,OC2,OC3,OC4)); 
end
[X, Y, Z] = sph2cart(Theta.vol, alpha1, mask);
vec3D = cat(4, X, Y, Z);
RGB = vec3D*255;
figure; imshow(squeeze(vec3D(:,:,1,:)));

mri = Psi;
mri.vol = RGB;
MRIwrite(mri,[filebase filesep 'Stacked_axis3D_RGB_diffusion2.nii'])