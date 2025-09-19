% load(['/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/vol_recon_proj_JW_2025/par_' filetag '_data.mat'])

alpha = pi/2-Psi_ObsLSQ;
alpha(alpha>pi/2)=alpha(alpha>pi/2)-pi;
alpha(alpha<-pi/2)=alpha(alpha<-pi/2)+pi;

I_B = (biref_ObsLSQ-prctile(biref_ObsLSQ(:),1))/(prctile(biref_ObsLSQ(:),20)) ;
I_B(I_B>1)=1;
I_B(I_B<0)=0;

[X, Y, Z] = sph2cart(Theta_ObsLSQ, alpha, I_B);

oct_vec_3d = cat(3, Y, Z, X);

figure; imshow(abs(oct_vec_3d));