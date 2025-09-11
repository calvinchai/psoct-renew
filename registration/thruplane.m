clear;
addpath ('/autofs/cluster/octdata2/users/Chao/code/demon_registration_version_8f');
addpath('/autofs/cluster/octdata2/users/Chao/code/telesto');
addpath ('/space/omega/1/users/3d_axis/PAPER/scripts');
parpool_num = 10;

basename1 = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/Preliminary_stitching_noComputeOverlap/';
basename2 = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/Preliminary_stitching_noComputeOverlap/';

disp('start parpool');
poolobj = parpool(parpool_num);
disp('parpool started');

for slice = 8
    % % read first set
    slice_id1 = slice*2-1;
    slice_id2 = slice*2;
    
    fixed1 = imread([basename1, 'mosaic_',sprintf('%03d',slice_id1),'_ori_stokes.tiff']); % orientation fixed
    moving1 = imread([basename2, 'mosaic_',sprintf('%03d',slice_id2),'_ori_stokes.tiff']);  % orientation moving

    fixed2 = imread([basename1,'mosaic_',sprintf('%03d',slice_id1),'_biref_stokes.tiff']);
    moving2 = imread([basename2, 'mosaic_',sprintf('%03d',slice_id2),'_biref_stokes.tiff']);

    fixed2 = imgaussfilt(fixed2,3);
    moving2 = imgaussfilt(moving2,3);
    
    clear moving1_temp moving2_temp
    f1_dims = size(fixed1);
    m1_dims = size(moving1);

    f1_xrange1 = 1;             
    f1_xrange2 = f1_dims(1)-39; 

    f1_yrange1 = 60;                  
    f1_yrange2 = m1_dims(2)+59; 


    m1_xrange1 = 40;          
    m1_xrange2 = m1_dims(1); 

    m1_yrange1 = 1;        
    m1_yrange2 = m1_dims(2); 

    fixed_o1 =  -fixed1(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
    fixed_bi1 =  fixed2(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
%
    moving_o1 = -moving1(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);
    moving_bi1 = moving2(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);
    %   figure, subplot(1,2,1);imagesc(fixed_bi1);subplot(1,2,2);imagesc(moving_bi1)

%     fixed_o1 = imresize (fixed_o1, [size(fixed_bi1,1) size(fixed_bi1,2)],'nearest');
%     moving_o1 = imresize (moving_o1, [size(moving_bi1,1) size(moving_bi1,2)],'nearest');
    




    thruplane_reg_optiz_tensor_XY_final_par_JW (fixed_bi1,moving_bi1, fixed_o1, moving_o1, -10, ['par_slice' num2str(slice)],parpool_num)

clear fixed_bi1 moving_bi1 fixed_o1 moving_o1
delete(poolobj)
end

% rest is done in thruplane_reg_optiz_tensor_XY_final_par_JW

