
function thruplane(basename, gamma)
%gamma = -10

addpath ('/autofs/cluster/octdata2/users/Chao/code/demon_registration_version_8f');
addpath('/autofs/cluster/octdata2/users/Chao/code/telesto');
addpath ('/space/omega/1/users/3d_axis/PAPER/scripts');
parpool_num = 24;

basename1 = basename;
basename2 = basename1;
slice_numbers = [10,11];

disp('start parpool');
% poolobj = parpool(parpool_num);
disp('parpool started');
for slice = slice_numbers
    % % read first set
    slice_id1 = slice*2-1;
    slice_id2 = slice*2;
    
    fixed1 = imread([basename1, 'mosaic_',sprintf('%03d',slice_id1),'_ori.tif']); % orientation fixed
    moving1 = imread([basename2, 'mosaic_',sprintf('%03d',slice_id2),'_ori.tif']);  % orientation moving

    fixed2 = imread([basename1,'mosaic_',sprintf('%03d',slice_id1),'_biref.tif']);
    moving2 = imread([basename2, 'mosaic_',sprintf('%03d',slice_id2),'_biref.tif']);
    
    fixed2 = imgaussfilt(fixed2,3);
    moving2 = imgaussfilt(moving2,3);
    
    clear moving1_temp moving2_temp
    f1_dims = size(fixed1);
    m1_dims = size(moving1);
    
%     aip1 = imread([basename1, 'mosaic_', sprintf('%03d', slice_id1), '_aip.tif']);
% aip2 = imread([basename2, 'mosaic_', sprintf('%03d', slice_id2), '_aip.tif']);
% 
% mask1 = aip1 > 60;
% mask2 = aip2 > 60;
% applyMask = @(img, m) ...
%     (ndims(img)==3) .* bsxfun(@times, img, cast(repmat(m,[1 1 size(img,3)]), class(img))) + ...
%     (ndims(img)~=3) .* (img .* cast(m, class(img)));
% 
% % --- apply masks ---
% fixed1  = applyMask(fixed1,  mask1);
% fixed2  = applyMask(fixed2,  mask1);   % same mask as fixed1 (same basename/slice)
% moving1 = applyMask(moving1, mask2);
% moving2 = applyMask(moving2, mask2);
    
%     f1_xrange1 = 1;             
%     f1_xrange2 = end; 
% 
%     f1_yrange1 = 1;                  
%     f1_yrange2 = end; 
% 
% 
%     m1_xrange1 = 1;          
%     m1_xrange2 = end; 
% 
%     m1_yrange1 = 1;        
%     m1_yrange2 = end; 
% 
%     fixed_o1 =  -fixed1(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
%     fixed_bi1 =  fixed2(f1_xrange1:f1_xrange2,f1_yrange1:f1_yrange2);
% %
%     moving_o1 = -moving1(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);
%     moving_bi1 = moving2(m1_xrange1:m1_xrange2,m1_yrange1:m1_yrange2);

    fixed_o1 =  -fixed1;
    fixed_bi1 =  fixed2;
%
    moving_o1 = -moving1;
    moving_bi1 = moving2;
    %   figure, subplot(1,2,1);imagesc(fixed_bi1);subplot(1,2,2);imagesc(moving_bi1)

%     fixed_o1 = imresize (fixed_o1, [size(fixed_bi1,1) size(fixed_bi1,2)],'nearest');
%     moving_o1 = imresize (moving_o1, [size(moving_bi1,1) size(moving_bi1,2)],'nearest');
    




    thruplane_reg_optiz_tensor_XY_final_par(fixed_bi1,moving_bi1, fixed_o1, moving_o1, gamma, [basename 'par_slice' num2str(slice)],parpool_num)

clear fixed_bi1 moving_bi1 fixed_o1 moving_o1
% delete(poolobj)
end

% rest is done in thruplane_reg_optiz_tensor_XY_final_par_JW

end