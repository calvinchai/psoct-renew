function [dBI3D,R3D,O3D] = view_dBI3D_from_complex(mosaic_num,tile_num,data_path)

filename = sprintf('mosaic_%03i_image_%04i_processed_cropped_focus.nii',mosaic_num,tile_num);
complex3D = niftiread(sprintf('%s/%s',data_path,filename));

if(mod(mosaic_num,2) == 1) % Odd mosaic nums --> Normal incidence
    num_x_px = 350;
elseif(mod(mosaic_num,2) == 0) % Even mosaic nums --> Tilted incidence
    num_x_px = 200;
end
Jones1_real = complex3D(0*num_x_px + 1:1*num_x_px,:,:);
Jones1_imag = complex3D(1*num_x_px + 1:2*num_x_px,:,:);
Jones2_real = complex3D(2*num_x_px + 1:3*num_x_px,:,:);
Jones2_imag = complex3D(3*num_x_px + 1:4*num_x_px,:,:);
Jones1 = Jones1_real + sqrt(-1)*Jones1_imag;
Jones2 = Jones2_real + sqrt(-1)*Jones2_imag;

IJones = abs(Jones1).^2+abs(Jones2).^2;
dBI3D = flip(10*log10(IJones),3);

R3D = flip(atan(abs(Jones1)./abs(Jones2))/pi*180,3);

offset = 100/180*pi;
phase1 = angle(Jones1);
phase2 = angle(Jones2);
phi = (phase1 - phase2) + offset*2;
index1 = (phi > pi);
phi(index1) = phi(index1) - 2*pi;
index2 = (phi < -pi);
phi(index2) = phi(index2) + 2*pi;
O3D = flip(phi/2/pi*180,3);
end