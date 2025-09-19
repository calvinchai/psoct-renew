mosaicnum = 16;

for tilenum = [520]%[309,310] %[312,313,314,317,318,319,320,275,276,277,278,279,280,281,266,267]
    disp(tilenum)
    % file_path_ = '/autofs/cluster/connects2/users/Nate/I80premotor_focustest_mineraloil_05092025/';
    file_path_ = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/RawData/'
    filename = [file_path_ ,sprintf('mosaic_%03d_image_%04d_processed_cropped_focus.nii',mosaicnum, tilenum)];
    aa=niftiread(filename);
    s = size(aa,1)/4;
    Jones1 = aa(1:s,:,:) + 1j*aa(s+1:s*2,:,:);
    Jones2 = aa(s*2+1:s*3,:,:) + 1j*aa(s*3+1:s*4,:,:);
    
    % dBI3D = flip(10*log10( abs(j1).^2+abs(j2).^2 ),3);

    R3D=flip(atan(abs(Jones1)./abs(Jones2))/pi*180,3); % tissue top should be dark
    offset = 0;%100./180*pi;
    phase1 = (angle(Jones1));
    phase2 = (angle(Jones2));
    phi=(phase1-phase2)+offset*2;
    index1=phi>pi;
    phi(index1)=phi(index1)-2*pi;
    index2=phi<-pi;
    phi(index2)=phi(index2)+2*pi;
    O3D = flip( (phi)/2/pi*180 ,3);
    
    % savefilename = sprintf('mosaic_%03d_image_%03d_processed_dBI.nii',mosaicnum,tilenum);
    % niftiwrite(dBI3D,[file_path_,savefilename])
    save(['/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/tempProcessed/',sprintf('mosaic_%03d_image_%04d_processed_O3D.mat',mosaicnum, tilenum)],'O3D')
    save(['/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/tempProcessed/',sprintf('mosaic_%03d_image_%04d_processed_R3D.mat',mosaicnum, tilenum)],'R3D')
end


