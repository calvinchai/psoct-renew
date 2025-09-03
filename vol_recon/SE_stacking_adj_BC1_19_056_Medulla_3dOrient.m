            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load('/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned//imregcorr_p.normal.mat','tran_idx','crop_idx','slice_idx')

load('/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned//imregcorr_p.all.mat','tran_idx','crop_idx','slice_idx')
%%
thru_planeRes_stack = 0.15;
in_planeRes = 0.01;
tran_idx
crop_idx
slice_idx

run1 = 1:3; run2 = 1:3; run3 = 1:13;
para.ProcDir = {'/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned'};


para.in = [run1 run2 run3 ];                                                      %sliceid from multiple runs
para.out = 1:length(para.in);             %
para.run = [ones(size(run1))*1 ones(size(run2))*2 ones(size(run3))*3];
sliceidx = [1:19;1:19;(1:19)*0+1];

outdir = ['/autofs/cluster/octdata3/users/BC1_19-056_Medulla/processed_cleaned/StackNii_normal']; %%%% output run2 folder
mkdir(outdir);

 %col_px = sz(2);  %

Modality_2D = {'Psi','Theta','biref'};
for m = 3
    modality = Modality_2D{m}; 
    fprintf('Stacking %s ',modality)
    Imag = [para.ProcDir{1} filesep sprintf('slice%i_data',19) '.mat'];
    load(Imag);
    S = whos('-file',Imag);
    I = eval([modality '_ObsLSQ']);
    I = rot90(I,-1);
    sz = size(I); %% xy size for stacked nifti
    stacknii = zeros([sz,length(sliceidx)],'single');
    
%% Stacking
for i = para.out%1:length(sliceid)
%     Info = sliceidx(:,i); %[in out run]
    i_in = sliceidx(1,i);
    i_out = sliceidx(2,i);
    i_run = sliceidx(3,i);
    
    ProcDir = para.ProcDir{i_run}; 
    indir = [ProcDir];
    
    fprintf('\tslice %d\n',i);
    
    Imag = [indir filesep sprintf('slice%i_data',i_in) '.mat'];
    load(Imag,[modality '_ObsLSQ']);
%     S = whos('-file',Imag);
    I = eval([modality '_ObsLSQ']);
    I = rot90(I,-1);
    sz = size(I);

        
    c1s=crop_idx(slice_idx(i),1);
    c1e=crop_idx(slice_idx(i),2);
    c2s=crop_idx(slice_idx(i),3);
    c2e=crop_idx(slice_idx(i),4);
    t1s=tran_idx(slice_idx(i),1)-1+c1s;
    t1e=tran_idx(slice_idx(i),1)-1+c1e;
    t2s=tran_idx(slice_idx(i),2)-1+c2s;
    t2e=tran_idx(slice_idx(i),2)-1+c2e;
    
    stacknii(t1s:t1e,t2s:t2e,i) = I(c1s:c1e,c2s:c2e);
    
end
%stacknii = rot90(stacknii); % commented out due to the need to imregcorr in dBI first
%(adj) and syn with dBI Z stitching
%% Saving

% waitbar(step/steps,f,sprintf('Saving Stacked %s.nii', modality))

name = [outdir filesep sprintf('Stacked_%s.nii', modality)];

fprintf(' xy res=%g mm\n z  res=%g mm\n',in_planeRes,thru_planeRes_stack);
resolution      = [in_planeRes in_planeRes thru_planeRes_stack];


disp(' - making hdr...');
% Make Nifti and header
colres = resolution(2); 
rowres = resolution(1); 
sliceres = resolution(3); 
% mri.vol = I;
mri.volres = [resolution(1) resolution(2) resolution(3)];
mri.xsize = rowres;
mri.ysize = colres;
mri.zsize = sliceres;
a = diag([-colres rowres sliceres 1]);

mri.vox2ras0 = a;
mri.volsize = size(stacknii);
mri.vol = stacknii;
% mri.vol = flip(mri.vol,1);
MRIwrite(mri,name,'float');
disp(' - END - ');
end