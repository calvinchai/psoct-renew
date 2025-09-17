function [] = do_stacking(ParameterFile, modality)
% 
% [] = do_stacking(ParameterFile, modality)
%   stack and save specified 2D/en-face contrast into nifti volume
% 
% USAGE:
% 
%  INPUTS:          
%       ParameterFile     =   Path to Parameters.mat file
%                           - expects to load Parameters and Stack: 
% 
%   ~2024-03-06~  
% 

load(ParameterFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING SLICE INPUT 
sliceid     = Stack.sliceid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir       = Stack.indir;
outdir      = Stack.outdir; if ~exist(outdir,'dir'); mkdir(outdir); end

fprintf(' - Input directory = \n%s\n',indir);
fprintf(' - Output directory = \n%s\n',outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING RESOLUTION 
resolution = Stack.Resolution;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set Waitbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = waitbar(0,'1','Name',['Stacking ' modality ' ...'],...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

%% Stacking
for i = 1:length(sliceid)
    
    Imag = [indir filesep  modality '_slice' sprintf('%03i',sliceid(i)) '.mat'];
    load(Imag);
    S = whos('-file',Imag);
    I = eval(S.name);
    if ~isa(I,'single');single(I);end
    
    if i == 1; stacknii = zeros([size(I,[1,2]),length(sliceid)],'like',I) ;end %fprintf('Stacking %s ',modality)
    
    stacknii(:,:,i) = I;
    
    % Update waitbar and message
    step = i; steps = length(sliceid);
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    waitbar(step/steps,f,sprintf('Slice %i/ -%i- /%i',sliceid(1), sliceid(i), sliceid(end)))
end
stacknii = rot90(stacknii);


%% Saving

fprintf(' xy res=%g mm\n z  res=%g mm\n',resolution(1),resolution(3));
waitbar(step/steps,f,sprintf('Saving Stacked_%s.nii.gz', modality))

name = sprintf('%s/Stacked_%s.nii.gz', outdir, modality);
write_mri(stacknii, name, resolution, 'float')

if strcmp(modality,'Orientation')
    stacknii = get_rgb_4D(stacknii, [-90, 90]); 
    stacknii = stacknii*255;
    name = sprintf('%s/Stacked_%s.RGB.nii.gz', outdir, modality);
    write_mri(stacknii, name, resolution, 'int')
end

delete(f)

end

function I_rgb = get_rgb_4D(ori, range)
r1 = range(1);
r2 = range(2);

sz = size(ori,[1,2,3]);
ori = double(ori);
ori(ori>r2)=ori(ori>r2)-180;
ori(ori<r1)=ori(ori<r1)+180;
ori = reshape(ori,sz(1),[]);

I_hsv = mat2gray(ori,range);
I_hsv(:,:,2:3) = 1;
I_rgb=hsv2rgb(I_hsv);

I_rgb = reshape (I_rgb,sz(1),sz(2),sz(3),3);

end

function write_mri(I, name, res, datatype)
    if ~exist("datatype","var")
        datatype = 'float';
    end

    disp(' - making hdr...');
    % Make Nifti and header
    colres = res(2); 
    rowres = res(1); 
    sliceres = res(3); 
    % mri.vol = I;
    mri.volres = [res(1) res(2) res(3)];
    mri.xsize = rowres;
    mri.ysize = colres;
    mri.zsize = sliceres;
    a = diag([-colres rowres sliceres 1]);
    mri.vox2ras0 = a;
    mri.volsize = size(I,[1,2,3]);
    mri.vol = I;
    % mri.vol = flip(mri.vol,1);
    MRIwrite(mri,name,datatype);
    disp(' - done - ');
end