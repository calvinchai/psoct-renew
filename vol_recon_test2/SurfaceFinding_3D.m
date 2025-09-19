function [f_Surface3D] = SurfaceFinding_3D(I,istart,refine_range)
tic

% istart = 10;
fprintf('Finding surface start from z = %i.\n',istart)
% refine_range = 70;
fprintf('Set surface finding range to %i at the refinement step.\n',refine_range)
sz_medfilt = 40;

h = fspecial3('average',10); %% in PBS data, 10 smooth noise to tissue

f = imfilter(I,h);
df = diff(f,[],3);
[~ , ss] = max(df(:,:,istart:end),[],3);
f_Surface = ss+istart-1;

t = round(median(f_Surface(:)));         % get median of SF_simple result
if t>700; fprintf('Tile surface median is %i. Hardcode this to 100.\n%s \n', t,datetime); t = 100; end
[~ , ss2] = max(df(:,:,istart:t+refine_range),[],3);    % get index of max diff again but above depth t+refine_range

f_Surface3D = round(medfilt1(medfilt1(ss2,sz_medfilt,[],1,'truncate'),sz_medfilt,[],2,'truncate'))+istart-1; % 1D median filter the result

toc
end