% Parameters
save_path = '/autofs/cluster/connects2/users/data/I80_premotor_slab_2025_05_13/ProcessedData/StitchingFiji/Stitching_diagram';
img_size = [350, 350];
num_tiles = 21*48;

for n = 1:num_tiles
    number_str = sprintf('%d',n);
    image = create_tile_with_number(img_size,number_str);
    imwrite(image, [save_path,sprintf('/number_img_%03d.tiff',n)]);
end

% figure(1);
% imshow(image, []);


%%

function img_gray = create_tile_with_number(img_size,number_str)

font_size = 60;

% Create a figure
f = figure('Visible', 'off', 'Position', [100 100 img_size(2) img_size(1)]);
axes('Position', [0 0 1 1]);  % Fill entire figure
axis off
xlim([0 img_size(2)])
ylim([0 img_size(1)])
set(gca, 'YDir', 'reverse')  % So text aligns top-left like an image

% Add centered text
text(img_size(2)/2, img_size(1)/2, number_str, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', font_size, ...
    'FontWeight', 'bold', ...
    'Color', 'black');

% Capture frame
frame = getframe(gca);
img_RGB = frame.cdata;

% Convert to grayscale
img_gray = rgb2gray(img_RGB);

% Resize if necessary (getframe may not exactly match 350x350)
img_gray = imresize(img_gray, img_size);




end