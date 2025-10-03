addpath ('/local_mount/space/megaera/1/users/kchai/code/psoct-data-processing/registration/');
source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)

source = '/vast/fiber/projects/20250920_CCtest3_15degrees_CCW/';
basename = '/homes/5/kc1708/project/EXP3/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)
stitch_section(basename, [9 5], [13 5])
thruplane(basename, gamma)
RGB_3Daxis(basename)

source = '/vast/fiber/projects/20250920_CCtest4_30degrees_CW/';
basename = '/homes/5/kc1708/project/EXP4/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)
stitch_section(basename, [13 5], [9 5])
thruplane(basename, gamma)
RGB_3Daxis(basename)

source = '/vast/fiber/projects/20250920_3d_axis_3_illuminations/';
basename = '/homes/5/kc1708/project/EXP5/';
gamma = 10;
batch_process_cc(source, [basename 'processed/'], 100)

stitch_section(basename, [6 5], [9 5], 'mosaic_001', 'mosaic_003')
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250920_3d_axis_3_illuminations/';
basename = '/homes/5/kc1708/project/EXP5-2/';
gamma = 25;
batch_process_cc(source, [basename 'processed/'], 100)

stitch_section(basename, [11 5], [8 5], 'mosaic_002', 'mosaic_003')
thruplane(basename, gamma)
RGB_3Daxis(basename)



source = '/vast/fiber/projects/20250919_megatome_surface_testing_X/';
basename = '/homes/5/kc1708/project/EXP6/';
gamma = 10;
batch_process_cc(source, [basename 'processed/'], 100, "","")

stitch_section(basename, [19 2], [30 2] )
thruplane(basename, gamma)
RGB_3Daxis(basename)

source = '/vast/fiber/projects/20250920_megatome_surface_testing_Y/';
basename = '/homes/5/kc1708/project/EXP7/';
gamma = 10;
batch_process_cc(source, [basename 'processed/'], 100, "","")

stitch_section(basename, [2 23], [2 23] )
thruplane(basename, gamma)
RGB_3Daxis(basename)

%% New Birefringence
source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-NewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100, "", "new")

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-NewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100,"", "new")
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)

%% New Orientation

source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-NewOri/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100, "new", "")

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-NewOri/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100,"new", "")
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)
%% New Orientation and Birefringence

source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-NewOriNewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100, "new", "new")

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-NewOriNewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100,"new", "new")
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)

%% UnWrap

source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-UnWrap/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100, "", "", true)

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-UnWrap/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100,"", "", true)
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


%% UnWrap Both New

source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-UnWrap-NewOri-NewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100, "new", "new", true)

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-UnWrap-NewOri-NewBiref/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100,"new", "new", true)
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)



%% 
addpath ('/local_mount/space/megaera/1/users/kchai/code/psoct-data-processing/registration/');
source = '/vast/fiber/projects/20250919_CCtest/';
basename = '/homes/5/kc1708/project/EXP1-100-200/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)

stitch_section(basename, [5 6], [8 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)


source = '/vast/fiber/projects/20250919_CCtest2_15degree/';
basename = '/homes/5/kc1708/project/EXP2-100-200/';
gamma = -15;
batch_process_cc(source, [basename 'processed/'], 100)
stitch_section(basename, [7 6], [5 6])
thruplane(basename, gamma)
RGB_3Daxis(basename)
