clear all;
close all;

%ABOUT THIS EXAMPLE
%This example demonstrates multiple features in a 2D model, including
%fluid-solid coupling, absorbing regions, defect generation (crack, 
%inclusion or void). All of these can be turned on or off as required.

params.include_fluid_region = 1;
params.include_absorbing_boundary = 1;
params.include_crack = 1;
params.include_scatterer = 1;
params.els_per_wavelength = 10;
params.element_shape = 'tri';
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = mod_2d_advanced(params)

%Show the mesh and stop if requested
if show_geom_only 
    figure; 
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    drawnow
    return
end

%--------------------------------------------------------------------------
%RUN THE MODEL

res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    display_options.draw_elements = 0; %makes it easier to see waves if element edges not drawn
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    anim_options.repeat_n_times = 1;
    anim_options.fld_time = res{1}.fld_time;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end
