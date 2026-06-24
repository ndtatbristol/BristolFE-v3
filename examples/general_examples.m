clear all;
close all;

%ABOUT THIS SCRIPT
%This runs any of the example models and shows results

%Uncomment one of the following model file names to determine which one
%will be used for comparison:
model_to_run = @mod_2d_basic;
% model_to_run = @mod_3d_basic;
model_to_run = @mod_2d_advanced;
% model_to_run = @mod_3d_advanced;
% model_to_run = @mod_3d_islaa;
model_to_run = @mod_2d_oblique;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 8;

%If you just want to see the model (without running it, set show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

if size(mod.nds, 2) == 3
    %Currently there is no plotting option for 3D field outputs so turn
    %them off for 3D models
    fe_options.field_output_every_n_frames = inf;
end

%Show the mesh and stop if requested
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'g.';
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

%--------------------------------------------------------------------------
%RUN THE MODEL

%This is where the model actually gets executed.
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Plot summed history output over monitoring nodes
figure;
if isfield(steps{1}.load, 'wts')
    ascan = steps{1}.load.wts .' * res{1}.dsps;
else
    ascan = sum(res{1}.dsps, 1);
end
plot(steps{1}.load.time, ascan);
xlabel('Time (s)')

%Animate field output result if available
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    anim_options.fld_time = res{1}.fld_time;
    % anim_options.norm_val = 20;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end