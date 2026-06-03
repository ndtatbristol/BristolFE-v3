clear all;
close all;

%ABOUT THIS EXAMPLE
%This is identical to ex_2d_basic.m except that the model is defined 
%parametrically in the separate function, mod_2d_basic.m. This is better
%practice for (a) keeping scripts tidy and (b) allowed automated parametric
%studies to be performed.

%Any of the default parameters (see top of mod_2d_basic.m for complete 
%list) can be overwritten here, e.g.
params.model_size = 10e-3;
params.el_typ_to_use_for_solid = 'CPE3'; 
params.els_per_wavelength = 15;

show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = mod_2d_basic(params);

%Show the mesh and stop if requested
if show_geom_only
    fig = figure;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
    display_options.node_sets_to_plot(2).col = 'g.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

%--------------------------------------------------------------------------
%RUN THE MODEL

%This is where the model actually gets executed.
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Plot history output at monitoring node
figure;
plot(steps{1}.load.time, res{1}.dsps);
xlabel('Time (s)')

%Animate field output result
figure;
h_patch = fn_show_geometry(mod, matls, el_types, display_options);
anim_options.fld_time = res{1}.fld_time;
fn_run_animation(h_patch, res{1}.fld, anim_options);
