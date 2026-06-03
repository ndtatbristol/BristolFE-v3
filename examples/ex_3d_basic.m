clear all
close all;

show_geom_only = 0;

%ABOUT THIS EXAMPLE
%This example is designed to show a simple 3D model executed using the BristoFE solver
%and Pogo

%PARAMETRIC DESCRIPTION OF MODEL
%The model is described in terms of a small number of parameters

params = [];

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = mod_3d_basic(params);

%Show the mesh and stop if requested
if show_geom_only 
    figure;
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
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
plot(steps{1}.load.time, -sum(res{1}.dsps, 1));
xlabel('Time (s)')

