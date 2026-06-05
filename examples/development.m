clear all;
close all;

%ABOUT THIS SCRIPT
%Development script!

%Uncomment one of the following model file names to determine which one
%will be used for comparison:
% model_to_run = @mod_2d_basic;
model_to_run = @mod_3d_basic;
model_to_run = @mod_2d_advanced;
% model_to_run = @mod_3d_advanced;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

pogo_path = 'C:\Program Files\Pogo\windows\new version';
pogo_matlab_path = 'C:\Program Files\Pogo\matlab';


%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 4;%13 is OK (775k els); 14 is out-of-memory (932k elements)
params.include_fluid_region = 1;
params.include_absorbing_boundary = 1;
params.include_crack = 1;
params.include_scatterer = 1;
params.scatterer_is_void = 1;
params.fe_options.field_output_every_n_frames = inf;



%If you just want to see the model (without running it, set 
%show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);
% i = mod.nds(:, 3) > 6e-3;
% mod.nds(i, 3) = (mod.nds(i, 3) - 6e-3) * 1.5 + 6e-3;

%Show the mesh and stop if requested
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'g.';
display_options.draw_elements = 1;
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

figure;
col = {'k', 'r--'};
ver = {'v6', 'v5'};
% ver = {'v6'};
for v = 1:numel(ver)
    fe_options.matrix_builder_version = ver{v};
    res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
    %Plot summed history output over monitoring nodes
    plot(steps{1}.load.time, sum(res{1}.dsps, 1), col{v});
    hold on;
end
xlabel('Time (s)')
