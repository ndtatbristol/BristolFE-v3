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
% model_to_run = @mod_3d_fastener_hole;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

pogo_path = 'C:\Program Files\Pogo\windows\new version';
pogo_matlab_path = 'C:\Program Files\Pogo\matlab';


%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 12;16;%13 is OK (775k els); 14 is out-of-memory (932k elements) with v4; %15 (1.16M elements) still works with v6
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
[mod, matls, el_types, steps, tmp_fe_options, params] = model_to_run(params);
fe_options{1} = tmp_fe_options;
fe_options{1}.solver_mode = 'predictor corrector';
% fe_options{2} = tmp_fe_options;
% fe_options{2}.solver_mode = 'implicit';


% i = mod.nds(:, 3) > 6e-3;
% mod.nds(i, 3) = (mod.nds(i, 3) - 6e-3) * 1.5 + 6e-3;

%Show the mesh and stop if requested
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'g.';
display_options.draw_elements = 0;
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

figure;
col = {'k', 'r--'};
ver = {'v7', 'v6'};
ver = {'v7'};
k = 1;
for v = 1:numel(ver)
    for o = 1:numel(fe_options);
        fe_options{o}.matrix_builder_version = ver{v};
        fe_options{o}.dynamic_solver_version = 'v7';
        res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options{o});
        %Plot summed history output over monitoring nodes
        plot(steps{1}.load.time, sum(res{1}.dsps, 1), col{k});
        k = k + 1;
        hold on;
    end
end
xlabel('Time (s)')
