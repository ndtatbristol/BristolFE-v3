clear all;
close all;

%ABOUT THIS SCRIPT
%This runs any of the subdomain example models and shows results.

%Uncomment one of the following model file names to determine which one
model_to_run = @mod_2d_subdomain;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 8;
params.fe_options.field_output_every_n_frames = inf;
params.fe_options.solver_mode = 'imp';

%This will also run full domain validation models to compare to subdomain
%results if set
run_validation_models = 1;

%If you just want to see the model (without running it, set show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'subdoms']))

%Add models subfolder to path
addpath([fileparts(mfilename('fullpath')), filesep, 'models']);

%Define the model
[main, fe_options, params] = model_to_run(params);

if size(main.mod.nds, 2) == 3
    %Currently there is no plotting option for 3D field outputs so turn
    %them off for 3D models
    fe_options.field_output_every_n_frames = inf;
end

%Show the mesh and stop if requested
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
display_options.node_sets_to_plot(1).col = 'g.';
if show_geom_only
    figure;
    h_patch = fn_show_geometry_with_subdomains(main, display_options);
    drawnow
    return
end

%Run main model
main = fn_run_main_model(main, fe_options);

%Run sub-domain model(s)
main = fn_run_subdomain_model(main, fe_options);

%Run validation models if requested
if run_validation_models
    fe_options.validation_mode = 1;
    main = fn_run_main_model(main, fe_options);
end

%Show results
figure;
for i = 1:numel(main.doms)
    subplot(numel(main.doms), 1, i)
    j = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
    mv = max(abs(sum(main.doms{i}.res.fmc.time_data(j:end,: ), 2)));
    plot(main.res.fmc.time, real(sum(main.res.fmc.time_data, 2)) / mv, 'Color', [1,1,1] * 0.75);
    hold on;
    plot(main.doms{i}.res.fmc.time, real(sum(main.doms{i}.res.fmc.time_data, 2)) / mv, 'k', 'LineWidth', 2);
    if run_validation_models
        plot(main.doms{i}.val.fmc.time, real(sum(main.doms{i}.val.fmc.time_data, 2)) / mv, 'g:', 'LineWidth', 2);
        ylim([-1,1]);
        yyaxis right
        plot(main.doms{i}.res.fmc.time, 20 * log10(abs(sum(main.doms{i}.res.fmc.time_data, 2) - sum(main.doms{i}.val.fmc.time_data, 2)) / mv));
        ylim([-60, 0]);
        legend('Pristine', 'Sub-domain method', 'Validation', 'Difference (dB)');
    else
        ylim([-1,1]);
        legend('Pristine', 'Sub-domain method');
    end
end

%Animate results if requested
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    anim_options.repeat_n_times = 1;
    anim_options.db_range = [-40, 0];
    anim_options.pause_value = 0.001;
    h_patches = fn_show_geometry_with_subdomains(main, anim_options);
    fn_run_subdomain_animations(main, h_patches, anim_options);
end

