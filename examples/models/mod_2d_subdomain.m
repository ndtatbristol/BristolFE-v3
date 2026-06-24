function [main, fe_options, params] = mod_2d_subdomain(params)

%This uses mod_2d_advanced (without either of the scatterers) as the 
%pristine model. The crack that would normally defined in mod_2d_advanced
%is then put in a sub-domain model.

default_params.fe_options.field_output_every_n_frames = 20;
%--------------------------------------------------------------------------
if isfield(params, 'fe_options') && isfield(default_params, 'fe_options')
    params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
else
    default_params.fe_options = [];
end
params = fn_set_default_fields(params, default_params);

params.include_crack = 0;
params.include_scatterer = 0;

[main.mod, main.matls, main.el_types, steps, fe_options, params] = mod_2d_advanced(params);
main.trans{1}.nds = steps{1}.load.frc_nds;
main.trans{1}.dfs = steps{1}.load.frc_dfs;

%Create a subdomain in the middle with a hole in surface as scatterer
scatterer_centre = (max(params.crack_vtcs) + min(params.crack_vtcs)) / 2;
min_dims = max(params.crack_vtcs) - min(params.crack_vtcs);
inner_bdry = [-1,-1;-1,1;1,1;1,-1] / 2 * (max(min_dims) + 2 * params.el_size) + scatterer_centre;

%First subdomain left empty to show what happens
main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);

%Second subdomain contains crack
main.doms{2}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);
cod = params.el_size / 10;
main.doms{2}.mod = fn_2d_add_crack(main.doms{2}.mod, main.el_types, params.crack_vtcs, [], cod);

%Input signal
main.inp.time = steps{1}.load.time;
main.inp.sig = steps{1}.load.frcs;
end

% %Show the mesh
% if show_geom_only %suppress graphics when running all scripts for testing
%     figure;
%     display_options.draw_elements = 0;
%     display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
%     display_options.node_sets_to_plot(1).col = 'r.';
%     h_patch = fn_show_geometry_with_subdomains(main, display_options);
%     drawnow
%     return
% end
% %--------------------------------------------------------------------------
% 
% %Run main model
% main = fn_run_main_model(main, fe_options);
% 
% %Demonstration of how sub-domain model can be run for multiple random scatterers
% figure;
% results = zeros(numel(main.inp.time), no_scatterers);
% for s = 1:no_scatterers
%     scat_pts =   fn_2d_create_smooth_random_blob(0.4, 3, 360) * scatterer_size / 2 + scatterer_centre;
%     main.doms{1}.mod = fn_2d_add_inclusion_or_void(empty_subdomain, main.el_types, scat_pts, 0, 0);
%     main = fn_run_subdomain_model(main, fe_options);
%     results(:,s) = sum(main.doms{1}.res.fmc.time_data, 2);
%     subplot(2, no_scatterers, s);
%     plot(main.doms{1}.res.fmc.time, real(results(:,s)))
%     subplot(2, no_scatterers, no_scatterers + s);
%     fn_show_geometry(main.doms{1}.mod, main.matls, main.el_types, []);
% end
% 
% %Animate results if requested
% if ~isinf(fe_options.field_output_every_n_frames)
%     figure;
%     anim_options.repeat_n_times = 1;
%     anim_options.db_range = [-40, 0];
%     anim_options.pause_value = 0.001;
%     h_patches = fn_show_geometry_with_subdomains(main, anim_options);
%     fn_run_subdomain_animations(main, h_patches, anim_options);
% end
% 
% %Run validation model if requested (just does it for last scatterer)
% if run_validation_models
%     fe_options.validation_mode = 1;
%     main = fn_run_main_model(main, fe_options);
% 
%     %View the time domain data and compare wih validation
%     figure;
%     i = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
%     mv = max(abs(sum(main.doms{1}.res.fmc.time_data(i:end,: ), 2)));
%     plot(main.doms{1}.res.fmc.time, real(sum(main.doms{1}.res.fmc.time_data, 2)) / mv, 'k', 'LineWidth', 2);
%     hold on;
%     plot(main.doms{1}.val.fmc.time, real(sum(main.doms{1}.val.fmc.time_data, 2)) / mv, 'g:', 'LineWidth', 2);
%     plot(main.res.fmc.time, real(sum(main.res.fmc.time_data, 2)) / mv, 'b');
%     ylim([-1,1]);
%     yyaxis right
%     plot(main.doms{1}.res.fmc.time, 20 * log10(abs(sum(main.doms{1}.res.fmc.time_data, 2) - sum(main.doms{1}.val.fmc.time_data, 2)) / mv));
%     ylim([-60, 0]);
%     legend('Sub-domain method', 'Validation', 'Pristine', 'Difference (dB)');
% 
%     %Animate validation results if requested
%     if ~isinf(fe_options.field_output_every_n_frames)
%         %Animate result
%         figure;
%         anim_options.repeat_n_times = 1;
%         anim_options.db_range = [-40, 0];
%         anim_options.pause_value = 0.001;
%         h_patches = fn_show_geometry(main.doms{1}.val_mod, main.matls, main.el_types, anim_options);
%         anim_options.fld_time = main.res.trans{1}.fld_time;
%         fn_run_animation(h_patches, main.doms{1}.val.trans{1}.fld, anim_options);
%     end
% end
