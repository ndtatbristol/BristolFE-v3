clear all;
close all;

rng(1);

%Test to isolate problems when subdom cuts free boundary (observed in 3D)

%Overall model geometry
xsize = 10e-3;
ysize = 10e-3;
clearance = 1e-3;
abs_bdry_thickness_in_wavelengths = 1;
els_per_wavelength = 5;
inc_angle_degs = 45; %wedge geometry will be calculated from this and inc_angle mode
inc_mode = 'S';

%Subdomain geometry
subdomain_size = 3e-3;

%Solid material properties
plate_matl_name = 'steel';
plate_matl_longitudinal_velocity = 5636;
plate_matl_shear_velocity = 3012;
plate_matl_density = 8900;

%Element types to use
el_typ_to_use_for_solid = 'CPE3';
% el_typ_to_use_for_solid = 'CPE4';

%Details of input signal
centre_freq = 5e6;
no_cycles = 5;

%Transducer details
trans_diam = 5e-3;

show_geom_only = 0; %Set to 1 to just show geometry without running model

run_validation_models = 1;
fe_options.field_output_every_n_frames = inf;1; %set to inf to suppress animations

fe_options.pogo_path = 'C:\Program Files\Pogo\windows\new version';
fe_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';
% fe_options.solver = 'pogo';
fe_options.sort_nds = 1;


%--------------------------------------------------------------------------
%END OF INPUTS

%Define the materials
%A cell array with an entry for each material used in the model is required.
plate_matl_i = 1;
main.matls{plate_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(plate_matl_name, plate_matl_longitudinal_velocity, plate_matl_shear_velocity, plate_matl_density);

%Define thickness and start of absorbing boundary region
max_wavelength = fn_estimate_max_min_wavelengths(main.matls, centre_freq);
abs_bdry_thickness = abs_bdry_thickness_in_wavelengths * max_wavelength;

%Maximum time
max_time = 1.5 * 2 * (ysize / plate_matl_longitudinal_velocity);

%Add path to BristolFE functions in case not already on path
addpath(genpath('../code'));
addpath(genpath('../subdoms'));
rmpath(genpath('../code/deprecated'));


%Define shape of model
xmin =0; xmax = xsize;
ymin =0; ymax = ysize;
bdry_pts = [
    xmin, 0
    xmax, 0
    xmax, ymax
    xmin, ymax];

%absorbing boundary for plate
abs_bdry_pts_plate = [
    xmin + abs_bdry_thickness, 0
    xmin + abs_bdry_thickness, ymax
    xmax - abs_bdry_thickness, ymax
    xmax - abs_bdry_thickness, 0];

%Other stuff
fe_options.dof_to_use = [1,2,3];

%Work out element size, create the nodes and elements of the mesh,
%determine elements in plate and wedge based on y value
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);
switch el_typ_to_use_for_solid
    case 'CPE3'
        main.mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
    case 'CPE4'
        main.mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size);
end

main.el_types = fn_2d_el_types();

main.mod.el_typ_i(:) = find(strcmp(main.el_types, el_typ_to_use_for_solid));

%Set materials
main.mod.el_mat_i(:) = plate_matl_i;

%Define the absorbing layer
main.mod = fn_2d_add_absorbing_layer(main.mod, abs_bdry_pts_plate, abs_bdry_thickness);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, main.matls, centre_freq);


%Lose the elements that need to go
el_ctrs = fn_calc_element_centres(main.mod.nds, main.mod.els);
els_to_go = (el_ctrs(:,1) - el_ctrs(:,2)) > 0 & el_ctrs(:,2) < ysize * 0.5;
% els_to_go = el_ctrs(:,1) < xsize / 2 & el_ctrs(:,2) < ysize / 2;

main.mod.els(els_to_go, :) = [];
main.mod.el_mat_i(els_to_go) = [];
main.mod.el_abs_i(els_to_go) = [];
main.mod.el_typ_i(els_to_go) = [];
[main.mod.nds, main.mod.els, ~, ~] = fn_remove_unused_nodes(main.mod.nds, main.mod.els);

%Define transducer
trans_end_pts = [
    xsize/2-trans_diam/2, ysize
    xsize/2+trans_diam/2, ysize];
nds = fn_find_nodes_nearest_to_line(main.mod.nds, trans_end_pts(1, :), trans_end_pts(2, :), el_size / 2);
main.trans{1}.nds = nds;
main.trans{1}.dfs = ones(size(nds)) * 2;

%Input signal
time_step = fn_get_suitable_time_step(main.matls, el_size);
time_pts = ceil(max_time / time_step);
main.inp.time = [0:time_pts - 1] * time_step;
main.inp.sig = fn_gaussian_pulse(main.inp.time, centre_freq, no_cycles);

%Create a subdomain in the middle with a hole in surface as scatterer
inner_bdry = [
    -subdomain_size / 2, -subdomain_size / 2
    -subdomain_size / 2, subdomain_size / 2% - el_size / 4
    subdomain_size / 2,  subdomain_size / 2% - el_size / 4
    subdomain_size / 2, -subdomain_size / 2] + [xsize, ysize] / 2;

inner_bdry(1,2) = inner_bdry(1,2) - 4 * el_size;

main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, abs_bdry_thickness);

%Show the mesh
    figure;
    subplot(1,2,1);
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(main.mod, main.matls, main.el_types, display_options);
    cols = 'rgbm';
    for i = 1:4
        b = main.doms{1}.mod.main_nd_i(main.doms{1}.mod.bdry_lyrs == i);
        plot(main.mod.nds(b,1),main.mod.nds(b,2), [cols(i), '.'], 'MarkerSize', 24);
    end

    subplot(1,2,2);
    display_options.node_sets_to_plot = [];
    h_patch = fn_show_geometry(main.doms{1}.mod, main.matls, main.el_types, display_options);
    cols = 'rgbm';
    for i = 1:4
        b = main.doms{1}.mod.main_nd_i(main.doms{1}.mod.bdry_lyrs == i);
        plot(main.mod.nds(b,1),main.mod.nds(b,2), [cols(i), '.'], 'MarkerSize', 24);
    end
if show_geom_only %suppress graphics when running all scripts for testing
    return
end
%--------------------------------------------------------------------------

%Run main model
main = fn_run_main_model(main, fe_options);

%Subdomain model with no scatterers
main = fn_run_subdomain_model(main, fe_options);


%Animate results if requested
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    anim_options.repeat_n_times = 1;
    anim_options.db_range = [-40, 0];
    anim_options.pause_value = 0.001;
    h_patches = fn_show_geometry_with_subdomains(main, anim_options);
    fn_run_subdomain_animations(main, h_patches, anim_options);
end

%Run validation model if requested (just does it for last scatterer)
if run_validation_models
    fe_options.validation_mode = 1;
    main = fn_run_main_model(main, fe_options);

    %View the time domain data and compare wih validation
    ref_subdom = 1;
    i = 1;
    mv = max(abs(sum(main.doms{ref_subdom}.res.fmc.time_data(i:end,: ), 2)));
    for subdom = 1:numel(main.doms)
        figure;
        plot(main.doms{subdom}.res.fmc.time, real(sum(main.doms{subdom}.res.fmc.time_data, 2)) / mv, 'k', 'LineWidth', 2);
        hold on;
        plot(main.doms{subdom}.val.fmc.time, real(sum(main.doms{subdom}.val.fmc.time_data, 2)) / mv, 'g:', 'LineWidth', 2);
        plot(main.res.fmc.time, real(sum(main.res.fmc.time_data, 2)) / mv, 'b');
        ylim([-1,1]);
        yyaxis right
        plot(main.doms{subdom}.res.fmc.time, 20 * log10(abs(sum(main.doms{subdom}.res.fmc.time_data, 2) - sum(main.doms{subdom}.val.fmc.time_data, 2)) / mv));
        ylim([-200, 0]);
        legend('Sub-domain method', 'Validation', 'Pristine', 'Difference (dB)');
    end

    if ~isinf(fe_options.field_output_every_n_frames)
        %Animate result
        figure;
        anim_options.repeat_n_times = 1;
        anim_options.db_range = [-40, 0];
        anim_options.pause_value = 0.001;
        h_patches = fn_show_geometry(main.doms{ref_subdom}.val_mod, main.matls, main.el_types, anim_options);
        fn_run_animation(h_patches, main.doms{ref_subdom}.val.trans{1}.fld, anim_options);
    end
end
