clear all;
close all;

rng(1);

%ABOUT THIS EXAMPLE
%This example is designed to show how to use the subdomain method for a
%classic NDE case: oblique incidence shear waves using a wedge transducer
%to detect surface-breaking cracks.

%PARAMETRIC DESCRIPTION OF MODEL
%The model is described in terms of a small number of parameters
%that are defined in the first part of the script. This is good practice as
%it makes it easier to alter.

%Overall model geometry
plate_thickness = 10e-3;
clearance = 1e-3;
abs_bdry_thickness_in_wavelengths = 1;
els_per_wavelength = 5;
inc_angle_degs = 45; %wedge geometry will be calculated from this and inc_angle mode
inc_mode = 'S';

%Subdomain geometry
subdomain_size_x = 5e-3;
subdomain_size_y = 5e-3;

%Solid material properties
plate_matl_name = 'steel';
plate_matl_longitudinal_velocity = 5636;
plate_matl_shear_velocity = 3012;
plate_matl_density = 8900;

%Wedge material properties
wedge_matl_name = 'perspex';
wedge_matl_longitudinal_velocity = 2730;
wedge_matl_shear_velocity = 1345;
wedge_matl_density = 1190;

%Element types to use
el_typ_to_use_for_solid = 'CPE3';

%Details of input signal
centre_freq = 5e6;
no_cycles = 5;

%Transducer details
trans_diam = 5e-3;

show_geom_only = 0; %Set to 1 to just show geometry without running model
run_validation_models = 1;
fe_options.field_output_every_n_frames = inf;10; %set to inf to suppress animations

%--------------------------------------------------------------------------
%END OF INPUTS

%Define the materials
%A cell array with an entry for each material used in the model is required.
plate_matl_i = 1;
main.matls{plate_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(plate_matl_name, plate_matl_longitudinal_velocity, plate_matl_shear_velocity, plate_matl_density);
wedge_matl_i = 2;
main.matls{wedge_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(wedge_matl_name, wedge_matl_longitudinal_velocity, wedge_matl_shear_velocity, wedge_matl_density, [255, 255, 224] / 256);

%Define thickness and start of absorbing boundary region
max_wavelength = fn_estimate_max_min_wavelengths(main.matls, centre_freq);
abs_bdry_thickness = abs_bdry_thickness_in_wavelengths * max_wavelength;

switch lower(inc_mode)
    case {'shear', 's', 't', 'transverse'}
        plate_vel = plate_matl_shear_velocity;
    case {'longitudinal',  'l'}
        plate_vel = plate_matl_longitudinal_velocity;
end
wedge_angle_degs = asind(wedge_matl_longitudinal_velocity / plate_vel * sind(inc_angle_degs));
force_angle_degs = wedge_angle_degs + 90; %force angles are w.r.t. x axis (i.e. usual convention)

%Geometry calculations. Defect is nominally at (0,0)
trans_rad = trans_diam / 2;
wedge_path_length = trans_rad * tand(wedge_angle_degs) + clearance;
wedge_axial_len = wedge_path_length + clearance + abs_bdry_thickness;
wedge_half_width = abs_bdry_thickness + trans_rad;
plate_path_length = plate_thickness / cosd(inc_angle_degs);

xmax = subdomain_size_x / 2 + clearance + abs_bdry_thickness;
xwedge_centre = -plate_thickness * tand(inc_angle_degs);
xmin_plate = xwedge_centre - trans_rad / cosd(inc_angle_degs) - clearance - abs_bdry_thickness;
xmin_wedge = xwedge_centre - wedge_axial_len * sind(wedge_angle_degs) - wedge_half_width * cosd(wedge_angle_degs);
xmin = min(xmin_wedge, xmin_plate);
ymin = 0;
ymax = wedge_axial_len * cosd(wedge_angle_degs) + wedge_half_width * sind(wedge_angle_degs) + plate_thickness;

%Maximum time
max_time = 1.5 * 2 * (wedge_path_length / wedge_matl_longitudinal_velocity + plate_path_length / plate_vel);

%Add path to BristolFE functions in case not already on path
addpath(genpath('../code'));
addpath(genpath('../subdoms'));
rmpath(genpath('../code/deprecated'));


%Define shape of model
bdry_pts = [
    xmin, 0
    xmax, 0
    xmax, ymax
    xmin, ymax];

%absorbing boundary for plate
abs_bdry_pts_plate = [
    xmin_plate + abs_bdry_thickness, 0
    xmin_plate + abs_bdry_thickness, plate_thickness
    xmax - abs_bdry_thickness, plate_thickness
    xmax - abs_bdry_thickness, 0];

%absorbing boundary for wedge - note these extend into plate material to
%make geometry easier, but absorption is only added to wedge material
tmp_len = wedge_path_length + clearance;
abs_bdry_pts_wedge = [
    xwedge_centre + tmp_len * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), plate_thickness - tmp_len * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)
    xwedge_centre + tmp_len * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), plate_thickness - tmp_len * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - tmp_len * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), plate_thickness + tmp_len * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - tmp_len * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), plate_thickness + tmp_len * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)];

%Other stuff
fe_options.dof_to_use = [1,2,3];

%Work out element size, create the nodes and elements of the mesh,
%determine elements in plate and wedge based on y value
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);
main.mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
main.el_types = fn_2d_el_types();
el_ctrs = fn_calc_element_centres(main.mod.nds, main.mod.els);
plate_els = el_ctrs(:, 2) <= plate_thickness;
wedge_els = el_ctrs(:, 2) > plate_thickness;
plate_els_to_go = plate_els & el_ctrs(:, 1) < xmin_plate;
main.mod.el_typ_i(:) = find(strcmp(main.el_types, el_typ_to_use_for_solid));

%Set materials
main.mod.el_mat_i(plate_els) = plate_matl_i;
main.mod.el_mat_i(wedge_els) = wedge_matl_i;

%Time step and max time
% main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size);
% main.mod.design_centre_freq = centre_freq;
% fe_options.time_pts = ceil(max_time / main.mod.max_safe_time_step);

%Define the absorbing layer
main.mod = fn_2d_add_absorbing_layer(main.mod, abs_bdry_pts_plate, abs_bdry_thickness, plate_els);
main.mod = fn_2d_add_absorbing_layer(main.mod, abs_bdry_pts_wedge, abs_bdry_thickness, wedge_els);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, main.matls, centre_freq);

wedge_els_to_go = wedge_els & (main.mod.el_abs_i == 1);

%Lose the elements that need to go
els_to_go = plate_els_to_go | wedge_els_to_go;
main.mod.els(els_to_go, :) = [];
main.mod.el_mat_i(els_to_go) = [];
main.mod.el_abs_i(els_to_go) = [];
main.mod.el_typ_i(els_to_go) = [];
[main.mod.nds, main.mod.els, ~, ~] = fn_remove_unused_nodes(main.mod.nds, main.mod.els);

%Define transducer
trans_end_pts = [
    xwedge_centre - wedge_path_length * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), plate_thickness + wedge_path_length * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - wedge_path_length * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), plate_thickness + wedge_path_length * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)];
nds = fn_find_nodes_nearest_to_line(main.mod.nds, trans_end_pts(1, :), trans_end_pts(2, :), el_size / 2);
main.trans{1}.nds = [
    nds
    nds];
main.trans{1}.dfs = [
    ones(size(nds)) * 1
    ones(size(nds)) * 2
    ];
main.trans{1}.wts = [
    ones(size(nds)) * cosd(force_angle_degs)
    ones(size(nds)) * sind(force_angle_degs)
    ];

%Input signal
time_step = fn_get_suitable_time_step(main.matls, el_size);
time_pts = ceil(max_time / time_step);
main.inp.time = [0:time_pts - 1] * time_step;
main.inp.sig = fn_gaussian_pulse(main.inp.time, centre_freq, no_cycles);

%Create a subdomain in the middle with a hole in surface as scatterer
inner_bdry = [
    -subdomain_size_x / 2, 0
    -subdomain_size_x / 2, subdomain_size_y
    subdomain_size_x / 2, subdomain_size_y
    subdomain_size_x / 2, 0];

empty_subdomain = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, abs_bdry_thickness);

%Sub-domain 1 contains a large notch like a corner echo as a reference
main.doms{1}.mod = empty_subdomain;
scat_pts = [
    0, 0
    subdomain_size_x / 2 - el_size, 0
    subdomain_size_x / 2 - el_size, subdomain_size_y - el_size
    0, subdomain_size_y - el_size];
main.doms{1}.mod = fn_2d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_pts, 0, 0);

%Sub-domain 2 contains a random crack
main.doms{2}.mod = empty_subdomain;
crack_len = subdomain_size_y;
crack_vtcs = fn_2d_random_walk(50, crack_len / 50, 0, pi / 2, pi / 8);
main.doms{2}.mod = fn_2d_add_crack(main.doms{2}.mod, main.el_types, crack_vtcs);



%Show the mesh
if show_geom_only %suppress graphics when running all scripts for testing
    figure;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry_with_subdomains(main, display_options);
    drawnow
    return
end
%--------------------------------------------------------------------------

%Run main model
main = fn_run_main_model(main, fe_options);

%Subdomain model with no scatterers
main = fn_run_subdomain_model(main, fe_options);


% %Demonstration of how sub-domain model can be run for multiple random scatterers
% results = zeros(numel(main.inp.time), no_scatterers);
% for s = 1:no_scatterers
%     scat_pts =   fn_2d_create_smooth_random_blob(0.4, 3, 360) * scatterer_size / 2 + scatterer_centre;
%     main.doms{1}.mod = fn_2d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_pts, 0, 0);
%     main = fn_run_subdomain_model(main, fe_options);
%     results(:,s) = sum(main.doms{1}.res.fmc.time_data, 2);
% end

% figure;
% plot(main.inp.time, real(results));

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
    i = min((find(main.res.fmc.time >  2 * (wedge_path_length / wedge_matl_longitudinal_velocity + plate_path_length / plate_vel))));
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
        ylim([-60, 0]);
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
