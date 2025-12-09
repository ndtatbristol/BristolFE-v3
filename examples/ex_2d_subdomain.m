clear all;
close all;
rng(1);
%ABOUT THIS EXAMPLE
%This example is designed to show how to use the subdomain method where a
%main model without scatterers is created and executed and then a subdomain
%model(s) containing a scatterer(s) of interest can be executed. In this example
%script, a validation is also performed by running the full model with the 
%scatterer in included directly. Note that generating field output for
%animations in sub-domain models requires an extra execution of the
%pristine model (because it needs to be executed with an impulse to get the
%transfer functions required for sub-domain modelling, with the output 
%convolved with the desired toneburst, but the animation is
%better if a toneburst input is used directly).

%PARAMETRIC DESCRIPTION OF MODEL
%The model is described in terms of a small number of parameters
%that are defined in the first part of the script. This is good practice as
%it makes it easier to alter.

%Overall model geometry
model_size = 10e-3;
fluid_thickness = 3e-3;
abs_bdry_thickness = 1e-3;
els_per_wavelength = 8;

%Subdomain and scatterer geometry
scatterer_size = 1e-3;
subdomain_size = scatterer_size + 0.1e-3;
scatterer_depth = 3e-3;

%Solid material properties
solid_matl_name = 'steel';
solid_matl_longitudinal_velocity = 5636;
solid_matl_shear_velocity = 3012;
solid_matl_density = 8900;

%Fluid material properties
fluid_name = 'water';
fluid_velocity = 1500;
fluid_density = 1000;

%Element types to use
el_typ_to_use_for_solid = 'CPE3'; 
el_typ_to_use_for_fluid = 'AC2D3'; 

%Details of input signal
centre_freq = 5e6;
no_cycles = 5;

%Transducer details
src_size = 3.5e-3;
src_dir = 4;


%Run the model for long enough to get the first back wall reflection plus a
%bit
max_time = 1.1 * 2 * (fluid_thickness / fluid_velocity + (model_size - fluid_thickness) / solid_matl_longitudinal_velocity);

show_geom_only = 0; %Set to 1 to just show geometry without running model
fe_options.field_output_every_n_frames = inf;20; %set to inf to suppress animations


%--------------------------------------------------------------------------
%THE ACTUAL CODE STARTS HERE
%--------------------------------------------------------------------------

%Add path to BristolFE functions in case not already on path
addpath(genpath('../code'));
addpath(genpath('../subdoms'));
rmpath(genpath('../code/deprecated'));


%Define shape of model
bdry_pts = [
    0, 0
    model_size, 0
    model_size, model_size
    0, model_size];

%Define region that will be fluid
fluid_bdry_pts = [
    0,          0
    model_size, 0
    model_size, fluid_thickness
    0,          fluid_thickness];

%Define start of absorbing boundary region and its thickness
abs_bdry_pts = [
    abs_bdry_thickness,                 abs_bdry_thickness
    model_size - abs_bdry_thickness,    abs_bdry_thickness
    model_size - abs_bdry_thickness,    model_size
    abs_bdry_thickness,                 model_size];

%Define the materials
%A cell array with an entry for each material used in the model is required.
solid_matl_i = 1;
main.matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(solid_matl_name, solid_matl_longitudinal_velocity, solid_matl_shear_velocity, solid_matl_density);
fluid_matl_i = 2;
main.matls{fluid_matl_i} = fn_matl_fluid_defined_by_velocity(fluid_name, fluid_velocity, fluid_density);

%Other stuff
fe_options.dof_to_use = [1,2,4];%x, y and pressure

%Work out element size and Create the nodes and elements of the mesh
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);
main.mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
main.el_types = fn_2d_el_types();

%First set material of all elements to solid ...
main.mod.el_mat_i(:) = solid_matl_i;
main.mod.el_typ_i(:) = find(strcmp(main.el_types, el_typ_to_use_for_solid));

%... then set elements inside fluid boundary material to fluid
els_in_fluid = fn_2d_find_elements_in_region(main.mod, fluid_bdry_pts);
main.mod.el_mat_i(els_in_fluid) = fluid_matl_i;
main.mod.el_typ_i(els_in_fluid) = find(strcmp(main.el_types, el_typ_to_use_for_fluid));

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
main.mod = fn_add_fluid_solid_interface_els(main.mod, main.el_types);

%Time step and max time
main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size);
main.mod.design_centre_freq = centre_freq;
fe_options.time_pts = ceil(max_time / main.mod.max_safe_time_step);

%Define the absorbing layer
main.mod = fn_2d_add_absorbing_layer(main.mod, abs_bdry_pts, abs_bdry_thickness);

%Define transducer
src_end_pts = [ model_size / 2 - src_size / 2, abs_bdry_thickness
                model_size / 2 + src_size / 2, abs_bdry_thickness];
[main.trans{1}.nds, s] = fn_find_nodes_nearest_to_line(main.mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
main.trans{1}.dfs = ones(size(main.trans{1}.nds)) * src_dir;

%Create a subdomain in the middle with a hole in surface as scatterer
scatterer_centre = [model_size / 2, fluid_thickness + scatterer_depth];
inner_bdry = [-1,-1;-1,1;1,1;1,-1] / 2 * subdomain_size + scatterer_centre;
scat_pts =   fn_2d_create_smooth_random_blob(0.4, 3, 360) * scatterer_size / 2 + scatterer_centre;

main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, abs_bdry_thickness);
main.doms{1}.mod = fn_2d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_pts, 0, 0);

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

%Run sub-domain model
main = fn_run_subdomain_model(main, fe_options);

%Animate results if requested
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    if ~isinf(fe_options.field_output_every_n_frames)
        figure;
        anim_options.repeat_n_times = 1;
        anim_options.db_range = [-40, 0];
        anim_options.pause_value = 0.001;
        h_patches = fn_show_geometry_with_subdomains(main, anim_options);
        fn_run_subdomain_animations(main, h_patches, anim_options);
    end
end

%Run validation model
fe_options.validation_mode = 1;
main = fn_run_main_model(main, fe_options);

%Animate validation results if requested
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    %View the time domain data and compare wih validation
    figure;
    i = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
    mv = max(abs(sum(main.doms{1}.res.fmc.time_data(i:end,: ), 2)));
    plot(main.doms{1}.res.fmc.time, real(sum(main.doms{1}.res.fmc.time_data, 2)) / mv, 'k', 'LineWidth', 2);
    hold on;
    plot(main.doms{1}.val.fmc.time, real(sum(main.doms{1}.val.fmc.time_data, 2)) / mv, 'g:', 'LineWidth', 2);
    plot(main.res.fmc.time, real(sum(main.res.fmc.time_data, 2)) / mv, 'b');
    ylim([-1,1]);
    yyaxis right
    plot(main.doms{1}.res.fmc.time, 20 * log10(abs(sum(main.doms{1}.res.fmc.time_data, 2) - sum(main.doms{1}.val.fmc.time_data, 2)) / mv));
    ylim([-60, 0]);
    legend('Sub-domain method', 'Validation', 'Pristine', 'Difference (dB)');

    if ~isinf(fe_options.field_output_every_n_frames)
        %Animate result
        figure;
        anim_options.repeat_n_times = 1;
        anim_options.db_range = [-40, 0];
        anim_options.pause_value = 0.001;
        h_patches = fn_show_geometry(main.doms{1}.val_mod, main.matls, main.el_types, anim_options);
        fn_run_animation(h_patches, main.doms{1}.val.trans{1}.fld, anim_options);
    end
end
