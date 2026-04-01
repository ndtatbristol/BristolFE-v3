clear all;
close all;
rng(2)
%ABOUT THIS EXAMPLE
%This example demonstrates multiple features in a 2D model, including
%fluid-solid coupling, absorbing regions, defect generation (crack, 
%inclusion or void). All of these can be turned on or off as required.

%PARAMETRIC DESCRIPTION OF MODEL

%What to include in model
include_fluid_region = 1;
include_absorbing_boundary = 1;
include_crack = 1;
include_scatterer = 1;
scatterer_is_void = 1;

%Geometric description of model
model_size_x = 10e-3;
model_size_y = 12e-3;
interface_pos_y = model_size_y / 2;
interface_angle_degs = 0;
els_per_wavelength = 10;
abs_bdry_thickness_in_wavelengths = 1;
transducer_width = model_size_x / 3;
crack_centre = [model_size_x * 0.4, model_size_y * 0.75];
crack_size = model_size_x * 0.1;

void_centre = [model_size_x * 0.6, model_size_y * 0.65];
void_size = model_size_x * 0.1;


%Solid material properties
solid_longitudinal_velocity = 6300;
solid_shear_velocity = 3150;
solid_density = 2700;
solid_name = 'aluminium';

%Fluid material properties
fluid_velocity = 1500;
fluid_density = 1000;
fluid_name = 'water';

%Details of input signal
centre_freq = 5e6;
no_cycles = 5;

%Run for long enough for longitudinal waves to travel 3 lengths of model
if include_fluid_region
    max_time = 3 * interface_pos_y / fluid_velocity + ...
        3 * (model_size_y - interface_pos_y) / solid_longitudinal_velocity;
else
    max_time = 3 * model_size_y / solid_longitudinal_velocity;
end

fe_options.field_output_every_n_frames = 10;

show_geom_only = 1; %Set to 1 to just show geometry without running model

%--------------------------------------------------------------------------
%THE ACTUAL CODE STARTS HERE
%--------------------------------------------------------------------------

%SET UP THE MODEL

%Define the materials in use
solid_matl_i = 1; 
matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(solid_name, solid_longitudinal_velocity, solid_shear_velocity, solid_density);

fluid_matl_i = 2;
matls{fluid_matl_i} = fn_matl_fluid_defined_by_velocity(fluid_name, fluid_velocity, fluid_density);

%Element types to use
el_typ_to_use_for_solid = 'CPE3'; 
el_typ_to_use_for_fluid = 'AC2D3'; 

%Define shape of model
bdry_pts = [
    0,              0 
    model_size_x,   0 
    model_size_x,   model_size_y 
    0,              model_size_y];

%SET UP THE MODEL

%Add path to BristolFE functions in case not already on path
addpath(genpath('../code'));

%PREPARE THE MESH

%Work out element and time-step size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, el_size);

%Create the nodes and elements of the mesh
% mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size);
el_types = fn_2d_el_types();

%First set material of all elements to steel ...
mod.el_mat_i(:) = solid_matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid));

%... then set elements inside water boundary material to water
if include_fluid_region
    %Define boundary of region that will be water
    water_bdry_pts = [
        0,              0
        model_size_x,   0
        model_size_x,   (1 + tand(interface_angle_degs)) * interface_pos_y
        0,              (1 - tand(interface_angle_degs)) * interface_pos_y];
    %Work out which elements are in that region and assign the material and
    %element type accordingly
    els_in_water = fn_2d_find_elements_in_region(mod, water_bdry_pts);
    mod.el_mat_i(els_in_water) = fluid_matl_i;
    mod.el_typ_i(els_in_water) = find(strcmp(el_types, el_typ_to_use_for_fluid));
    
    %Add interface elements - this is crucial otherwise there will be no
    %coupling between fluid and solid
    mod = fn_add_fluid_solid_interface_els(mod, el_types);
    
    %Set source direction
    transducer_dof = 4; %volumetric source if fluid region is used as source will be in fluid
else
    %Set source direction
    transducer_dof = 2; %forcing in y direction if no fluid is used as source will be on surface of solid
end

if include_crack
    n = 10;
    cod = el_size / 10;
    crack_pts = fn_2d_random_walk(n,crack_size / n, 0, 0, 0, 0.4);
    crack_vtcs = crack_centre + crack_pts - mean(crack_pts);
    mod = fn_2d_add_crack(mod, el_types, crack_vtcs, [], cod);
end

if include_scatterer
    min_rad_frac = 0.5;
    complexity = 3;
    no_pts = 200;
    if scatterer_is_void
        scat_matl_i = 0;
        scat_el_typ_i = 0;
    else
        scat_matl_i = fluid_matl_i;
        scat_el_typ_i = find(strcmp(el_types, el_typ_to_use_for_fluid));
    end
    pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * void_size / 2 + void_centre;
    mod = fn_2d_add_inclusion_or_void(mod, el_types, pts, scat_matl_i, scat_el_typ_i);
end


%Define the absorbing layer
if include_absorbing_boundary
    max_wavelength = fn_estimate_max_min_wavelengths(matls, centre_freq);
    abs_bdry_thickness = abs_bdry_thickness_in_wavelengths * max_wavelength;
    %Define start of absorbing boundary region and its thickness (this may or may not be used)
    abs_bdry_pts = [
        abs_bdry_thickness,                 abs_bdry_thickness
        model_size_x - abs_bdry_thickness,  abs_bdry_thickness
        model_size_x - abs_bdry_thickness,  model_size_y - abs_bdry_thickness
        abs_bdry_thickness,                 model_size_y - abs_bdry_thickness];

    mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);
    [fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq);
end

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
%Define a line along which sources will be placed to excite waves and
%monitoring points will monitor them (in this case they are the same, so it
%is effectively a pulse-echo transducer
transducer_end_pts = [
    model_size_x / 2 + transducer_width / 2, abs_bdry_thickness
    model_size_x / 2 - transducer_width / 2, abs_bdry_thickness];
steps{1}.load.frc_nds = fn_find_nodes_nearest_to_line(mod.nds, transducer_end_pts(1, :), transducer_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * transducer_dof;

steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
if show_geom_only %suppress graphics when running all scripts for testing
    figure; 
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    drawnow
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
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    display_options.draw_elements = 0; %makes it easier to see waves if element edges not drawn
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    anim_options.repeat_n_times = 1;
    anim_options.fld_time = res{1}.fld_time;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end
