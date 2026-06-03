function [mod, matls, el_types, steps, fe_options, params] = mod_2d_advanced(params)

%What to include in model
default_params.include_fluid_region = 1;
default_params.include_absorbing_boundary = 1;
default_params.include_crack = 1;
default_params.include_scatterer = 1;
default_params.scatterer_is_void = 1;

%Geometric description of model
default_params.model_size = [10e-3, 12e-3];
default_params.interface_pos_fraction = 1 / 2;
default_params.interface_angle_degs = 0;
default_params.els_per_wavelength = 10;
default_params.abs_bdry_thickness_in_wavelengths = 1;
default_params.transducer_width_fraction = 1 / 3;

default_params.crack_centre_fraction = [0.4, 0.75];
default_params.crack_size = 1e-3;

default_params.void_centre_fraction = [0.6, 0.65];
default_params.void_size = 1e-3;

%Element types to use
default_params.element_shape = 'tri';

%Solid material properties
default_params.solid_longitudinal_velocity = 6300;
default_params.solid_shear_velocity = 3150;
default_params.solid_density = 2700;
default_params.solid_name = 'aluminium';

%Fluid material properties
default_params.fluid_velocity = 1500;
default_params.fluid_density = 1000;
default_params.fluid_name = 'water';

%Details of input signal
default_params.centre_freq = 5e6;
default_params.no_cycles = 5;
%Run for long enough for longitudinal waves to travel this many lengths of model
default_params.max_time_multiplier = 3; 


default_params.fe_options.field_output_every_n_frames = 10;
default_params.random_seed = 1;

%--------------------------------------------------------------------------
params = fn_set_default_fields(params, default_params);
rng(params.random_seed); 

fe_options = params.fe_options;
el_types = fn_2d_el_types();

%Define the materials in use
solid_matl_i = 1; 
matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(params.solid_name, params.solid_longitudinal_velocity, params.solid_shear_velocity, params.solid_density);

fluid_matl_i = 2;
matls{fluid_matl_i} = fn_matl_fluid_defined_by_velocity(params.fluid_name, params.fluid_velocity, params.fluid_density);

%Define shape of model
bdry_pts = [
    0,              0 
    params.model_size(1),   0 
    params.model_size(1),   params.model_size(2) 
    0,              params.model_size(2)];

%Work out element and time-step size
el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, el_size);

%Create the nodes and elements of the mesh
switch params.element_shape
    case {'tri', 'triangular'}
        el_typ_to_use_for_solid = 'CPE3';
        el_typ_to_use_for_fluid = 'AC2D3';
    case {'quad', 'rect', 'quadrilateral', 'rectangle', 'rectangular'}
        el_typ_to_use_for_solid = 'CPE4';
        el_typ_to_use_for_fluid = 'AC2D4';
    otherwise
        error('Unknown element shape')
end
mod = fn_2d_structured_mesh(bdry_pts, el_size, el_typ_to_use_for_solid);
el_types = fn_2d_el_types();

%First set material of all elements to steel ...
mod.el_mat_i(:) = solid_matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid));

%... then set elements inside water boundary material to water
if params.include_fluid_region
    interface_pos_y = params.model_size(2) * params.interface_pos_fraction;
    %Define boundary of region that will be water
    water_bdry_pts = [
        0,              0
        params.model_size(1),   0
        params.model_size(1),   (1 + tand(params.interface_angle_degs)) * interface_pos_y
        0,              (1 - tand(params.interface_angle_degs)) * interface_pos_y];
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

if params.include_crack
    crack_centre = params.model_size .* params.crack_centre_fraction;
    n = 10;
    cod = el_size / 10;
    crack_pts = fn_2d_random_walk(n, params.crack_size / n, 0, 0, 0, 0.4);
    crack_vtcs = crack_centre + crack_pts - mean(crack_pts);
    mod = fn_2d_add_crack(mod, el_types, crack_vtcs, [], cod);
end

if params.include_scatterer
    void_centre = params.model_size .* params.void_centre_fraction;
    min_rad_frac = 0.5;
    complexity = 3;
    no_pts = 200;
    if params.scatterer_is_void
        scat_matl_i = 0;
        scat_el_typ_i = 0;
    else
        scat_matl_i = fluid_matl_i;
        scat_el_typ_i = find(strcmp(el_types, el_typ_to_use_for_fluid));
    end
    pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * params.void_size / 2 + void_centre;
    mod = fn_2d_add_inclusion_or_void(mod, el_types, pts, scat_matl_i, scat_el_typ_i);
end


%Define the absorbing layer
if params.include_absorbing_boundary
    max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
    abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;
    %Define start of absorbing boundary region and its thickness (this may or may not be used)
    abs_bdry_pts = [
        abs_bdry_thickness,                 abs_bdry_thickness
        params.model_size(1) - abs_bdry_thickness,  abs_bdry_thickness
        params.model_size(1) - abs_bdry_thickness,  params.model_size(2) - abs_bdry_thickness
        abs_bdry_thickness,                 params.model_size(2) - abs_bdry_thickness];

    mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);
    [fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, params.centre_freq);
end

%Run for long enough for longitudinal waves to travel 3 lengths of model
if params.include_fluid_region
    max_time = params.max_time_multiplier * interface_pos_y / params.fluid_velocity + ...
        params.max_time_multiplier * (params.model_size(2) - interface_pos_y) / params.solid_longitudinal_velocity;
else
    max_time = params.max_time_multiplier * params.model_size(2) / params.solid_longitudinal_velocity;
end

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
%Define a line along which sources will be placed to excite waves and
%monitoring points will monitor them (in this case they are the same, so it
%is effectively a pulse-echo transducer
transducer_width = params.transducer_width_fraction * params.model_size(1);
transducer_end_pts = [
    params.model_size(1) / 2 + transducer_width / 2, abs_bdry_thickness
    params.model_size(1) / 2 - transducer_width / 2, abs_bdry_thickness];
steps{1}.load.frc_nds = fn_find_nodes_nearest_to_line(mod.nds, transducer_end_pts(1, :), transducer_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * transducer_dof;

steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end