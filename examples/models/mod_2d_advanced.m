function [mod, matls, el_types, steps, fe_options, params] = mod_2d_advanced(params)

%By default, this model uses predictor-corrector (pc) solver mode which is
%the fastest, but not guaranteed to be stable.

%What to include in model
default_params.include_fluid_region = 1;
default_params.include_absorbing_boundary = 1;
default_params.include_crack = 1;
default_params.include_scatterer = 1;
default_params.scatterer_is_void = 1;

%Geometric description of model
default_params.model_size = [8e-3, 12e-3];
default_params.interface_pos_fraction = 1 / 3;
default_params.interface_angle_degs = 0;
default_params.els_per_wavelength = 10;
default_params.abs_bdry_thickness_in_wavelengths = 1;
default_params.transducer_width_fraction_of_x = 1 / 3;

default_params.crack_centre_fraction = [0.4, 0.6];
default_params.crack_size = 1e-3;

default_params.scatterer_centre_fraction = [0.6, 0.75];
default_params.scatterer_size = 1e-3;

%Element types to use
default_params.element_shape = 'tri';

%Solid material properties of part where scatterers are
default_params.solid1_name = 'aluminium';

%Fluid material properties
default_params.fluid_name = 'water';

%Solid material properties of part where sources are (if fluid not used)
default_params.solid2_name = 'perspex';

%Details of input signal
default_params.centre_freq = 5e6;
default_params.no_cycles = 5;
%Run for long enough for longitudinal waves to travel this many lengths of model
default_params.max_time_multiplier = 3; 
default_params.safety_factor = 1.5;

default_params.random_seed = 1;
default_params.fe_options.field_output_every_n_frames = 20;
default_params.fe_options.solver_mode = 'pc';

%--------------------------------------------------------------------------
params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
params = fn_set_default_fields(params, default_params);
rng(params.random_seed); 

fe_options = params.fe_options;
el_types = fn_2d_el_types();

%Define the materials in use
matl1_i = 1; 
matls{matl1_i} = fn_material_library(params.solid1_name);

matl2_i = 2;
if params.include_fluid_region
    matls{matl2_i} = fn_material_library(params.fluid_name);
else
    matls{matl2_i} = fn_material_library(params.solid2_name);
end

%Define shape of model
bdry_pts = [
    0,              0 
    params.model_size(1),   0 
    params.model_size(1),   params.model_size(2) 
    0,              params.model_size(2)];

%Work out element and time-step size
params.el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, params.el_size, params.safety_factor);

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
mod = fn_2d_structured_mesh(bdry_pts, params.el_size, el_typ_to_use_for_solid);
el_types = fn_2d_el_types();

%First set material of all elements to solid material 1
mod.el_mat_i(:) = matl1_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid));

%Then identify all elements in second material (either solid2 or fluid)
interface_pos_y = params.model_size(2) * params.interface_pos_fraction;
%Define boundary of region that will be water
matl2_bdry_pts = [
    0,              0
    params.model_size(1),   0
    params.model_size(1),   (1 + tand(params.interface_angle_degs)) * interface_pos_y
    0,              (1 - tand(params.interface_angle_degs)) * interface_pos_y];
%Work out which elements are in that region and assign the material and
%element type accordingly
els_in_matl2 = fn_2d_find_elements_in_region(mod, matl2_bdry_pts);
mod.el_mat_i(els_in_matl2) = matl2_i;
if params.include_fluid_region
    %Set element type to fluid
    mod.el_typ_i(els_in_matl2) = find(strcmp(el_types, el_typ_to_use_for_fluid));
    
    %Add interface elements - this is crucial otherwise there will be no
    %coupling between fluid and solid
    mod = fn_add_fluid_solid_interface_els(mod, el_types);
    
    %Set source direction for fluid
    transducer_dof = 4; %volumetric source if fluid region is used as source will be in fluid
else
    %Set source direction for solid
    transducer_dof = 2; %forcing in y direction if no fluid is used as source will be on surface of solid
end

crack_centre = params.model_size .* params.crack_centre_fraction;
n = 10;
cod = params.el_size / 10;
crack_pts = fn_2d_random_walk(n, params.crack_size / n, 0, 0, 0, 0.1);
params.crack_vtcs = crack_centre + crack_pts - mean(crack_pts);
if params.include_crack
    mod = fn_2d_add_crack(mod, el_types, params.crack_vtcs, [], cod);
end

void_centre = params.model_size .* params.scatterer_centre_fraction;
min_rad_frac = 0.5;
complexity = 3;
no_pts = 200;
if params.scatterer_is_void
    scat_matl_i = 0;
    scat_el_typ_i = 0;
else
    scat_matl_i = matl2_i;
    scat_el_typ_i = find(strcmp(el_types, el_typ_to_use_for_fluid));
end
params.void_pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * params.scatterer_size / 2 + void_centre;
if params.include_scatterer
    mod = fn_2d_add_inclusion_or_void(mod, el_types, params.void_pts, scat_matl_i, scat_el_typ_i);
end


%Define the absorbing layer
max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
params.abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;
if params.include_absorbing_boundary
    %Define start of absorbing boundary region and its thickness (this may or may not be used)
    abs_bdry_pts = [
        params.abs_bdry_thickness,                 params.abs_bdry_thickness
        params.model_size(1) - params.abs_bdry_thickness,  params.abs_bdry_thickness
        params.model_size(1) - params.abs_bdry_thickness,  params.model_size(2) - 0
        params.abs_bdry_thickness,                 params.model_size(2) - 0];

    mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, params.abs_bdry_thickness);
    [fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(params.abs_bdry_thickness, matls, params.centre_freq);
else
    params.abs_bdry_thickness = 0;
end

%Run for long enough for longitudinal waves to travel 3 lengths of model
[vel1, ~] = fn_estimate_max_min_vels(matls{matl1_i});
[vel2, ~] = fn_estimate_max_min_vels(matls{matl2_i});
max_time = params.max_time_multiplier * interface_pos_y / vel2 + ...
    params.max_time_multiplier * (params.model_size(2) - interface_pos_y) / vel1;

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
%Define a line along which sources will be placed to excite waves and
%monitoring points will monitor them (in this case they are the same, so it
%is effectively a pulse-echo transducer
transducer_width = params.transducer_width_fraction_of_x * params.model_size(1);
transducer_end_pts = [
    params.model_size(1) / 2 + transducer_width / 2, params.abs_bdry_thickness
    params.model_size(1) / 2 - transducer_width / 2, params.abs_bdry_thickness];
steps{1}.load.frc_nds = fn_find_nodes_nearest_to_line(mod.nds, transducer_end_pts(1, :), transducer_end_pts(2, :), params.el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * transducer_dof;

steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end