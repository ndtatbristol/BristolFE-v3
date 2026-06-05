function [mod, matls, el_types, steps, fe_options, params] = mod_3d_advanced(params)

%What to include in model
default_params.include_fluid_region = 1;
default_params.include_absorbing_boundary = 1;
default_params.include_crack = 1;
default_params.include_scatterer = 1;
default_params.scatterer_is_void = 1;

%Solid material properties of part where scatterers are
default_params.solid1_name = 'aluminium';

%Fluid material properties
default_params.fluid_name = 'water';

%Solid material properties of part where sources are (if fluid not used)
default_params.solid2_name = 'perspex';

%Shape of elements to use
default_params.element_shape = 'hex';

%Size of model
default_params.model_size = [8, 4, 12] * 1e-3;
default_params.interface_pos_fraction = 1 / 2;
default_params.interface_angle_degs = 0; %will tilt interface in x-z plane
default_params.abs_bdry_thickness_in_wavelengths = 1;

%Scatterer somewhere in top half
default_params.scatterer_centre_fraction = [0.6, 0.5, 0.75];
default_params.scatterer_size = 1e-3;

%Crack somewhere in top half
default_params.crack_centre_fraction = [0.4, 0.5, 0.75];
default_params.crack_size = 1e-3;
default_params.crack_angle_degs = 30; %will tilt crack in x-z plane


%Source is a circular disk of excitation in centre of top (max z) surface
default_params.transducer_width_fraction_of_x = 1 / 3;

%Elements per wavelength (higher = more accurate and higher computational cost)
default_params.els_per_wavelength = 6;

% default_params.src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
default_params.centre_freq = 5e6;
default_params.no_cycles = 5;
%Run for long enough for longitudinal waves to travel this many lengths of model
default_params.max_time_multiplier = 3; 

default_params.fe_options.field_output_every_n_frames = 10;

%--------------------------------------------------------------------------
%PREPARE THE MESH

params = fn_set_default_fields(params, default_params);

fe_options = params.fe_options;
el_types = fn_3d_el_types();

switch params.element_shape
    case {'hex', 'hexahedral'}
        el_typ_to_use_for_solid = 'C3D8';
        el_typ_to_use_for_fluid = 'AC3D8';
    case {'tet', 'tetrahedral'}
        error('Not implemented yet')
        %el_typ_to_use_for_solid = 'C3D4';
        %el_typ_to_use_for_fluid = 'AC3D4';
    otherwise
        error('Unknown element shape')
end

%Define the materials in use
matl1_i = 1; 
matls{matl1_i} = fn_material_library(params.solid1_name);

matl2_i = 2;
if params.include_fluid_region
    matls{matl2_i} = fn_material_library(params.fluid_name);
else
    matls{matl2_i} = fn_material_library(params.solid2_name);
end

%Work out element size
el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);

%Create the nodes and elements of the mesh - need similar general function
%as for 2D meshes at some point
crnr_pts = [
    0, 0, 0
    params.model_size];
switch params.element_shape
    case {'hex', 'hexahedral'}
        mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
end

el_types = fn_3d_el_types(); 

solid_el_typ_i = find(strcmp(el_types, el_typ_to_use_for_solid));

mod.el_typ_i = ones(size(mod.el_typ_i)) * solid_el_typ_i;
mod.el_mat_i = ones(size(mod.el_typ_i)) * matl1_i;

%Then identify all elements in second material (either solid2 or fluid)
interface_pos_z = params.model_size(3) * params.interface_pos_fraction;

%Define boundary of region that will be second material
matl2_bdry_pts = [
    0,              0, 0 
    params.model_size(1),   0, 0 
    params.model_size(1),   params.model_size(2), 0 
    0,   params.model_size(2), 0 
    0,              0, (1 - tand(params.interface_angle_degs)) * interface_pos_z; 
    params.model_size(1),   0, (1 + tand(params.interface_angle_degs)) * interface_pos_z 
    params.model_size(1),   params.model_size(2), (1 + tand(params.interface_angle_degs)) * interface_pos_z 
    0,   params.model_size(2), (1 - tand(params.interface_angle_degs)) * interface_pos_z];
[vtmatl2_bdry_ptscs, matl2_bdry_fcs] = fn_3d_hexahedral_surface(matl2_bdry_pts);

%Work out which elements are in that region and assign the material and
%element type accordingly
els_in_matl2 = fn_3d_find_elements_in_region(mod, matl2_bdry_pts, matl2_bdry_fcs);
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
    transducer_dof = 3; %forcing in z direction if no fluid is used as source will be on surface of solid
end

if params.include_crack
    crack_centre = params.model_size .* params.crack_centre_fraction;
    [crack_vtcs, crack_fcs] = fn_3d_elliptical_surface(crack_centre, [-sind(params.crack_angle_degs), 0, cosd(params.crack_angle_degs)], [0,1,0], params.crack_size / 2, params.crack_size / 2);
    mod = fn_3d_add_crack(mod, el_types, crack_vtcs, crack_fcs, 0);
end

if params.include_scatterer
    void_cent = params.model_size .* params.scatterer_centre_fraction;
    [scat_vtcs, scat_fcs] =  fn_3d_spherical_surface(void_cent, params.scatterer_size /2);
    scat_matl_i = 0;
    scat_el_typ_i = [];
    mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i);
end

%Define the absorbing layer
if params.include_absorbing_boundary
    max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
    abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;
    pt1 = [0, 0, 0] + abs_bdry_thickness;
    pt2 = params.model_size - abs_bdry_thickness;
    [abs_bdry_pts, abs_bdry_fcs] = fn_3d_rectalinear_surface(pt1, pt2);
    mod = fn_3d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_fcs, abs_bdry_thickness);
    [fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, params.centre_freq);
else
    abs_bdry_thickness = 0;
end

%Source - a disk of forces / pressures in lower material
src_centre = [0.5 * params.model_size(1), 0.5 * params.model_size(2), abs_bdry_thickness];
src_normal = [0,0,1]; 
src_radius = params.transducer_width_fraction_of_x * params.model_size(1) / 2;
[bdry_nds, bdry_fcs] = fn_3d_disk_surface(src_centre, src_centre - src_normal, src_radius);
d = fn_3d_signed_dist_to_bdry(mod.nds, bdry_nds, bdry_fcs);
tmp = find(abs(d) < el_size / 2);
steps{1}.load.frc_nds = tmp(:);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * transducer_dof;

%Run for long enough for longitudinal waves to travel 3 lengths of model
[vel1, ~] = fn_estimate_max_min_vels(matls{matl1_i});
[vel2, ~] = fn_estimate_max_min_vels(matls{matl2_i});
max_time = params.max_time_multiplier * interface_pos_z / vel2 + ...
    params.max_time_multiplier * (params.model_size(3) - interface_pos_z) / vel1;

time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end