function [mod, matls, el_types, steps, fe_options, params] = mod_3d_basic(params)

%Solid material properties
default_params.solid_name = 'aluminium';

%Shape of elements to use
default_params.element_shape = 'hex';

%Size of model
default_params.model_size = [10, 11, 12] * 1e-3;

%Spherical void somewhere inside
default_params.include_void = 1;
default_params.void_centre_fraction = [0.6, 0.65, 0.5];
default_params.void_size = 1e-3;

%Source is a circular disk of excitation in centre of top (max z) surface
default_params.src_radius_fraction = 1 / 6;

%Elements per wavelength (higher = more accurate and higher computational cost)
default_params.els_per_wavelength = 6;

default_params.src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
default_params.centre_freq = 2e6;
default_params.no_cycles = 4;
%Run for long enough for longitudinal waves to travel this many lengths of model
default_params.max_time_multiplier = 3; 

default_params.fe_options.field_output_every_n_frames = 10;

%--------------------------------------------------------------------------
%PREPARE THE MESH
if isfield(params, 'fe_options') && isfield(default_params, 'fe_options')
    params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
else
    default_params.fe_options = [];
end
params = fn_set_default_fields(params, default_params);

fe_options = params.fe_options;
el_types = fn_3d_el_types();

switch params.element_shape
    case {'hex', 'hexahedral'}
        el_typ_to_use_for_solid = 'C3D8';
    case {'tet', 'tetrahedral'}
        error('Not implemented yet')
        %el_typ_to_use_for_solid = 'C3D4';
    otherwise
        error('Unknown element shape')
end

%Define the materials in use
solid_matl_i = 1; 
matls{solid_matl_i} = fn_material_library(params.solid_name);

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
mod.el_mat_i = ones(size(mod.el_typ_i)) * solid_matl_i;

if params.include_void
    void_cent = params.model_size .* params.void_centre_fraction;
    [scat_vtcs, scat_fcs] =  fn_3d_spherical_surface(void_cent, params.void_size /2);
    scat_matl_i = 0;
    scat_el_typ_i = [];
    mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i);
end

%Source - a disk of forces on top surface
src_centre = [0.5 * params.model_size(1), 0.5 * params.model_size(2), params.model_size(3)];
src_radius = params.src_radius_fraction * min(params.model_size);

steps{1}.load.frc_nds = fn_find_node_nearest_to_point(mod.nds, src_centre, el_size);
steps{1}.load.frc_nds = find(...
    abs(mod.nds(:, 3) - src_centre(3)) < el_size / 2 & ...
    sqrt(sum( (mod.nds(:, 1:2) - src_centre(1:2)) .^ 2, 2 )) < src_radius ...
    );


steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * params.src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
%Run model for long enough to see first few echoes depending on params.max_time_multiplier
[vel, ~] = fn_estimate_max_min_vels(matls{solid_matl_i});

max_time = 2 * params.max_time_multiplier * params.model_size(3) / vel;
time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end