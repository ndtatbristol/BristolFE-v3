function [mod, matls, el_types, steps, fe_options, params] = mod_waveguide(params)
%This is a parametric model description of a 1D waveguide (default =
%circular cross-section)

default_params.els_per_wavelength = 40;
% default_params.abs_bdry_thickness_in_wavelengths = 1;

%Material properties
default_params.matl_name = 'aluminium';

%Details of input signal applied at source
default_params.centre_freq = 0.1e6; %used to determine element size - final result is in frequency domain anyway
default_params.no_cycles = 5;
default_params.source_direction = 2;

%Shape and length of cross section - can be fully specified in 2D, or 2D perimeter
%can be specified or it can just be a rod of specified diameter
default_params.nds_2d = [];
default_params.els_2d = [];
default_params.perimeter_vtcs = [];
default_params.diam = 10e-3;
default_params.z_max = 500e-3;


%Run for long enough for longitudinal waves to travel this many times
%across model (while scattered signal rings down)
default_params.max_time_multiplier = 6; 

%Element shape to use (tri or quad since this is 2d model)
default_params.element_shape = 'tri'; 

%Solver options - specify how ofter field output is produced to use in
%animation
default_params.fe_options.field_output_every_n_frames = inf;

%--------------------------------------------------------------------------
params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
params = fn_set_default_fields(params, default_params);
fe_options = params.fe_options;
el_types = [fn_2d_el_types(), fn_3d_el_types()];

switch params.element_shape
    case {'tri', 'triangular'}
        el_typ_to_use_for_solid_2d = 'CPE3';
        el_typ_to_use_for_fluid_2d = 'AC2D3';
        el_typ_to_use_for_solid_3d = 'C3D6';
        el_typ_to_use_for_fluid_3d = 'AC3D6';
    case {'quad', 'rect', 'quadrilateral', 'rectangle', 'rectangular'}
        el_typ_to_use_for_solid_2d = 'CPE4';
        el_typ_to_use_for_fluid_2d = 'AC2D4';
        el_typ_to_use_for_solid_3d = 'C3D8';
        el_typ_to_use_for_fluid_3d = 'AC3D8';
    otherwise
        error('Unknown element shape')
end

%BUILD THE MODEL USING THE PARAMETERS GIVEN ABOVE

matl_i = 1; 
matls{matl_i} = fn_material_library(params.matl_name);


max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
% abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;

%Define the material
%A cell array with an entry for each material used in the model is required.
%In this case there is only one material, but the index of this material in 
%the cell array is given a name, matl_i, so you can see where it is used
%later when elements in the model are assigned to elements
matl_i = 1; %material index is given a name so you can see where it appears later
matls{matl_i} = fn_material_library(params.matl_name);

%Work out element size and time step
params.el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, params.el_size);

if ~isempty(params.els_2d) && ~isempty(params.nds_2d)
    mod.nds = params.nds_2d;
    mod.els = params.els_2d;
    mod.el_mat_i = zeros(size(mod.els, 1), 1);
    mod.el_typ_i = ones(size(mod.els, 1), 1);
    mod.el_abs_i = zeros(size(mod.els, 1), 1);
else
    if ~isempty(params.perimeter_vtcs)
        bdry_pts = params.perimeter_vtcs;
    else
        
        a = linspace(0, 2 * pi, 361)';
        bdry_pts = params.diam / 2 * [cos(a), sin(a)];
    end
    mod = fn_2d_structured_mesh(bdry_pts, params.el_size, el_typ_to_use_for_solid_2d);
end

%Create the nodes and elements of the mesh
if isfield(params, 'z_pts')
    z_pts = params.z_pts;
else
    z_pts = linspace(0, params.z_max, ceil(params.z_max / params.el_size));
end
mod = fn_extrude_2d_mesh(mod, z_pts);


mod.el_mat_i(:) = matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid_3d)); %extracts the index of the chosen element type from the cell array of possible element types

%Identify all nodes at z = 0
steps{1}.load.frc_nds = find(abs(mod.nds(:, 3)) < params.el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * params.source_direction;

%Provide the time signal for the loading
params.model_size = max(mod.nds(:,3)) - min(mod.nds(:,3));
[vel, ~] = fn_estimate_max_min_vels(matls{matl_i});
max_time = params.model_size / vel * params.max_time_multiplier;
steps{1}.load.time = 0: time_step: max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Say where the displacement should be monitored
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end