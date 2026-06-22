function [mod, matls, el_types, steps, fe_options, params] = mod_2d_oblique(params)

%This example is designed to show a common NDE set up.

%What to include in the geometry
default_params.include_crack = 1;
default_params.reference_notch = 0; %this takes precedence over crack if set

%Overall model geometry
default_params.plate_thickness = 10e-3;
default_params.clearance = 1e-3;
default_params.abs_bdry_thickness_in_wavelengths = 1;
default_params.els_per_wavelength = 10;
default_params.inc_angle_degs = 45; %wedge geometry will be calculated from this and inc_angle mode
default_params.inc_mode = 'S';

%Bit where crack will go
default_params.defect_region_size = [5e-3, 5e-3];

%Solid material properties
default_params.plate_matl_name = 'steel';

%Wedge material properties
default_params.wedge_matl_name = 'perspex';

%Element types to use
default_params.element_shape = 'tri';

%Details of input signal
default_params.centre_freq = 5e6;
default_params.no_cycles = 5;

%Transducer details
default_params.trans_diam = 5e-3;

default_params.random_seed = 1;

default_params.fe_options.field_output_every_n_frames = 20; %set to inf to suppress animations

%--------------------------------------------------------------------------
if isfield(params, 'fe_options') && isfield(default_params, 'fe_options')
    params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
else
    default_params.fe_options = [];
end
params = fn_set_default_fields(params, default_params);
fe_options = params.fe_options;
rng(params.random_seed); 

%Define the materials
%A cell array with an entry for each material used in the model is required.
plate_matl_i = 1;
matls{plate_matl_i} = fn_material_library(params.plate_matl_name);
[plate_matl_longitudinal_velocity, plate_matl_shear_velocity] = fn_estimate_max_min_vels(matls{plate_matl_i});

wedge_matl_i = 2;
matls{wedge_matl_i} = fn_material_library(params.wedge_matl_name);
[wedge_matl_longitudinal_velocity, wedge_matl_shear_velocity] = fn_estimate_max_min_vels(matls{wedge_matl_i});

%Define thickness and start of absorbing boundary region
max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
params.abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;

switch lower(params.inc_mode)
    case {'shear', 's', 't', 'transverse'}
        plate_vel = plate_matl_shear_velocity;
    case {'longitudinal',  'l'}
        plate_vel = plate_matl_longitudinal_velocity;
end
wedge_angle_degs = asind(wedge_matl_longitudinal_velocity / plate_vel * sind(params.inc_angle_degs));
force_angle_degs = wedge_angle_degs + 90; %force angles are w.r.t. x axis (i.e. usual convention)

%Geometry calculations. Defect is nominally at (0,0)
trans_rad = params.trans_diam / 2;
wedge_path_length = trans_rad * tand(wedge_angle_degs) + params.clearance;
wedge_axial_len = wedge_path_length + params.clearance + params.abs_bdry_thickness;
wedge_half_width = params.abs_bdry_thickness + trans_rad;
plate_path_length = params.plate_thickness / cosd(params.inc_angle_degs);

xmax = params.defect_region_size(1) / 2 + params.clearance + params.abs_bdry_thickness;
xwedge_centre = -params.plate_thickness * tand(params.inc_angle_degs);
xmin_plate = xwedge_centre - trans_rad / cosd(params.inc_angle_degs) - params.clearance - params.abs_bdry_thickness;
xmin_wedge = xwedge_centre - wedge_axial_len * sind(wedge_angle_degs) - wedge_half_width * cosd(wedge_angle_degs);
xmin = min(xmin_wedge, xmin_plate);
ymin = 0;
ymax = wedge_axial_len * cosd(wedge_angle_degs) + wedge_half_width * sind(wedge_angle_degs) + params.plate_thickness;

%Maximum time
max_time = 1.5 * 2 * (wedge_path_length / wedge_matl_longitudinal_velocity + plate_path_length / plate_vel);

%Define shape of model
bdry_pts = [
    xmin, 0
    xmax, 0
    xmax, ymax
    xmin, ymax];

%absorbing boundary for plate
abs_bdry_pts_plate = [
    xmin_plate + params.abs_bdry_thickness, 0
    xmin_plate + params.abs_bdry_thickness, params.plate_thickness
    xmax - params.abs_bdry_thickness, params.plate_thickness
    xmax - params.abs_bdry_thickness, 0];

%absorbing boundary for wedge - note these extend into plate material to
%make geometry easier, but absorption is only added to wedge material
tmp_len = wedge_path_length + params.clearance;
abs_bdry_pts_wedge = [
    xwedge_centre + tmp_len * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), params.plate_thickness - tmp_len * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)
    xwedge_centre + tmp_len * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), params.plate_thickness - tmp_len * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - tmp_len * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), params.plate_thickness + tmp_len * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - tmp_len * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), params.plate_thickness + tmp_len * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)];

%Other stuff
fe_options.dof_to_use = [1,2,3];

%Work out element size, create the nodes and elements of the mesh,
%determine elements in plate and wedge based on y value
params.el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
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

el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
plate_els = el_ctrs(:, 2) <= params.plate_thickness;
wedge_els = el_ctrs(:, 2) > params.plate_thickness;
plate_els_to_go = plate_els & el_ctrs(:, 1) < xmin_plate;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid));

%Set materials
mod.el_mat_i(plate_els) = plate_matl_i;
mod.el_mat_i(wedge_els) = wedge_matl_i;

%Define the absorbing layer
mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts_plate, params.abs_bdry_thickness, plate_els);
mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts_wedge, params.abs_bdry_thickness, wedge_els);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(params.abs_bdry_thickness, matls, params.centre_freq);

wedge_els_to_go = wedge_els & (mod.el_abs_i == 1);

%Lose the elements that need to go
els_to_go = plate_els_to_go | wedge_els_to_go;
mod.els(els_to_go, :) = [];
mod.el_mat_i(els_to_go) = [];
mod.el_abs_i(els_to_go) = [];
mod.el_typ_i(els_to_go) = [];
[mod.nds, mod.els, ~, ~] = fn_remove_unused_nodes(mod.nds, mod.els);

%Define crack and notch geometry (may or may not be included in model
%defined, but points are used anyway in subdomain examples
crack_len = params.defect_region_size(2);
crack_pts = 50;
params.crack_vtcs = fn_2d_random_walk(crack_pts, crack_len / crack_pts, 0, pi / 2, 0, pi / 4 / sqrt(crack_pts));
params.notch_pts = [
    0, 0
    params.defect_region_size(1) / 2, 0
    params.defect_region_size(1) / 2, params.defect_region_size(2)
    0, params.defect_region_size(2)];
if params.reference_notch
    mod = fn_2d_add_inclusion_or_void(mod, el_types, params.notch_pts, 0, 0);
else
    if params.include_crack
        mod = fn_2d_add_crack(mod, el_types, params.crack_vtcs);
    end
end

%Define transducer
trans_end_pts = [
    xwedge_centre - wedge_path_length * sind(wedge_angle_degs) + trans_rad * cosd(wedge_angle_degs), params.plate_thickness + wedge_path_length * cosd(wedge_angle_degs) + trans_rad * sind(wedge_angle_degs)
    xwedge_centre - wedge_path_length * sind(wedge_angle_degs) - trans_rad * cosd(wedge_angle_degs), params.plate_thickness + wedge_path_length * cosd(wedge_angle_degs) - trans_rad * sind(wedge_angle_degs)];
nds = fn_find_nodes_nearest_to_line(mod.nds, trans_end_pts(1, :), trans_end_pts(2, :), params.el_size / 2);
steps{1}.load.frc_nds = [
    nds
    nds];
steps{1}.load.frc_dfs = [
    ones(size(nds)) * 1
    ones(size(nds)) * 2
    ];
steps{1}.load.wts = [
    ones(size(nds)) * cosd(force_angle_degs)
    ones(size(nds)) * sind(force_angle_degs)
    ];

steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Input signal
time_step = fn_get_suitable_time_step(matls, params.el_size);
time_pts = ceil(max_time / time_step);
steps{1}.load.time = [0:time_pts - 1] * time_step;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);
end
