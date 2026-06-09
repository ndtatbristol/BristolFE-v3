function [mod, matls, el_types, steps, fe_options, params] = mod_3d_fastener_hole(params)
default_params.els_per_wavelength = 4;
default_params.corner_calibration_sample = 0;
default_params.trans_x_offset = 0;
default_params.max_time_safety_factor = 2;
default_params.absorption_in_component = 1;
default_params.absorption_in_wedge = 1;

params = fn_set_default_fields(params, default_params);

%PARAMETRIC DESCRIPTION OF MODEL
thickness = 8e-3;
trans_diam = 5e-3;
hole_diam = 5e-3;
clearance_x = 0.5e-3;
clearance_y = 0.25e-3;
clearance_z = 0.25e-3;
clearance_w = 0.25e-3;
% clearance_y = 4e-3;

sub_dom_xsize = 0.5e-3;
sub_dom_ysize = 4.25e-3;
sub_dom_zsize = 4.25e-3;
sub_dom_y_overlap_with_hole = 0.35e-3;
sub_dom_z_clearance = 0.25e-3;

trans_y_offset = hole_diam / 2; %How far off hole c/line
abs_bdry_thickness = 1e-3;
inc_mode = 'shear';
inc_angle_degs = 45;


show_section = 0;

%Details of input signal
centre_freq = 5e6;
no_cycles = 6;

%Component properties (Aluminium)
component_mat_i = 1;
component_c_L = 6300;
component_c_T = 3150;
component_density = 2700;

%Wedge properties (Acrylic)
wedge_mat_i = 2;
wedge_c_L = 2730;
wedge_c_T = 1345;
wedge_density = 1190;


solid_element_type = 'C3D8'; %C3D8 is an 8-noded brick

%--------------------------------------------------------------------------
switch inc_mode
    case 'shear'
        vel = component_c_T;
    case 'longitudinal'
        vel = component_c_L;
end

wedge_angle_degs = asind(wedge_c_L / vel * sind(inc_angle_degs));

trans_nom_pos = -thickness * tand(inc_angle_degs);
trans_x_pos = trans_nom_pos + params.trans_x_offset;%note this is actual position (at a -ve x value)

wedge_axial_len = trans_diam / 2 * tand(wedge_angle_degs) + clearance_w + abs_bdry_thickness;
wedge_radius = abs_bdry_thickness + trans_diam / 2;
xmin_trans = trans_x_pos - wedge_axial_len * sind(wedge_angle_degs) - wedge_radius * cosd(wedge_angle_degs);
xmin_comp = -abs_bdry_thickness - clearance_x - trans_diam / 2 / cosd(inc_angle_degs) + trans_nom_pos;
xmax_comp = hole_diam / 2 + clearance_x + abs_bdry_thickness;
ymax_comp = trans_y_offset + max(trans_diam / 2, sub_dom_ysize) + clearance_y + abs_bdry_thickness;
ymin_comp = trans_y_offset - trans_diam / 2 - clearance_y - abs_bdry_thickness;
wedge_max_thickness = wedge_axial_len * cosd(wedge_angle_degs) + wedge_radius * sind(wedge_angle_degs) + clearance_z;
crnr_pts = [
    min(xmin_trans, xmin_comp), ymin_comp, 0 %upper left back
    xmax_comp, ymax_comp, thickness + wedge_max_thickness]; %bottom right front

matls{component_mat_i} = fn_matl_isotropic_solid_defined_by_velocities('Aluminium', component_c_L, component_c_T, component_density);
matls{wedge_mat_i} = fn_matl_isotropic_solid_defined_by_velocities('Wedge', wedge_c_L, wedge_c_T, wedge_density, [255, 255, 224] / 256);

%Run model for long enough to see first echo with a bit spare
nom_path_len = 2 * sqrt(thickness ^ 2 + trans_x_pos ^ 2);
nom_arrival_time = nom_path_len / vel;
max_time = params.max_time_safety_factor * nom_arrival_time;

fe_options.pogo_path = 'C:\Program Files\Pogo\windows\new version';
fe_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';
fe_options.pogo_number_of_diff_absorbing_matls = 50;
fe_options.solver = 'pogo';
fe_options.dof_to_use = [1,2,3];


%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, params.els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);

el_types = fn_3d_el_types(); 

mod.el_typ_i = ones(size(mod.el_typ_i)) * find(strcmp(el_types, solid_element_type));
% mod.el_mat_i = ones(size(mod.el_typ_i));

%Identify which elements aer in wedge and component based on height
el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
component_els = el_ctrs(:, 3) <= thickness;
wedge_els = el_ctrs(:, 3) > thickness;

%Set element materials
mod.el_mat_i(component_els) = component_mat_i;
mod.el_mat_i(wedge_els) = wedge_mat_i;

%Absorbing layers for component
comp_abs_bdry_corner_nds = [
    xmin_comp + abs_bdry_thickness, ymin_comp + abs_bdry_thickness, 0         - clearance_z
    xmax_comp - abs_bdry_thickness, ymax_comp - abs_bdry_thickness, thickness + clearance_z
    ];
[comp_abs_bdry_nds, comp_abs_bdry_fcs] = fn_3d_rectalinear_surface(comp_abs_bdry_corner_nds(1,:), comp_abs_bdry_corner_nds(2,:));

mod = fn_3d_add_absorbing_layer(mod, comp_abs_bdry_nds, comp_abs_bdry_fcs, abs_bdry_thickness, component_els);
if ~params.absorption_in_component
    mod.el_abs_i(component_els) = 0;
end

%Identify elements of component to be keep LH end
component_els_in_use = component_els & (el_ctrs(:, 1) >= xmin_comp);

%Other cut-outs depending on whether it is corner calibration or hole
if params.corner_calibration_sample
    %Cut step out of lower half of model at x > 0
    component_els_in_use = component_els_in_use & ~(el_ctrs(:, 1) > 0 & el_ctrs(:, 3) < thickness / 2);
else
    %Make the hole
    component_els_in_use = component_els_in_use & (sqrt(sum((el_ctrs(:, 1:2) - [0, 0]) .^ 2, 2)) > hole_diam / 2);
end

%Define wedge
src_centre = [trans_x_pos, trans_y_offset, thickness];
dummy_pt1 = src_centre;
dummy_pt2 = src_centre;
trans_centre = src_centre;
wedge_path_len = clearance_z / cosd(wedge_angle_degs) + trans_diam / 2 * tand(wedge_angle_degs);
len_to_start_of_abs_layer = wedge_path_len + clearance_w;
dummy_pt1(1) = dummy_pt1(1) - sind(wedge_angle_degs) * len_to_start_of_abs_layer;
dummy_pt1(3) = dummy_pt1(3) + cosd(wedge_angle_degs) * len_to_start_of_abs_layer;
dummy_pt2(1) = dummy_pt2(1) + 2 * sind(wedge_angle_degs) * len_to_start_of_abs_layer; %dummy point on transducer axis deep in component
dummy_pt2(3) = dummy_pt2(3) - 2 * cosd(wedge_angle_degs) * len_to_start_of_abs_layer;

[wedge_abs_bdry_nds, wedge_abs_bdry_fcs] = fn_3d_cylindrical_surface(dummy_pt1, dummy_pt2, trans_diam / 2);
mod = fn_3d_add_absorbing_layer(mod, wedge_abs_bdry_nds, wedge_abs_bdry_fcs, abs_bdry_thickness, wedge_els);
wedge_els_in_use = wedge_els & (mod.el_abs_i < 1);
if ~params.absorption_in_wedge
    mod.el_abs_i(wedge_els) = 0;
end

in_use = component_els_in_use | wedge_els_in_use;
if show_section
    in_use = in_use & (el_ctrs(:,2) > trans_y_offset);
end

%Lose the elements and nodes that need to go
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
[mod.nds, mod.els, ~, ~] = fn_remove_unused_nodes(mod.nds, mod.els);

%Identify nodes on transducer face in wedge
trans_centre(1) = trans_centre(1) - sind(wedge_angle_degs) * wedge_path_len;
trans_centre(3) = trans_centre(3) + cosd(wedge_angle_degs) * wedge_path_len;
trans_centre_delay = wedge_path_len / wedge_c_L;

[bdry_nds, bdry_fcs] = fn_3d_disk_surface(trans_centre, trans_centre - dummy_pt2, trans_diam / 2);
d = fn_3d_signed_dist_to_bdry(mod.nds, bdry_nds, bdry_fcs);
tmp = find(abs(d) < el_size / 2);
steps{1}.load.frc_nds = [tmp(:); tmp(:)];
steps{1}.load.frc_dfs = [ones(size(tmp(:))); ones(size(tmp(:))) * 3]; %x and z DoF
steps{1}.load.wts = [ones(size(tmp(:))) * sind(wedge_angle_degs); -ones(size(tmp(:))) * cosd(wedge_angle_degs)];

% %create the subdomain
% subdomain_crns = [
%     -sub_dom_xsize / 2, -sub_dom_y_overlap_with_hole + hole_diam / 2,                 -sub_dom_z_clearance
%      sub_dom_xsize / 2, -sub_dom_y_overlap_with_hole + hole_diam / 2 + sub_dom_ysize,  sub_dom_zsize 
%     ];
% [inner_bndry_vtcs, inner_bndry_fcs] = fn_3d_rectalinear_surface(subdomain_crns(1, :), subdomain_crns(2, :));
% subdomain_mod = fn_create_subdomain(mod, el_types, inner_bndry_vtcs, inner_bndry_fcs, abs_bdry_thickness);

time_step = fn_get_suitable_time_step(matls, el_size);
time_pts = ceil(max_time / time_step);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq);

steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

end

