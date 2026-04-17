clear all
close all;

addpath(genpath('../code'));
show_geom_only = 0;

%ABOUT THIS EXAMPLE
%This example is designed to show a 3D model executed using the BristoFE solver
%(should be same result as ex_3d_pogo.m)

%PARAMETRIC DESCRIPTION OF MODEL
%The model is described in terms of a small number of parameters

%Material properties
solid_mat_i = 1;
solid_c_L = 6300;
solid_c_T = 3150;
solid_density = 2700;
matls{solid_mat_i} = fn_matl_isotropic_solid_defined_by_velocities('Aluminium', solid_c_L, solid_c_T, solid_density);

solid_element_type = 'C3D8'; %C3D8 is an 8-noded brick

%Define shape of model
model_size_x = 10e-3;
model_size_y = 11e-3;
model_size_z = 12e-3;

%Spherical void in the middle
scat_cent = [model_size_x, model_size_y, model_size_z] / 2;
scat_size = model_size_x / 10;

%Source is a circular disk of excitation in centre of top (max z) surface
src_radius = 3e-3;

abs_layer_thickness = 1e-3;
%corner bts
crnr_pts = [
    0, 0, 0
    model_size_x, model_size_y, model_size_z];

%Define a line along which sources will be placed to excite waves
src_centre = [0.5 * model_size_x, 0.5 * model_size_y, model_size_z];

src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;

%Run model for long enough to see first two echoes
max_time = 1.2 * (2 * model_size_z) / solid_c_L;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 6;

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
% mod.el_types = el_types;
el_types = fn_3d_el_types(); 

mod.el_typ_i = ones(size(mod.el_typ_i)) * find(strcmp(el_types, solid_element_type));
mod.el_mat_i = ones(size(mod.el_typ_i)) * solid_mat_i;

[scat_vtcs, scat_fcs] =  fn_3d_spherical_surface(scat_cent, scat_size /2);
scat_matl_i = 0;
scat_el_typ_i = [];
mod = fn_3d_add_inclusion_or_void(mod, el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i);


% %Quick test of absorbing layers
% el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
% mod.el_abs_i = zeros(size(mod.el_mat_i));
% tmp = (abs_layer_thickness - el_ctrs(:,1)) / abs_layer_thickness;
% mod.el_abs_i = tmp .* (tmp > 0);

%Quick test of making a void
% i = sqrt(sum((el_ctrs - [0.5, 0.5, 0] * model_size) .^ 2, 2)) < 2e-3;
% mod.els(i, :) = [];
% mod.el_mat_i(i) = [];
% mod.el_abs_i(i) = [];
% [mod.nds, mod.els, ~, ~] = fn_remove_unused_nodes(mod.nds, mod.els);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_node_nearest_to_point(mod.nds, src_centre, el_size);
steps{1}.load.frc_nds = find(...
    abs(mod.nds(:, 3) - src_centre(3)) < el_size / 2 & ...
    sqrt(sum( (mod.nds(:, 1:2) - src_centre(1:2)) .^ 2, 2 )) < src_radius ...
    );


steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
if show_geom_only %suppress graphics when running all scripts for testing
    figure;
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end
%--------------------------------------------------------------------------
%RUN THE MODEL
fe_options = [];
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, -sum(res{1}.dsps, 1));
xlabel('Time (s)')

