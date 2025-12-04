clear all
close all

addpath(genpath('..\code'));

fe_options.pogo_path = 'C:\Program Files\Pogo\windows\new version';
fe_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

show_geom_only = 1;

abs_bdry_thickness = 1e-3;

%Material properties
solid_matl_i = 1;
main.matls(solid_matl_i) = fn_matl_isotropic_solid_defined_by_stiffness('Steel', 210e9, 0.3, 8900);

% main.el_types = fn_3d_el_types(); %C3D8 8 noded brick
solid_element_type = 'C3D8';

%Define shape of model
model_size_x = 10e-3;
model_size_y = 10e-3;
model_size_z = 12e-3;

src_radius = 3e-3;

abs_layer_thickness = 1e-3;
%corner bts
crnr_pts = [
    0, 0, 0
    model_size_x, model_size_y, model_size_z];

subdom_centre = [model_size_x, model_size_y, model_size_z] / 2;
subdom_rad = 1e-3;

%Define a line along which sources will be placed to excite waves
src_centre = [0.5 * model_size_x, 0.5 * model_size_y, model_size_z];

src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 8;

fe_options.solver = 'pogo';
%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
main.mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
main.el_types = fn_3d_el_types()
main.mod.el_typ_i = ones(size(main.mod.el_typ_i)) * find(strcmp(main.el_types, solid_element_type));
main.mod.el_mat_i = ones(size(main.mod.el_typ_i)) * solid_matl_i;

%Quick test of absorbing layers
% el_ctrs = fn_calc_element_centres(main.mod.nds, main.mod.els);
% main.mod.el_abs_i = zeros(size(main.mod.el_mat_i));
% tmp = (abs_layer_thickness - el_ctrs(:,1)) / abs_layer_thickness;
% main.mod.el_abs_i = tmp .* (tmp > 0);


%Quick test of making a void
% i = sqrt(sum((el_ctrs - [0.5, 0.5, 0] * model_size) .^ 2, 2)) < 2e-3;
% mod.els(i, :) = [];
% mod.el_mat_i(i) = [];
% mod.el_abs_i(i) = [];
% [mod.nds, mod.els, ~, ~] = fn_remove_unused_nodes(mod.nds, mod.els);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
main.trans{1}.nds = find(...
    abs(main.mod.nds(:, 3) - src_centre(3)) < el_size / 2 & ...
    sqrt(sum( (main.mod.nds(:, 1:2) - src_centre(1:2)) .^ 2, 2 )) < src_radius ...
    );
main.trans{1}.dfs = ones(size(main.trans{1}.nds)) * src_dir;


%Create the subdomain

[inner_bndry_vtcs, inner_bndry_fcs] = fn_3d_spherical_surface(subdom_centre, subdom_rad);
main.doms{1}.mod = fn_3d_create_subdomain(main.mod, inner_bndry_vtcs, inner_bndry_fcs, abs_bdry_thickness);

% main.doms{1}.mod = fn_2d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_pts, 0);


% %Also provide the time signal for the loading (if this is a vector, it will
% %be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
% %of different time signals for each frc_nds/frc_dfs
% time_step = fn_get_suitable_time_step(matls, el_size);
% steps{1}.load.time = 0: time_step:  max_time;
% steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);
% 
% %Also record displacement history at same points (NB there is no reason why
% %these have to be same as forcing points)
% steps{1}.mon.nds = steps{1}.load.frc_nds;
% steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
if show_geom_only %suppress graphics when running all scripts for testing
    figure;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry_with_subdomains(main, display_options);
    drawnow
    % figure;
    % display_options.transparency = 0.5;
    % display_options.draw_elements = 0;
    % display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    % display_options.node_sets_to_plot(1).col = 'r.';
    % h_patch = fn_show_geometry(mod, matls, display_options);
    return
end
%--------------------------------------------------------------------------
%RUN THE MODEL

% [res, mats] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps, 1));
xlabel('Time (s)')

