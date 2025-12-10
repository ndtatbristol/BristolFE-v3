clear all
close all

addpath(genpath('..\code'));
addpath(genpath('..\subdoms'));

fe_options.pogo_path = 'C:\Program Files\Pogo\windows\new version';
fe_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

show_geom_only = 0;

abs_bdry_thickness = 1e-3;

%Material properties
solid_matl_i = 1;
main.matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities('Aluminium', 6300, 3150, 2700);

% main.el_types = fn_3d_el_types(); %C3D8 8 noded brick
solid_element_type = 'C3D8';

%Define shape of model
model_size_x = 10e-3;
model_size_y = 10e-3;
model_size_z = 12e-3;

src_radius = 3e-3;

%corner bts
crnr_pts = [
    0, 0, 0
    model_size_x, model_size_y, model_size_z];

subdom_centre = [model_size_x, model_size_y, model_size_z] / 2;
subdom_rad = 1e-3;

scat_cent = subdom_centre;
scat_rad = subdom_rad / 2;

%Define a line along which sources will be placed to excite waves
src_centre = [0.5 * model_size_x, 0.5 * model_size_y, model_size_z];

src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
fe_options.number_of_cycles = 5;
fe_options.max_time = 1.5 * 2 * model_size_z / 6300;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 5;

fe_options.solver = 'pogo';
fe_options.dof_to_use = [1,2,3];
%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
main.mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
main.mod.design_centre_freq = centre_freq;
main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size);

main.el_types = fn_3d_el_types();
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
main.doms{1}.mod = fn_3d_create_subdomain(main.mod, main.el_types, inner_bndry_vtcs, inner_bndry_fcs, abs_bdry_thickness);

[scat_vtcs, scat_fcs] =  fn_3d_spherical_surface(scat_cent, scat_rad);
scat_matl_i = 0;
scat_el_typ_i = [];
main.doms{1}.mod = fn_3d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_vtcs, scat_fcs, scat_matl_i, scat_el_typ_i);


%Show the mesh
if show_geom_only %suppress graphics when running all scripts for testing
    figure;
    subplot(1,2,1);
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(main.mod, main.matls, main.el_types, display_options);
    aa = axis;
    title('Main')
    subplot(1,2,2);
    display_options = [];
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    % display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    % display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(main.doms{1}.mod, main.matls, main.el_types, display_options);
    title('Subdomain')
    axis(aa);
    return
end
%--------------------------------------------------------------------------
%RUN THE MODEL
%Run main model
main = fn_run_main_model(main, fe_options);

%Run sub-domain model
main = fn_run_subdomain_model(main, fe_options);

%Run validation model
fe_options.validation_mode = 1;
main = fn_run_main_model(main, fe_options);

figure;
i = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
mv = max(abs(sum(main.doms{1}.res.fmc.time_data(i:end,: ), 2)));
plot(main.doms{1}.res.fmc.time, real(sum(main.doms{1}.res.fmc.time_data, 2)) / mv, 'k', 'LineWidth', 2);
hold on;
plot(main.doms{1}.val.fmc.time, real(sum(main.doms{1}.val.fmc.time_data, 2)) / mv, 'g:', 'LineWidth', 2);
plot(main.res.fmc.time, real(sum(main.res.fmc.time_data, 2)) / mv, 'b');
ylim([-1,1]);
yyaxis right
plot(main.doms{1}.res.fmc.time, 20 * log10(abs(sum(main.doms{1}.res.fmc.time_data, 2) - sum(main.doms{1}.val.fmc.time_data, 2)) / mv));
ylim([-60, 0]);
legend('Sub-domain method', 'Validation', 'Pristine', 'Difference (dB)');

