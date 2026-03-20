clear
close all

addpath(genpath('../code'));
addpath(genpath('../subdoms'));


%test of reinsertion of sub-domain model into pristing model

%Overall model geometry
xsize = 10e-3;
ysize = 10e-3;
src_size = 3.5e-3;
abs_bdry_thickness = 1e-3;
el_size = 0.1e-3;

%Subdomain geometry
subdomain_size = 3e-3;

%Solid material properties
plate_matl_name = 'steel';
plate_matl_longitudinal_velocity = 5636;
plate_matl_shear_velocity = 3012;
plate_matl_density = 8900;

%Element types to use
el_typ_to_use_for_solid = 'CPE3';
% el_typ_to_use_for_solid = 'CPE4';

%--------------------------------------------------------------------------
%END OF INPUTS

%Define the materials
%A cell array with an entry for each material used in the model is required.
plate_matl_i = 1;
main.matls{plate_matl_i} = fn_matl_isotropic_solid_defined_by_velocities(plate_matl_name, plate_matl_longitudinal_velocity, plate_matl_shear_velocity, plate_matl_density);

%Define shape of model
xmin =0; xmax = xsize;
ymin =0; ymax = ysize;
bdry_pts = [
    -1, -1
    -1,  1
    1,  1
    1, -1] .* [xsize, ysize] / 2;

%absorbing boundary for plate
abs_bdry_pts_plate = [
    -1, -1
    -1,  1
    1,  1
    1, -1] .* ([xsize, ysize] / 2 - [1, 1] * abs_bdry_thickness);

switch el_typ_to_use_for_solid
    case 'CPE3'
        main.mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
    case 'CPE4'
        main.mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size);
end

main.el_types = fn_2d_el_types();

main.mod.el_typ_i(:) = find(strcmp(main.el_types, el_typ_to_use_for_solid));

%Set materials
main.mod.el_mat_i(:) = plate_matl_i;



%Define the absorbing layer
main.mod = fn_2d_add_absorbing_layer(main.mod, abs_bdry_pts_plate, abs_bdry_thickness);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, main.matls, 5e6);


%Lose the elements that need to go
el_ctrs = fn_calc_element_centres(main.mod.nds, main.mod.els);
% els_to_go = (el_ctrs(:,1) - el_ctrs(:,2)) > 0 & el_ctrs(:,2) < ysize * 0.5;
els_to_go = el_ctrs(:,1) < 0 & el_ctrs(:,2) < 0;
% els_to_go = [];
main.mod.els(els_to_go, :) = [];
main.mod.el_mat_i(els_to_go) = [];
main.mod.el_abs_i(els_to_go) = [];
main.mod.el_typ_i(els_to_go) = [];
[main.mod.nds, main.mod.els, ~, ~] = fn_remove_unused_nodes(main.mod.nds, main.mod.els);

src_end_pts = [- src_size / 2, ysize / 2 - abs_bdry_thickness
               + src_size / 2, ysize / 2 - abs_bdry_thickness];
[main.trans{1}.nds, s] = fn_find_nodes_nearest_to_line(main.mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
main.trans{1}.dfs = ones(size(main.trans{1}.nds)) * 1;


%Create a subdomain in the middle with a hole in surface as scatterer
inner_bdry = [
    -1, -1
    -1,  1
    1,  1
    1, -1] .* subdomain_size / 2;

scatterer_size = subdomain_size;
scat_pts =   fn_2d_create_smooth_random_blob(0.4, 3, 360) * scatterer_size / 2;
main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, abs_bdry_thickness);
main.doms{1}.mod = fn_2d_add_inclusion_or_void(main.doms{1}.mod, main.el_types, scat_pts, 0, 0);

[val_mod, old_nds, new_nds] = fn_insert_subdomain_model_into_main(main.mod, main.doms{1}.mod);

%Show the mesh
figure;
subplot(2,2,1);
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(main.mod, main.matls, main.el_types, display_options);
cols = 'rgbm';
for i = 1:4
    b = main.doms{1}.mod.main_nd_i(main.doms{1}.mod.bdry_lyrs == i);
    plot(main.mod.nds(b,1),main.mod.nds(b,2), [cols(i), '.']);
end

subplot(2,2,2);
display_options.node_sets_to_plot = [];
h_patch = fn_show_geometry(main.doms{1}.mod, main.matls, main.el_types, display_options);
cols = 'rgbm';
for i = 1:4
    b = main.doms{1}.mod.main_nd_i(main.doms{1}.mod.bdry_lyrs == i);
    plot(main.mod.nds(b,1),main.mod.nds(b,2), [cols(i), '.']);
end

subplot(2,2,3);
display_options.node_sets_to_plot = [];
h_patch = fn_show_geometry(val_mod, main.matls, main.el_types, display_options);

