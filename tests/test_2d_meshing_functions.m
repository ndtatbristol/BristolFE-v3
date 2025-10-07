clear
close all
%restoredefaultpath

addpath(genpath('../code'));

test_basic_mesh_gen = 0;

%Steel
steel_matl_i = 1;
matls{steel_matl_i}.rho = 8900; %Density
matls{steel_matl_i}.D = fn_isotropic_stiffness_matrix(210e9, 0.3);
matls{steel_matl_i}.col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls{steel_matl_i}.name = 'Steel';

%Water
water_matl_i = 2;
matls{water_matl_i}.rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calculated here from ultrasonic velocity (1500) and density (1000)
matls{water_matl_i}.D = 1500 ^ 2 * 1000;
matls{water_matl_i}.col = hsv2rgb([0.6,0.5,0.8]);
matls{water_matl_i}.name = 'Water';

el_typ_solid = 'CPE3';
el_typ_fluid = 'AC2D3';
el_types = {el_typ_solid, el_typ_fluid};
bdry_pts = [0.02,0.01; 0.02, 0.7; 1.03,0.7; 1.03, 0.01];
water_pts = [0.5, 0; 0.6, 1; 2, 1; 2, 0];
options.draw_elements = 1;

%--------------------------------------------------------------------------
%Basic mesh generation
if test_basic_mesh_gen
    el_size = 0.1;

    figure;
    force_reg_elements = 1;
    subplot(2,2,1);
    mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, force_reg_elements);
    fn_show_geometry(mod, matls, options);
    hold on; fn_plot_line(bdry_pts, 'r', 1)
    subplot(2,2,2);
    mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size, force_reg_elements);
    fn_show_geometry(mod, matls, options);
    hold on; fn_plot_line(bdry_pts, 'r', 1)
    force_reg_elements = 0;
    subplot(2,2,3);
    mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, force_reg_elements);
    fn_show_geometry(mod, matls, options);
    hold on; fn_plot_line(bdry_pts, 'r', 1)
    subplot(2,2,4);
    mod2 = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size, force_reg_elements);
    fn_show_geometry(mod2, matls, options);
    hold on; fn_plot_line(bdry_pts, 'r', 1)
end
%--------------------------------------------------------------------------
%Change some material, add fluid-solid interface, and add voids
el_size = 0.02;
force_reg_elements = 0;
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, force_reg_elements);

%Add water region and deal with element types
els_in_water = fn_elements_in_region(mod, water_pts);
mod.el_typ_i = ones(size(mod.els, 1), 1) * find(strcmp(el_types, el_typ_solid));
mod.el_mat_i(els_in_water) = water_matl_i;
mod.el_typ_i(els_in_water) = find(strcmp(el_types, el_typ_fluid));


rad = 0.2;
min_rad_frac = 0.5;
complexity = 3;
no_pts = 200;
cent = [0.5, 0.35];
scat_pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * rad + cent;
scat_matl = 0;
[mod, el_types] = fn_add_fluid_solid_interface_els(mod, el_types);
mod = fn_2d_add_inclusion_or_void(mod, matls, el_types, scat_pts, scat_matl);

figure;
options.draw_elements = 0;
fn_show_geometry(mod, matls, options);
hold on; fn_plot_line(bdry_pts, 'r', 1)
hold on; fn_plot_line(scat_pts, 'r', 1);
