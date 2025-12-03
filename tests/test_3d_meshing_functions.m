clear
close all
addpath(genpath('../code'));

%Simple mesh
model_size = 10;
crnr_pts = [
    0, 0, 0
    1, 1, 1] * model_size;
el_size = 0.2;



%Materials
solid_matl_i = 1;
matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities('Al', 6300, 3150, 2700);
solid_matl2_i = 2;
matls{solid_matl2_i} = fn_matl_isotropic_solid_defined_by_velocities('St', 5900, 3150, 8900);
fluid_matl_i = 3;
matls{fluid_matl_i} = fn_matl_fluid_defined_by_velocity('water', 1500, 1000);

el_typ_solid = 'C3D8R'; %C3D8 8 noded brick
el_typ_fluid = 'AC3D8'; %C3D8 8 noded brick
el_types = fn_3d_el_types();

mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
mod.el_typ_i(:) = find(strcmp(el_typ_solid, el_types));
mod.el_mat_i(:) = solid_matl_i;

cent = [0,0,0];
rad = model_size / 2;
n_sub_divisions = 1;
[vtcs, fcs] = fn_3d_spherical_surface(cent, rad, n_sub_divisions);
els_in_mat2 = fn_3d_find_elements_in_region(mod, vtcs, fcs);
mod.el_mat_i(els_in_mat2) = solid_matl2_i;
mod.el_typ_i(els_in_mat2) = find(strcmp(el_typ_solid, el_types));

cent = [model_size,0,0];
rad = model_size / 3;
n_sub_divisions = 2;
[vtcs, fcs] = fn_3d_spherical_surface(cent, rad, n_sub_divisions);
els_in_fluid = fn_3d_find_elements_in_region(mod, vtcs, fcs);
mod.el_mat_i(els_in_fluid) = fluid_matl_i;
mod.el_typ_i(els_in_fluid) = find(strcmp(el_typ_fluid, el_types));

figure;
display_options.transparency = 0.5;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, el_types, display_options);

axis on;
xlabel('x')
ylabel('y')
zlabel('z')