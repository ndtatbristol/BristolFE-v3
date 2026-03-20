clear
close all
addpath(genpath('../code'));

%Simple mesh
model_size = 10;
el_size = 0.2;
cod = el_size / 10;

host_matl_i = 1;
matls{host_matl_i} = fn_matl_isotropic_solid_defined_by_velocities('Al', 6300, 3150, 2700);

scat_matl_i = 2
matls{scat_matl_i} = fn_matl_isotropic_solid_defined_by_stiffness('Steel', 210e9, 0.3, 8900);

el_typ_host = 'CPE3'; 
el_typ_scat = 'CPE3'; 

%Create the mesh
bdry_pts = [
    0,  0
    0,  1
    1, 1
    1, 0] *model_size;
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
el_types = fn_2d_el_types();
mod.el_typ_i(:) = find(strcmp(el_typ_host, el_types));
mod.el_mat_i(:) = host_matl_i;

%CRACKS
%Defined crack vertices - a nice branched crack with some roughness
crack_pts1 = 20;
crack_len1 = model_size * 0.5;
crack_pts2 = 10;
crack_len2 = model_size * 0.25;
crack_start1 = [-el_size,  model_size * 0.5];
%top part
crack_vtcs1 = crack_start1 + fn_2d_random_walk(crack_pts1, crack_len1/crack_pts1, 0, 0, 0, 0.4);
%branch start (mid way along top branch)
crack_start2 = crack_vtcs1(round(crack_pts1 / 2), :); 
crack_vtcs2 = crack_start2 + fn_2d_random_walk(crack_pts2, crack_len2/crack_pts2, 0, 0, 0, 0.4);
%Add the cracks
mod = fn_2d_add_crack(mod, el_types, crack_vtcs1, [], cod);
mod = fn_2d_add_crack(mod, el_types, crack_vtcs2, [], cod);

%INCLUSIONS
scat_centre1 = [model_size * 0.6, model_size * 0.75];
scat_size1 = model_size * 0.1;
scat_el_typ_i1 = find(strcmp(el_typ_scat, el_types));
scat_min_rad_frac1 = 0.5;
scat_complexity1 = 3;
scat_pts1 = 200;
scat_vtcs1 = fn_2d_create_smooth_random_blob(scat_min_rad_frac1, scat_complexity1, scat_pts1) * scat_size1 / 2 + scat_centre1;
scat_matl_i1 = 0;

scat_centre2 = [model_size * 0.6, model_size * 0.25];
scat_size2 = model_size * 0.2;
scat_el_typ_i2 = find(strcmp(el_typ_scat, el_types));
scat_min_rad_frac2 = 0.2;
scat_complexity2 = 4;
scat_pts2 = 200;
scat_vtcs2 = fn_2d_create_smooth_random_blob(scat_min_rad_frac2, scat_complexity2, scat_pts2) * scat_size2 / 2 + scat_centre2;
scat_matl_i2 = scat_matl_i;

mod = fn_2d_add_inclusion_or_void(mod, el_types, scat_vtcs1, scat_matl_i1, scat_el_typ_i1);
mod = fn_2d_add_inclusion_or_void(mod, el_types, scat_vtcs2, scat_matl_i2, scat_el_typ_i2);

%Plot result
figure;
options.draw_elements = 1;
fn_show_geometry(mod, matls, fn_2d_el_types(), options);
hold on;
plot(crack_vtcs1(:, 1), crack_vtcs1(:, 2), 'r:')
plot(crack_vtcs2(:, 1), crack_vtcs2(:, 2), 'r:')
