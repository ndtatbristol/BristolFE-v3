clear
close all
addpath(genpath('../code'));

%Simple mesh
bdry_pts = [
    0,  0
    0,  10
    10, 10
    10, 0];
el_size = 0.1;
cod = el_size / 10;

steel_mat_i = 1;
matls{steel_mat_i}.rho = 8900; %Density
matls{steel_mat_i}.D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls{steel_mat_i}.col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls{steel_mat_i}.name = 'Steel';

el_typ_solid = 'CPE3'; 

%Defined crack vertices - a nice branched crack with some roughness
pts1 = 20;
len1 = 5;
pts2 = 10;
len2 = 2.5;

min_rad_frac = 0.5;
complexity = 3;
no_pts = 200;
scat_matl_i = 0;
scat_el_typ_i = 0;
blob_rad = 3;
blob_centre = [5,5];


%Create the mesh
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
el_types = {el_typ_solid};
mod.el_typ_i(:) = find(strcmp(el_typ_solid, el_types));
mod.el_mat_i(:) = steel_mat_i;

%Add the inclusion
pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * blob_rad + blob_centre;
mod = fn_2d_add_inclusion_or_void(mod, el_types, pts, scat_matl_i, scat_el_typ_i);

%Plot result
figure;
options.draw_elements = 1;
fn_show_geometry(mod, matls, options);
