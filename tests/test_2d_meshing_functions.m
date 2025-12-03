clear
close all

addpath(genpath('../code'));

%Materials
solid_matl_i = 1;
matls{solid_matl_i} = fn_matl_isotropic_solid_defined_by_velocities('Al', 6300, 3150, 2700);
solid_matl2_i = 2;
matls{solid_matl2_i} = fn_matl_isotropic_solid_defined_by_velocities('St', 5900, 3150, 8900);
fluid_matl_i = 3;
matls{fluid_matl_i} = fn_matl_fluid_defined_by_velocity('water', 1500, 1000);

el_typ_solid = 'CPE3';
el_typ_fluid = 'AC2D3';
el_types = fn_2d_el_types();
abs_bdry_thickness = 0.05;
bdry_pts = [
    0.02, 0.01
    0.02, 0.7
    1.03, 0.7
    1.03, 0.01];
abs_bdry_pts = [
    0.02 + abs_bdry_thickness, 0.01 + abs_bdry_thickness
    0.02 + abs_bdry_thickness, 0.7 - abs_bdry_thickness
    1.03 - abs_bdry_thickness, 0.7 - abs_bdry_thickness
    1.03 - abs_bdry_thickness, 0.01 + abs_bdry_thickness];

water_pts = [0.5, 0; 0.6, 1; 2, 1; 2, 0];
el_size = 0.1;
options.draw_elements = 1;

test_basic_mesh_gen = 0;
test_advance_mesh_gen = 1;

%--------------------------------------------------------------------------
%Basic mesh generation
if test_basic_mesh_gen


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
if test_advance_mesh_gen
    rad = 0.1;
    min_rad_frac = 0.5;
    complexity = 3;
    no_pts = 200;
    cent = [0.5, 0.25];
    scat_pts1 = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * rad + cent;
    scat_matl1 = 0;

    cent = [0.6, 0.5];
    scat_pts2 = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * rad + cent;
    scat_matl2 = solid_matl2_i;


    for i = 1:2

        %Change some material, add fluid-solid interface, and add voids
        el_size = 0.02;
        force_reg_elements = 0;
        switch i
            case 1
                mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, force_reg_elements);
                el_typ_solid = 'CPE3';
                el_typ_fluid = 'AC2D3';
            case 2
                mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size, force_reg_elements);
                el_typ_solid = 'CPE4';
                el_typ_fluid = 'AC2D4';
        end


        %Add water region and deal with element types
        els_in_water = fn_2d_find_elements_in_region(mod, water_pts);
        mod.el_typ_i = ones(size(mod.els, 1), 1) * find(strcmp(el_types, el_typ_solid));
        mod.el_mat_i(els_in_water) = fluid_matl_i;
        mod.el_typ_i(els_in_water) = find(strcmp(el_types, el_typ_fluid));


        mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

        %Add some scatterers (NB the add fluid-solid interface is called
        %automatically in these)
        scat_el_typ1 = 0;
        mod = fn_2d_add_inclusion_or_void(mod, el_types, scat_pts1, scat_matl1, scat_el_typ1);
        scat_el_typ2 = find(strcmp(el_types, el_typ_solid));
        mod = fn_2d_add_inclusion_or_void(mod, el_types, scat_pts2, scat_matl2, scat_el_typ2);

        figure;
        fn_show_geometry(mod, matls, options);
        hold on; fn_plot_line(bdry_pts, 'r', 1)
        hold on; fn_plot_line(scat_pts1, 'r', 1);
        hold on; fn_plot_line(scat_pts2, 'r', 1);
    end
end