function [main, fe_options, params] = mod_2d_subdomain_oblique(params)

%This uses mod_2d_oblique without either of the scatterers) as the 
%pristine model. The crack and reference notch that would normally defined 
%in mod_2d_oblique are then put into two sub-domain models

default_params.fe_options.field_output_every_n_frames = 20;
%--------------------------------------------------------------------------
params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
params = fn_set_default_fields(params, default_params);

params.include_crack = 0;

[main.mod, main.matls, main.el_types, steps, fe_options, params] = mod_2d_oblique(params);
params.fe_options = fe_options;
main.trans{1}.nds = steps{1}.load.frc_nds;
main.trans{1}.dfs = steps{1}.load.frc_dfs;

%Create a subdomain based on params.defect_region_size
inner_bdry = [-1,0;-1,2;1,2;1,0] / 2 .* (params.defect_region_size + 2 * params.el_size);

%First subdomain left empty to show what happens
main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);

%Second subdomain contains crack
main.doms{2}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);
cod = params.el_size / 10;
main.doms{2}.mod = fn_2d_add_crack(main.doms{2}.mod, main.el_types, params.crack_vtcs, [], cod);

%Third subdomain contains reference notch
main.doms{3}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);
main.doms{3}.mod = fn_2d_add_inclusion_or_void(main.doms{3}.mod, main.el_types, params.notch_pts, 0, 0);

%Input signal
main.inp.time = steps{1}.load.time;
main.inp.sig = steps{1}.load.frcs;
end

