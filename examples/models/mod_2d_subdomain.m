function [main, fe_options, params] = mod_2d_subdomain(params)

%This uses mod_2d_advanced (without either of the scatterers) as the 
%pristine model. The crack that would normally defined in mod_2d_advanced
%is then put in a sub-domain model.

default_params.fe_options_field_output_every_n_frames = 20;
%--------------------------------------------------------------------------
params = fn_set_default_fields(params, default_params);
fe_options = fn_set_fe_options_from_params(params);


params.include_crack = 0;
params.include_scatterer = 0;

[main.mod, main.matls, main.el_types, steps, fe_options, params] = mod_2d_advanced(params);
main.trans{1}.nds = steps{1}.load.frc_nds;
main.trans{1}.dfs = steps{1}.load.frc_dfs;

%Create a subdomain in the middle with a hole in surface as scatterer
scatterer_centre = (max(params.crack_vtcs) + min(params.crack_vtcs)) / 2;
min_dims = max(params.crack_vtcs) - min(params.crack_vtcs);
inner_bdry = [-1,-1;-1,1;1,1;1,-1] / 2 * (max(min_dims) + 2 * params.el_size) + scatterer_centre;

%First subdomain left empty to show what happens
main.doms{1}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);

%Second subdomain contains crack
main.doms{2}.mod = fn_2d_create_subdomain(main.mod, main.el_types, inner_bdry, params.abs_bdry_thickness);
cod = params.el_size / 10;
main.doms{2}.mod = fn_2d_add_crack(main.doms{2}.mod, main.el_types, params.crack_vtcs, [], cod);

%Input signal
main.inp.time = steps{1}.load.time;
main.inp.sig = steps{1}.load.frcs;
end


