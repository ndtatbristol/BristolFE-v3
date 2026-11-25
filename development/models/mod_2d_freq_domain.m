function [mod, matls, el_types, steps, fe_options, params] = mod_2d_freq_domain(params)
%Simple 2D freq domain model of embedded scatterer

default_params.incident_angs_degs = [0, 10];
default_params.freq = 1e6;

default_params.scatterer_size = 1e-3;
default_params.clearance = 1e-3;
default_params.abs_bdry_thickness_in_wavelengths = 1;

default_params.els_per_wavelength = 10;

default_params.fluid_matl_name = 'water';

default_params.include_scatterer = 1;
default_params.scatterer_matl_name = 'glycerin';
default_params.scatterer_shape = 'blob';

%Element type to use
default_params.element_shape = 'tri';

%Used for time-domain simulations only - kept in to allow comparison
default_params.safety_factor = 1.5;

default_params.random_seed = 1;

%--------------------------------------------------------------------------
params = fn_set_default_fields(params, default_params);
fe_options = fn_set_fe_options_from_params(params);

ang_pts = 361;
ang_degs = linspace(0, 360, ang_pts)';
cos_sin = [cosd(ang_degs), sind(ang_degs)];

rng(params.random_seed); 
el_types = fn_2d_el_types();

%Define the materials in use
fluid_matl_i = 1;
matls{fluid_matl_i} = fn_material_library(params.fluid_matl_name);

scat_matl_i = 2;
matls{scat_matl_i} = fn_material_library(params.scatterer_matl_name);

%El size and time step
params.el_size = fn_get_suitable_el_size(matls, max(params.freq), params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, params.el_size, params.safety_factor);

switch params.element_shape
    case {'tri', 'triangular'}
        el_typ_to_use_for_solid = 'CPE3';
        el_typ_to_use_for_fluid = 'AC2D3';
    case {'quad', 'rect', 'quadrilateral', 'rectangle', 'rectangular'}
        el_typ_to_use_for_solid = 'CPE4';
        el_typ_to_use_for_fluid = 'AC2D4';
    otherwise
        error('Unknown element shape')
end

switch fn_matl_type(matls{scat_matl_i})
    case 'solid'
        scat_el_typ = el_typ_to_use_for_solid;
    case 'fluid'
        scat_el_typ = el_typ_to_use_for_fluid;
end
scat_el_typ_i = find(strcmp(el_types, scat_el_typ));

max_wavelength = fn_estimate_max_min_wavelengths(matls, min(params.freq));
params.abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;

model_rad = params.scatterer_size + 1.5 * params.clearance + params.abs_bdry_thickness;

bdry_pts = cos_sin * model_rad;

mod = fn_2d_structured_mesh(bdry_pts, params.el_size, el_typ_to_use_for_fluid);
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_fluid));
mod.el_mat_i(:) = fluid_matl_i;

%absorbing layer
params.abs_bdry_pts = cos_sin * (model_rad - params.abs_bdry_thickness);
mod = fn_2d_add_absorbing_layer(mod, params.abs_bdry_pts, params.abs_bdry_thickness);
[fe_options.max_damping, fe_options.damping_power_law , fe_options.max_stiffness_reduction]= fn_optimum_absorbing_bdry_properties(params.abs_bdry_thickness, matls, params.freq);


if params.include_scatterer
    switch params.scatterer_shape
        case 'blob'
            min_rad_frac = 0.5;
            complexity = 3;
            no_pts = 200;
            params.scatterer_pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts) * params.scatterer_size;
        case 'void'
            params.scatterer_pts = cos_sin * params.scatterer_size;
            scat_matl_i = 0;
    end
    % scat_matl_i = 0;
    mod = fn_2d_add_inclusion_or_void(mod, el_types, params.scatterer_pts, scat_matl_i, scat_el_typ_i);
end

%get the boundary nodes using subdomain functions
inner_bdry_vtcs = cos_sin * (model_rad - params.abs_bdry_thickness);
els_in_region = fn_2d_find_elements_in_region(mod, inner_bdry_vtcs);
[mod.bdry_lyrs, el_used] = fn_subdomain_bdry_layers(mod.nds, mod.els, els_in_region);

%should do better subdomain implementation so it starts outside these bdry
%nodes - see line 45 etc in fn_subdomain_bdry_layers

params.bdry_nd_indices = find(mod.bdry_lyrs);
params.fluid_vel = fn_vels_from_matl_def(matls{fluid_matl_i});
% for s = 1:numel(params.incident_angs_degs)
% 
%     steps{s}.load.frc_nds = bdry_nd_indices;
% 
%     steps{s}.load.freq = params.freq;
%     for f = 1:numel(params.freq)
%         % steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * transducer_dof;
%         % steps{1}.load.frc_wts = ones(1, numel(steps{1}.load.frc_nds));
%         % steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);
% 
%         %Also record displacement history at same points (NB there is no reason why
%         %these have to be same as forcing points)
%         steps{s}.mon.dsp_nds = steps{1}.load.frc_nds;
%         % steps{1}.mon.dsp_dfs = steps{1}.load.frc_dfs;
%         % steps{1}.mon.dsp_wts = steps{1}.load.frc_wts;
%     end
% end
steps = [];

end


