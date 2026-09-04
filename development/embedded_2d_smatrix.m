clear
close all

%This is work in progress!!

%Basic idea (see note below for where it's at):

%Run defect-free sub-domain model for each incident direction to get
%transfer functions between boundary displacemenets and incident/scattered
%plane waves. Then run defect model for each incident direction (probably
%subtracted defect-free scattered results first (see note at end) and 
%project measured scattered field back to each scattering direction using
%reciprocity. There need to be a bit of Fourier and phase shifting action
%and probably some scaling to get freq-domain results in exact same form as
%usual S-matrix definition.

%30/6/26 - incident wave works ok if toneburst used, but there is still 
%some leakage due to analytical model of incident wave not matching FE 
%exactly. In case with free surface this is always going to be case, so 
%probably solution is to always run defect-free case and subtract result of
%this from any defect ones.

%The projection back to scattered directions is not yet done. Probably this
%can be achieved by adapting the subdomain code as that deals with all the
%reciprocity stuff. You would end up with FMC-like data where each
%transmitter is an incident plane wave angle and each receiver is a plane
%wave receiver angle.

model_to_run = @mod_2d_embedded_smatrix;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 20;

params.surface_model = 0; %test to see what heppens with free surface in model

%If you just want to see the model (without running it, set show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'subdoms']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model - subdomain code is used for this in order to get the
%nodes around the edge
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

%Show the subdomain mesh and stop if requested
col = 'rgbc';
for i = 1:4
    display_options.node_sets_to_plot(i).nd = find(mod.bdry_lyrs == i);
    display_options.node_sets_to_plot(i).col = [col(i), '.'];
end
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

%compute displacements on bdry nodes for certain incident wave
params.wave_start_distance = params.wave_start_distance * 1.5;
inc_mode = 'L';
inc_angle = 30 * pi / 180;
p = [params.wave_start_distance,  params.wave_start_distance
     params.wave_start_distance, -params.wave_start_distance];
p = [cos(inc_angle), -sin(inc_angle)
     sin(inc_angle),  cos(inc_angle)] * p;

bdry_nds = find(mod.bdry_lyrs ~= 0);
bdry_dists = fn_dist_point_to_line(mod.nds(bdry_nds,:), p(:,1)', p(:, 2)');

%Function to approximate dirac delta at arbitrary time by vector at discrete times
fn_discrete_delta = @(t, td) (1 - abs(t-td)) .* (abs(t-td) < abs(t(2)-t(1)));

no_cycles = 3;
beta = 10 ^ (-60 / 20);
T = no_cycles / (2 * params.centre_freq * sqrt(log(1/beta)));
fn_toneburst = @(t, td) exp(-((t - td) / T) .^ 2) .* sin(2 * pi * params.centre_freq * (t - td));

forcing_function = fn_toneburst;
% forcing_function = fn_discrete_delta;

switch lower(inc_mode)
    case {'l', 'long', 'longitudinal'}
        tmp_dsps = forcing_function(params.time, bdry_dists / params.long_vel);
        bdry_dsps = [-tmp_dsps * cos(inc_angle); -tmp_dsps * sin(inc_angle); zeros(size(tmp_dsps))];
    case {'s', 'shear', 't', 'transverse'}
    otherwise
        error('Unknown mode')
end

dfs = 1:3;
bdry_dfs = kron(dfs(:), ones(size(bdry_nds)));
bdry_nds = kron(ones(numel(dfs), 1), bdry_nds);

%Get the matrices from the global modal
fe_options.return_matrices_only = 1;
mats = fn_FE_entry_point(mod, matls, el_types, [], fe_options);

gl_i = mats.gl_lookup(sub2ind(size(mats.gl_lookup), bdry_nds, bdry_dfs));

%be mindful of these lines
%[mn_res_i, gl_i, bdry_nds, bdry_dfs, bdry_lyrs] = fn_get_main_subdomain_mappings(main.doms{d}.mod, main.res, main.res.mats.gl_lookup);
%K_sub = main.res.mats.K(gl_i, gl_i);

%this is the crux line to get to
fe_options.solver_mode = 'exp';
[frcs, frce_set] = fn_convert_disps_to_forces_v2(...
            mats.K(gl_i, gl_i), mats.C(gl_i, gl_i), mats.M(gl_i, gl_i), ...
            params.time(2) - params.time(1), bdry_dsps, mod.bdry_lyrs(bdry_nds), 'in', fe_options.solver_mode);
steps{1}.load.frc_nds = bdry_nds(frce_set);
steps{1}.load.frc_dfs = bdry_dfs(frce_set);
steps{1}.load.time = params.time;
steps{1}.load.frcs = frcs;
steps{1}.mon.dsp_nds = bdry_nds;
steps{1}.mon.dsp_dfs = bdry_dfs;

fe_options.return_matrices_only = 0;
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

if ~isempty(res{1}.fld)
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    anim_options.fld_time = res{1}.fld_time;
    anim_options.pause_value = 0.05;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end





