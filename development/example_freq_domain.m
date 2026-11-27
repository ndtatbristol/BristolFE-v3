clear 
% close all
% restoredefaultpath

%Frequency domain solver

addpath(genpath(['..', filesep, 'code']));

addpath('models');

model_to_run = @mod_2d_freq_domain;
params = [];
params.include_scatterer = 1;
params.scatterer_shape = 'void';
params.els_per_wavelength = 20;
params.scatterer_size = 1e-3;

show_geom_only = 1;
show_field_output = 1;
animate_field_output = 0;
no_inc_angles = 180;
%--------------------------------------------------------------------------
%Define the model
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);
display_options = [];
col = 'rgbkm';
if show_geom_only
    figure;
    for i = 1:5
        display_options.node_sets_to_plot(i).nd = find(mod.bdry_lyrs == i);
        display_options.node_sets_to_plot(i).col = [col(i), '.'];
        h_patch = fn_show_geometry(mod, matls, el_types, display_options);
        hold on;
    end
    return
end

%incident wave vectors
inc_angles_degs = linspace(0, 360, no_inc_angles + 1)';
inc_angles_degs = inc_angles_degs(1:end - 1);
omega = params.freq * 2 * pi;
k = [cosd(inc_angles_degs), sind(inc_angles_degs)] * params.freq * 2 * pi / params.fluid_vel;

%Convert these to displacements at the boundary nodes
nds_reduced = mod.nds(mod.bdry_lyrs ~= 0, :);
nd_i  = (1:size(mod.nds, 1))';
nd_i = nd_i(mod.bdry_lyrs ~= 0);
df = ones(size(nd_i)) * 4;
lyrs = mod.bdry_lyrs(mod.bdry_lyrs ~= 0);
lyrs12 = lyrs == 1  | lyrs == 2;
lyrs23 = lyrs == 2  | lyrs == 3;
lyrs34 = lyrs == 3  | lyrs == 4;
u_plane_waves = exp(1i * nds_reduced * k'); %dot product handled in matrix mult

%Use special option to apply equiv forces for excitation directly in solver
steps{1}.load.freq = params.freq;
%Displacements specified on layers 1 and 2
steps{1}.load.dsp_nds = nd_i(lyrs12);
steps{1}.load.dsp_dfs = df(lyrs12);
steps{1}.load.dsps = u_plane_waves(lyrs12, :);
%Will be converted in solver to forces applied on layers 2 and 3
steps{1}.load.frc_nds = nd_i(lyrs23);
steps{1}.load.frc_dfs = df(lyrs23);
fe_options.input_dsps_as_equiv_frcs = 1;

%Use special option to get equiv forces for reception directly from solver
%Displacements monitored on layers 3 and 4
steps{1}.mon.dsp_nds = nd_i(lyrs34);
steps{1}.mon.dsp_dfs = df(lyrs34);
%Will be converted in solver to equivalent forces applied on layers 2 and 3
steps{1}.mon.frc_nds = nd_i(lyrs23);
steps{1}.mon.frc_dfs = df(lyrs23);
fe_options.output_dsps_as_equiv_frcs = 1;

fe_options.freq_domain_fld = 1;
fe_options.field_output_type = 'sum(u)';
[res, mats] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%Output is forces to generate scattered field on layers 2 and 3, so use
%reciprocity to convert back to plane waves
S = u_plane_waves(lyrs23, :).' * res{1}.frcs;

S_shifted = zeros(size(S));
for j = 1:size(S,2)
    S_shifted(:,j) = circshift(S(:,j), 1-j);
end

%Plot S-matrix and comparison with theory if relevant
figure;subplot(1,2,1);imagesc(inc_angles_degs, inc_angles_degs, abs(S)); axis xy; axis equal; axis tight; colorbar; %clim([0, 1e-4]);
subplot(1,2,2); 
cla;
h1 = plot(inc_angles_degs, abs(S_shifted) * pi, 'Color', [1,1,1] * 0.3);
if strcmp(params.scatterer_shape, 'void')
    hold on;
    S_theory = fn_s_matrix_circ_void(inc_angles_degs + 180, params.freq * 2 * pi / params.fluid_vel, params.scatterer_size);
    h2 = plot(inc_angles_degs, abs(S_theory), 'r');
    legend([h1(1), h2], 'FE', 'Theory');
end

if show_field_output
    if animate_field_output
        omega_t = exp(-1i * linspace(0, 2*pi, 12));
        fld = real(kron(res{1}.fld, omega_t));
        anim_options.repeat_n_times = 1;
        anim_options.fld_time = [];
    else
        fld = real(res{1}.fld);
        anim_options.repeat_n_times = 1;
        anim_options.fld_time = [];
    end
    figure;
    anim_options.clear_at_end = 0;
    anim_options.db_or_linear = 'linear';
    h_patch = fn_show_geometry(mod, matls, el_types, anim_options);
    fn_run_animation(h_patch, fld, anim_options);
end

