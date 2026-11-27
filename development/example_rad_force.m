clear 
close all
% restoredefaultpath

%Frequency domain solver - particle in fluid, calculation of radiation
%force

addpath(genpath(['..', filesep, 'code']));

addpath('models');

model_to_run = @mod_2d_freq_domain;
params = [];
params.include_scatterer = 1;
params.scatterer_shape = 'void';
params.els_per_wavelength = 10;
params.scatterer_size = 1e-3; %this is radius!
params.no_bdry_lyrs = 5;%note different to S-matrix style model

show_geom_only = 0;
show_field_output = 1;
animate_field_output = 0;
no_inc_angles = 90;
input_pressure = 1; %this is pressure of incident plane wave
%--------------------------------------------------------------------------
%Define the model
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

%incident wave vectors
inc_angles_degs = linspace(0, 360, no_inc_angles + 1)';
inc_angles_degs = inc_angles_degs(1:end - 1);
omega = params.freq * 2 * pi;
k = [cosd(inc_angles_degs), sind(inc_angles_degs)] * omega / params.fluid_vel;

%Convert these to displacements at the boundary nodes
nds_reduced = mod.nds(mod.bdry_lyrs ~= 0, :);
nd_i  = (1:size(mod.nds, 1))';
nd_i = nd_i(mod.bdry_lyrs ~= 0);
df = ones(size(nd_i)) * 4;
lyrs = mod.bdry_lyrs(mod.bdry_lyrs ~= 0);
lyrs_inc_dsp = lyrs == 3 | lyrs == 4;
lyrs_inc_frc = lyrs == 4 | lyrs == 5;
lyrs_mon_dsp = lyrs == 1 | lyrs == 2 | lyrs == 3;

%
P_plane_wave = input_pressure * exp(1i * nds_reduced * k') / (1i * omega); 

%Use special option to apply equiv forces for excitation directly in solver
steps{1}.load.freq = params.freq;
%Displacements specified on layers 1 and 2
steps{1}.load.dsp_nds = nd_i(lyrs_inc_dsp);
steps{1}.load.dsp_dfs = df(lyrs_inc_dsp);
steps{1}.load.dsps = P_plane_wave(lyrs_inc_dsp, :);
%Will be converted in solver to forces applied on layers 2 and 3
steps{1}.load.frc_nds = nd_i(lyrs_inc_frc);
steps{1}.load.frc_dfs = df(lyrs_inc_frc);
fe_options.input_dsps_as_equiv_frcs = 1;

%Record scattered field on bndr nodes
steps{1}.mon.dsp_nds = nd_i(lyrs_mon_dsp);
steps{1}.mon.dsp_dfs = df(lyrs_mon_dsp);

%Force integration contour facets
[facet_nds, facet_cent, facet_area, facet_norm] = fn_facets_between_el_sets(mod.nds, mod.els,  mod.el_typ_i, mod.el_set_i, el_types, 2, 3);
contour_nds = unique(facet_nds(:));

% contour_nds = nd_i(lyrs == 2);
% [facet_nds, facet_cent, facet_area, facet_norm] = fn_facets_from_nodes(mod.nds, mod.els, mod.el_typ_i, el_types, contour_nds);

if show_field_output
    fe_options.freq_domain_fld = 1;
else
    fe_options.freq_domain_fld = 0;
end
fe_options.field_output_type = 'sum(u)';

display_options = [];
bl_col = 'rb';
un_bl = unique(mod.bdry_lyrs);
un_bl = un_bl(un_bl > 0);
if show_geom_only
    figure;
    % mod.el_abs_i(mod.el_abs_i > 0) = 1;
    % mod.el_abs_i = mod.el_used / max(mod.el_used(:));
    for i = 1:numel(un_bl)
        display_options.node_sets_to_plot(i).nd = find(mod.bdry_lyrs == un_bl(i));
        display_options.node_sets_to_plot(i).col = [bl_col(rem(i,2) + 1), '.'];
    end
    i = i+1;
    display_options.node_sets_to_plot(i).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(i).col = 'ro';
    i = i+1;
    display_options.node_sets_to_plot(i).nd = steps{1}.mon.dsp_nds;
    display_options.node_sets_to_plot(i).col = 'go';
    i = i+1;
    display_options.node_sets_to_plot(i).nd = contour_nds
    display_options.node_sets_to_plot(i).col = 'bx';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    hold on;
    plot(facet_cent(:,1), facet_cent(:,2), 'b.');
    plot([mod.nds(facet_nds(:,1), 1), mod.nds(facet_nds(:,2), 1)]', ...
         [mod.nds(facet_nds(:,1), 2), mod.nds(facet_nds(:,2), 2)]', 'b');
    plot([facet_cent(:,1), facet_cent(:,1) + facet_norm(:, 1) * params.el_size]', ...
         [facet_cent(:,2), facet_cent(:,2) + facet_norm(:, 2) * params.el_size]', 'm');
    return
end

[res, mats] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%Compute radiation force at all incident angles

fn_console_output('Radiation force calculation ...');
t1 = clock;

force = zeros(no_inc_angles, 2);
G = fn_grad(mod.nds, mod.els, contour_nds); %gradient calc matrices
gradP = zeros(numel(contour_nds),numel(G));
P = zeros(size(mod.nds, 1), 1);
V = zeros(size(mod.nds, 1), 2);
for i = 1:no_inc_angles
    %Total pressure field at monitoring nodes
    P(steps{1}.mon.dsp_nds) = res{1}.dsps(:, i) * (1i * omega);

    %Total pressure gradient on contour
    for k = 1:numel(G)
        gradP(:,k) = G{k} * P;
    end
    
    %Total velocity field on contour
    V(contour_nds, :) = -1 / (1i * omega * params.fluid_density) * gradP;

    %debugging figure
    % if i == 1
    %     figure
    %     theta = atan2d(mod.nds(:,1), mod.nds(:,2));
    %     subplot(1,2,1);
    %     plot(theta(contour_nds), abs(P(contour_nds)), 'r.')
    %     subplot(1,2,2);
    %     plot(theta(contour_nds), abs(V(contour_nds)), 'b.')
    % end

    %Qtys at facet centres
    Vf = (V(facet_nds(:, 1), :) + V(facet_nds(:, 2), :)) / 2;
    Pf = (P(facet_nds(:, 1), :) + P(facet_nds(:, 2), :)) / 2;

    %Integrate over boundary
    % force(i, :) = sum((...
    %     params.fluid_density * sum(Vf .* conj(Vf), 2) / 2 + ...
    %     abs(Pf) .^ 2 /  (4 * params.fluid_density * params.fluid_vel ^ 2) ...
    %     ) .* facet_norm .* facet_area);
    force(i, :) = -sum(...
        (...
            (...
                (...
                + abs(Pf) .^ 2         / (4 * params.fluid_density * params.fluid_vel ^ 2) ... %pressure ^ 2 - OK
                - sum(abs(Vf) .^ 2, 2) * (params.fluid_density / 4) ... %vel ^ 2
                ) .* facet_norm ...
            ) ...
        + real(sum(Vf .* facet_norm, 2) .* conj(Vf)) * (params.fluid_density / 2) ... %momentum flux
        ) ...
    .* facet_area);
end
fn_console_output(sprintf(' completed in %.2f secs\n', etime(clock, t1)), [], 0);

figure;
plot(inc_angles_degs, force);
xlabel('Incident angle (^o)');
ylabel('Radiation force (N/m)');
legend('X', 'Y');

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

