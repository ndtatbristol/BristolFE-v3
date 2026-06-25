%test of waveguide code
clear
close all


model_to_run = @mod_contact_wire;

params.els_per_wavelength = 20;
params.fe_options = [];
params.z_max = 750e-3;
params.z_min = -250e-3;
% params.z_max = 75e-3;
% params.z_min = -25e-3;


show_geom_only = 0;

%--------------------------------------------------------------------------
%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

display_options.node_sets_to_plot(1).nd = steps{1}.mon.dsp_nds;
display_options.node_sets_to_plot(1).col = 'bo';
display_options.draw_elements = 1;
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    for s = 1:numel(steps)
        xyz = mod.nds(steps{s}.load.frc_nds(1),:);
        xyz = [xyz; xyz + steps{s}.load.frc_wts * params.el_size * 5];

        plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'r.-')
    end
    return
end

%Run model
res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%Parse to fmc_data
fmc_data = zeros(numel(steps{1}.load.time), numel(steps), numel(steps));
for s = 1:numel(steps)
    fmc_data(:, s, :) = permute(res{s}.dsps, [2, 3, 1]);
end

figure;
plot(steps{1}.load.time, fmc_data(:, 1, 2))
hold on;
plot(steps{1}.load.time, fmc_data(:, 2, 1), 'r--')
hold on;
xlabel('Time (s)')
title(sprintf('%i EPW', params.els_per_wavelength))