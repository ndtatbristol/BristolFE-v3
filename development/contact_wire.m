%test of waveguide code
clear
% close all


model_to_run = @mod_contact_wire;

params.els_per_wavelength = 40;
params.fe_options = [];
params.z_max = 1000e-3;

show_geom_only = 0;

%--------------------------------------------------------------------------
%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'g.';
display_options.draw_elements = 1;
if show_geom_only
    figure;
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

figure;
ascan = sum(res{1}.dsps, 1);
plot(steps{1}.load.time, ascan);
hold on;
xlabel('Time (s)')
title(sprintf('%i EPW', params.els_per_wavelength))