%test of waveguide code
clear
% close all


model_to_run = @mod_waveguide;

params.fe_options = [];

show_geom_only = 0;
 
n_pts = 24; radius = 5e-3;

%create unstructured mesh
model = createpde;
theta = linspace(0, 2*pi, n_pts+1)';
theta(end) = [];  % remove duplicate point
p = radius * [cos(theta), sin(theta)];
gd = [2; size(p,1); p(:,1); p(:,2)];   % polygon geometry
ns = char('P1')';
sf = 'P1';
g = decsg(gd, sf, ns);



geometryFromEdges(model, g);
Hmax = 2 * pi * radius / n_pts;
Hmin = 2 * pi * radius / n_pts;
mesh = generateMesh(model, ...
    'GeometricOrder', 'linear', ...
    'Hmax', Hmax, ...
    'Hmin', Hmin);

params.nds_2d = mesh.Nodes';       % (n_nodes x 2)
params.els_2d = mesh.Elements';    % (n_elements x 3)

%--------------------------------------------------------------------------
%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
display_options.node_sets_to_plot(2).nd = steps{1}.mon.nds;
display_options.node_sets_to_plot(2).col = 'g.';
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
