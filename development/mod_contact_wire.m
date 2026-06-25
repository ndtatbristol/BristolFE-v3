function [mod, matls, el_types, steps, fe_options, params] = mod_contact_wire(params)
%This is a parametric model description of a 1D waveguide (default =
%circular cross-section)

default_params.els_per_wavelength = 40;
% default_params.abs_bdry_thickness_in_wavelengths = 1;

%Material properties
default_params.matl_name = 'copper';
default_params.section_fname = 'contact_wire107.mat';
default_params.array_fname = 'contact_wire_array.mat';

%Details of input signal applied at source
default_params.centre_freq = 50e3; %used to determine element size - final result is in frequency domain anyway
default_params.no_cycles = 6;
default_params.source_direction = 3;

%Shape and length of cross section - can be fully specified in 2D, or 2D perimeter
%can be specified or it can just be a rod of specified diameter
default_params.z_max =  500e-3;
default_params.z_min = -250e-3;

%Run for long enough for longitudinal waves to travel this many times
%the length of the waveguide
default_params.max_time_multiplier = 6; 

default_params.safety_factor = 3;

%Solver options - specify how ofter field output is produced to use in
%animation
default_params.fe_options.field_output_every_n_frames = inf;

%--------------------------------------------------------------------------
params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
params = fn_set_default_fields(params, default_params);
fe_options = params.fe_options;
el_types = [fn_2d_el_types(), fn_3d_el_types()];

el_typ_to_use_for_solid_2d = 'CPE3';
el_typ_to_use_for_fluid_2d = 'AC2D3';
el_typ_to_use_for_solid_3d = 'C3D6';
el_typ_to_use_for_fluid_3d = 'AC3D6';

%BUILD THE MODEL USING THE PARAMETERS GIVEN ABOVE

matl_i = 1; 
matls{matl_i} = fn_material_library(params.matl_name);

max_wavelength = fn_estimate_max_min_wavelengths(matls, params.centre_freq);
% abs_bdry_thickness = params.abs_bdry_thickness_in_wavelengths * max_wavelength;

%Define the material
%A cell array with an entry for each material used in the model is required.
%In this case there is only one material, but the index of this material in 
%the cell array is given a name, matl_i, so you can see where it is used
%later when elements in the model are assigned to elements
matl_i = 1; %material index is given a name so you can see where it appears later
matls{matl_i} = fn_material_library(params.matl_name);

%Work out element size and time step
params.el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, params.el_size, params.safety_factor);

tmp = load(params.section_fname);
mod2d.nds = tmp.mod.nds;
mod2d.els = tmp.mod.els;

params.array = load(params.array_fname)
params.array = params.array.array;

%Correct element size so that rows are at integer spacing (assumes rows are
%equi-spaced!
tmp = mean(abs(diff(params.array.row_pos)));
params.el_size = tmp / ceil(tmp / params.el_size);

%------DEBUG
params.array.row_pos = 0;
params.array.trans_pos = [1,2];
params.array.trans_row = [1,1];
params.array.trans_pos_orientations = [
    0,1,0
    0,0,1];
params.array.trans_node_list = [379 373];
%------

%find XY positions of transducers from original section mesh
trans_pos_xy = [mod2d.nds(params.array.trans_node_list,1),mod2d.nds(params.array.trans_node_list,2)];

mod2d = fn_remesh(mod2d, params.el_size);

mod2d.el_mat_i = zeros(size(mod2d.els, 1), 1);
mod2d.el_typ_i = ones(size(mod2d.els, 1), 1);
mod2d.el_abs_i = zeros(size(mod2d.els, 1), 1);

%Create the nodes and elements of the mesh
if isfield(params, 'z_pts')
    z_pts = params.z_pts;
else
    z_pts = linspace(params.z_min, params.z_max, ceil((params.z_max - params.z_min) / params.el_size));
end
mod = fn_extrude_2d_mesh(mod2d, z_pts);


mod.el_mat_i(:) = matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid_3d)); %extracts the index of the chosen element type from the cell array of possible element types

%Identify the transducer nodes
no_trans = numel(params.array.trans_pos);
trans_nds = zeros(no_trans, 1);
trans_z = params.array.row_pos(params.array.trans_row(:));
trans_pos_xyz = [trans_pos_xy(params.array.trans_pos(:), :), trans_z(:)];
for i = 1:no_trans
    trans_nds(i) = fn_find_node_nearest_to_point(mod.nds, trans_pos_xyz(i, :), params.el_size);
end

trans_nds = repmat(trans_nds, [1,3]);
trans_dfs = repmat([1:3], [no_trans, 1]);
trans_wts = params.array.trans_pos_orientations;

% figure;plot3(mod.nds(:,1), mod.nds(:,2), mod.nds(:,3), 'k.');
% hold on;plot3(mod.nds(trans_nds,1), mod.nds(trans_nds,2), mod.nds(trans_nds,3), 'ro');

%Provide the time signal for the loading
params.model_size = max(mod.nds(:,3)) - min(mod.nds(:,3));
[vel, ~] = fn_estimate_max_min_vels(matls{matl_i});
max_time = params.model_size / vel * params.max_time_multiplier;
in_time = 0: time_step: max_time;
in_signal = fn_gaussian_pulse(in_time, params.centre_freq, params.no_cycles);

%create the load steps
for s = 1:no_trans
    steps{s}.load.frc_nds = trans_nds(s, :);
    steps{s}.load.frc_dfs = trans_dfs(s, :);
    steps{s}.load.frc_wts = trans_wts(s, :);
    steps{s}.load.time = in_time;
    steps{s}.load.frcs = in_signal;

    %Where the displacement should be monitored (same for all steps)
    steps{s}.mon.dsp_nds = trans_nds(:);
    steps{s}.mon.dsp_dfs = trans_dfs(:);
    steps{s}.mon.dsp_wts = zeros(no_trans, numel(trans_nds));
    no_dfs = size(trans_wts, 2);
    for t = 1:no_trans
        cols = ([1:no_dfs]-1) * no_trans + t;
        steps{s}.mon.dsp_wts(t, cols) = trans_wts(t, :);
    end
end

%---------debug (reciprocity test)
% steps{1}.load.frc_nds = [137, 489];
% steps{1}.load.frc_dfs = [2, 3];
% steps{1}.load.frc_wts = randn(1,2);
% steps{1}.load.time = in_time;
% steps{1}.load.frcs = in_signal;
% steps{1}.mon.nds = [956, 1024];
% steps{1}.mon.dfs = [1, 3];
% steps{1}.mon.wts = randn(1,2);
% 
% steps{2}.load.frc_nds = steps{1}.mon.nds;
% steps{2}.load.frc_dfs = steps{1}.mon.dfs;
% steps{2}.load.wts = steps{1}.mon.wts;
% steps{2}.load.time = in_time;
% steps{2}.load.frcs = in_signal;
% steps{2}.mon.nds = steps{1}.load.frc_nds;
% steps{2}.mon.dfs = steps{1}.load.frc_dfs;
% steps{2}.mon.wts = steps{1}.load.wts;

%---------debug


end

function mod = fn_remesh(mod, el_size)
free_ed = fn_find_free_edges(mod.els);
%TODO - need to get perim nodes in order; make dummy axis s of distance around perimeter and use interp1 to interpolate x and
%y with desired number of new nodes in s. Then remesh using this perimeter.
ed_nds = fn_order_perimeter(free_ed);
ed_nds = [ed_nds(:); ed_nds(1)];
s = [0; cumsum(sqrt(sum((mod.nds(ed_nds(2:end), :) - mod.nds(ed_nds(1:end-1), :)) .^ 2, 2)))];
s2 = linspace(s(1), s(end), ceil((max(s)-min(s)) / el_size))';
ed_nds2 = [interp1(s, mod.nds(ed_nds, 1), s2), interp1(s, mod.nds(ed_nds, 2), s2)];
ed_nds2 = ed_nds2(1:end-1,:);

model = createpde;
gd = [2; size(ed_nds2, 1); ed_nds2(:, 1); ed_nds2(:, 2)];   % polygon geometry
ns = char('P1')';
sf = 'P1';
g = decsg(gd, sf, ns);

geometryFromEdges(model, g);
Hmax = el_size;
Hmin = el_size;
mesh = generateMesh(model, ...
    'GeometricOrder', 'linear', ...
    'Hmax', Hmax, ...
    'Hmin', Hmin);

mod.nds = mesh.Nodes';       % (n_nodes x 2)
mod.els = mesh.Elements';    % (n_elements x 3)

end

function nodes = fn_order_perimeter(edges)
% edges: n x 2 array of node indices
% nodes: n x 1 ordered node list around the loop

    % Build adjacency list
    maxNode = max(edges(:));
    adj = cell(maxNode,1);
    for i = 1:size(edges,1)
        a = edges(i,1);
        b = edges(i,2);
        adj{a}(end+1) = b;
        adj{b}(end+1) = a;
    end

    % Start from first edge
    nodes = zeros(size(edges,1),1);
    nodes(1) = edges(1,1);
    nodes(2) = edges(1,2);

    % Walk around the loop
    for k = 3:length(nodes)
        prev = nodes(k-2);
        curr = nodes(k-1);
        nbrs = adj{curr};
        
        % pick the neighbour that is not the previous node
        if nbrs(1) == prev
            nodes(k) = nbrs(2);
        else
            nodes(k) = nbrs(1);
        end
    end
end