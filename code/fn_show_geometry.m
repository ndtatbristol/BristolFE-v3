function h_patch = fn_show_geometry(mod, matls, el_types, options)
%SUMMARY
%   Plots the model geometry and returns handle to patches used for each
%   element that can be used as an input argument for subsequent field 
%   animations if desired.

default_options.draw_elements = 0;
default_options.element_edge_color = [1,1,1] * 0.95;
default_options.element_face_colour = [1,1,1] * 0.75;
default_options.mesh_edge_color = 'k';
default_options.draw_mesh_edges = 1;
default_options.node_sets_to_plot = [];
default_options.scale_factor = [];
default_options.show_abs = 1;
default_options.matl_cols = [];
default_options.interface_el_col = [0,0,1];
default_options.transparency = 0.5;
default_options.scale = 1;
default_options.offset = 0;

%--------------------------------------------------------------------------

if isempty(matls)
    %if materials not specified make up some default colours
    %no_unique_mats = numel(unique(mod.el_mat_i));
    default_options.matl_cols = linspace(0, 1, max(mod.el_mat_i))' * ones(1,3);
else
    default_options.matl_cols = zeros(numel(matls), 3);
    for i = 1:numel(matls)
        default_options.matl_cols(i, :) = matls{i}.col;
    end
end

% addpath(genpath('..'));
default_options.offset = 0;
default_options.scale = 1;
default_options.norm_value = []; %empty to normalise to max

options = fn_set_default_fields(options, default_options);

h_patch = fn_display_result_v3(mod, el_types, options);
drawnow;
end

