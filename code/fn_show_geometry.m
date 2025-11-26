function h_patch = fn_show_geometry(mod, matls, options)
%SUMMARY
%   Plots the model geometry and returns handle to patches used for each
%   element that can be used as an input argument for subsequent field 
%   animations if desired.

if isempty(matls)
    %if materials not specified make up some default colours
    no_unique_mats = numel(unique(mod.el_mat_i));
    default_options.matl_cols = linspace(0, 1, no_unique_mats)' * ones(1,3);
else
    %Deal with legacy code where materials is structure array rather than cell array
    if isstruct(matls)
        matls = arrayfun(@(x) x, matls, 'UniformOutput', false);
    end
    default_options.matl_cols = zeros(numel(matls), 3);
    for i = 1:numel(matls)
        if isempty(matls{i})
            continue
        end
        default_options.matl_cols(i, :) = matls{i}.col;
    end

end

if ~isfield(mod, 'el_abs_i')
     mod.el_abs_i = zeros(size(mod.els, 1), 1);
end

% addpath(genpath('..'));
default_options.offset = 0;
default_options.scale = 1;
default_options.norm_value = []; %empty to normalise to max

options = fn_set_default_fields(options, default_options);



% if isfield(mod, 'el_mat_i')
%     options.el_mat_i = mod.el_mat_i;
% else
%     options.el_mat_i = ones(size(mod.els, 1), 1);
% end


% if isfield(mod, 'el_abs_i')
%     options.el_abs_i = mod.el_abs_i;
% else
%     options.el_abs_i = zeros(size(mod.els, 1), 1);
% end

% h_patch = fn_display_result_v2(mod.nds * options.scale + options.offset, mod.els, options);
h_patch = fn_display_result_v3(mod, options);
drawnow;
end

