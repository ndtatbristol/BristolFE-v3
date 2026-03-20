function [h_patches, scale, offsets] = fn_show_geometry_with_subdomains(main, display_options)

default_options.sep_frac = 0.1;
default_options.subdomain_cols = 'rgbcmy';
default_options.subdomain_plot_positions = 'south';
default_options.subdomain_plot_horizontal_alignment = 'left';
default_options.subdomain_plot_vertical_alignment = 'top';
display_options = fn_set_default_fields(display_options, default_options);

%Get dimensions of main model and max dimensions of any subdomain
no_doms = numel(main.doms);
max_dom_size = 0;
for d = 1:no_doms
    max_dom_size = max(max_dom_size, max(main.doms{d}.mod.nds) - min(main.doms{d}.mod.nds));
end
main_min = min(main.mod.nds);
main_max = max(main.mod.nds);
main_dims = main_max - main_min;
sep = sqrt(sum((main_dims .^ 2))) * display_options.sep_frac;

%General form of coordinate transform for sudomain d (d = 1...no_doms) is
%(xd, yd) = offset + step * d + ((x, y) - min(nds)) * scale
%         = (offset + step * d - min(nds) * scale) + (x, y) * scale
%         = offsets(d, :) + (x, y) * scale

%Work out scale factor for sub-domain plots
switch lower(display_options.subdomain_plot_positions)
    case {'south', 'north'}
        total_dom_plot_size = no_doms * max_dom_size(1);
        scale = min((main_dims(1) - (no_doms - 1) * sep) / total_dom_plot_size, 1);
    case {'east', 'west'}
        total_dom_plot_size = no_doms * max_dom_size(2);
        scale = min((main_dims(2) - (no_doms - 1) * sep) / total_dom_plot_size, 1);
end

%work out step between subdomain plots 
switch lower(display_options.subdomain_plot_positions)
    case {'south', 'north'}
        step = [max_dom_size(1), 0] * scale + [1, 0] * sep;
    case {'east', 'west'}
        step = [0, -max_dom_size(2)] * scale + [0, -1] * sep; %-ve because these go down from top (i.e. decreasing y)
end

%work out offset for first subdomain plot (effectively position for s = 0)
switch lower(display_options.subdomain_plot_positions)
    case 'south'
        offset = [main_min(1), main_min(1)] - [0, sep] - [0, max_dom_size(2)] * scale - step;
    case 'north'
        offset = [main_min(1), main_max(2)] + [0, sep] - step;
    case 'west'
        offset = [main_min(1), main_max(2)] - [sep, 0] - [max_dom_size(1), max_dom_size(2)] * scale - step;
    case 'east'
        offset = [main_max(1), main_max(2)] + [sep, 0] - [0, max_dom_size(2)] * scale - step;
end

%Work out offsets for individual plots
offsets = zeros(no_doms, 2);
for d = 1:no_doms
    subdomain_origin = [0, 0];
    switch lower(display_options.subdomain_plot_horizontal_alignment)
        case 'left'
            subdomain_origin(1) = min(main.doms{d}.mod.nds(:, 1));
        case 'right'
            subdomain_origin(1) = max(main.doms{d}.mod.nds(:, 1)) - max_dom_size(1);
    end
    switch lower(display_options.subdomain_plot_vertical_alignment)
        case 'top'
            subdomain_origin(2) = max(main.doms{d}.mod.nds(:, 2)) - max_dom_size(2);
        case 'bottom'
            subdomain_origin(2) = min(main.doms{d}.mod.nds(:, 2));
    end
    offsets(d, :) = offset + step * d - subdomain_origin * scale;
end

%Main
p = 1;
display_options.scale = 1;
display_options.offset = 0;
h_patches{p} = fn_show_geometry(main.mod, main.matls, main.el_types, display_options);
display_options.node_sets_to_plot = [];
hold on;

%Subdomains
display_options.scale = scale;
for d = 1:no_doms
    p = p + 1;
    display_options.offset = offsets(d, :);
    h_patches{p} = fn_show_geometry(main.doms{d}.mod, main.matls, main.el_types, display_options);
    col = display_options.subdomain_cols(rem(d - 1, numel(display_options.subdomain_cols)) + 1);
    tmp = main.doms{d}.mod.inner_bndry_pts;
    tmp = [tmp; tmp(1,:)];
    xy = fn_dom_coord(tmp, display_options.scale, display_options.offset);
    plot(xy(:, 1), xy(:, 2), col);
    xy = fn_main_coord(tmp + min(main.mod.nds), -min(main.mod.nds));
    plot(xy(:, 1), xy(:, 2), col);
end

end

function xy = fn_dom_coord(p, scale, offset)
xy = p * scale + offset;
end

function xy = fn_main_coord(p, offset)
xy = p + offset;
end 
