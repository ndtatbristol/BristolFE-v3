function h = fn_label_subplots(options)
%SUMMARY
%   Labels subplots (e.g. (a), (b) , (c) ...) on current figure
%AUTHOR
%   Paul Wilcox, 2007
%USAGE
%   fn_label_subplots(options)
%INPUTS
%   options - structured variable containing fields
%       options.left_right_first - [1] whether images are numbered
%       horizontally then vertically [default] or vice versa
%       options.font - ['Times'] font
%       options.font_size - [10] font size in points
%       options.first_label - ['a'] single character for first label (e.g. 'a',
%       'A', '1'), rest will be obtained by incrementing ASCII code
%       options.prepend -['('] string to put before label
%       options.append - [')'] string to put after label
%       options.xoffset - [NaN] shift x positions
%       options.yoffset - [NaN] shift y positions

default_options.left_right_first = 1;
default_options.font = 'Helvetica';
default_options.font_size = 10;
default_options.first_label = 'a';
default_options.prepend = '(';
default_options.append = ')';
default_options.xoffset = NaN;
default_options.yoffset = NaN;
options = fn_set_default_fields(options, default_options);

debug = 0;
%find positions [left, bottom, width, height] of all subplots
ch = findobj(gcf, 'Type', 'axes', '-not', 'Tag', 'legend', '-not', 'Tag', 'Colorbar');
if length(ch) < 2
    warning('No subplots found');
end
pos = [];
min_gap = inf;
for ii = 1:length(ch)
    set(ch(ii), 'Units', 'normalized');
    pos = [pos; get(ch(ii), 'Position')];
    if ii > 1
        gaps = pos(ii, 1) - [0; pos(1:ii - 1, 1) + pos(1:ii - 1, 3)];
        min_gap = min([min_gap, min(gaps(gaps > 0))]);
    end
end



%debugging
if debug
    for ii = 1:size(pos, 1)
        h = annotation('textbox', [pos(ii, 1), pos(ii, 2), 0, 0]);
        set(h, 'LineStyle', 'none');
        set(h, 'String', 'POS')
        set(h, 'VerticalAlignment', 'middle');
        set(h, 'HorizontalAlignment', 'center');
        set(h, 'FontName', options.font);
        set(h, 'FontSize', options.font_size);
    end
end

%calc coordinates of label centres
lab_pos = zeros(size(pos, 1), 2);
for ii = 1:size(pos, 1)
    lab_pos(ii, 1) = pos(ii, 1);
    lab_pos(ii, 2) = pos(ii, 2) + pos(ii, 4);
end

lab_pos(:, 1) = lab_pos(:, 1) - min_gap / 2;

if debug
    for ii = 1:size(pos, 1)
        h = annotation('textbox', [lab_pos(ii, 1), lab_pos(ii, 2), 0, 0]);
        set(h, 'LineStyle', 'none');
        set(h, 'String', 'LABPOS')
        set(h, 'VerticalAlignment', 'middle');
        set(h, 'HorizontalAlignment', 'center');
        set(h, 'FontName', options.font);
        set(h, 'FontSize', options.font_size);
    end
end

%sort into correct order for labelling
if options.left_right_first
    sort_order = [2, 1];
else
    sort_order = [1, 2];
end
lab_pos(:, 2) = -lab_pos(:, 2);
lab_pos = sortrows(lab_pos, sort_order);
lab_pos(:, 2) = -lab_pos(:, 2);

if ~isnan(options.xoffset)
    lab_pos(:,1) = lab_pos(:,1) + options.xoffset;
end;
if ~isnan(options.yoffset)
    lab_pos(:,2) = lab_pos(:,2) + options.yoffset;
end
lab_pos(lab_pos < 0) = 0;
lab_pos(lab_pos > 1) = 1;

%fill in labels
for ii = 1: size(lab_pos, 1)
    str = [options.prepend, char(double(options.first_label) + ii - 1), options.append];
    h(ii) = annotation('textbox',[lab_pos(ii, 1), lab_pos(ii, 2), 0, 0]);
    set(h(ii), 'LineStyle', 'none');
    set(h(ii), 'String', str)
    set(h(ii), 'VerticalAlignment', 'middle');
    set(h(ii), 'HorizontalAlignment', 'center');
    set(h(ii), 'FontName', options.font);
    set(h(ii), 'FontSize', options.font_size);
end

end