function mod = fn_add_fluid_solid_interface_els(mod, el_types)
%USAGE
%   mod = fn_add_fluid_solid_interface_els(mod, el_types)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Adds the necessary interface elements between all solid and fluid
%   elements in a model. Without these there is no coupling between the
%   solid and fluid domains. Only needs to be called once for a given model
%   after all features are added (i.e. it should be the last step)
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els.
%OUTPUT
%   mod - modified model with additional interface elements


%Deal with legacy v2 calls where args are mod, matls, varargin - THIS WILL
%NO LONGER WORK 26/11/25!!!!
if isstruct(el_types) && isfield(el_types, 'rho')
    matls = el_types;
    if isstruct(matls)
        matls = arrayfun(@(x) x, matls, 'UniformOutput', false);
    end
    [mod, el_types] = fn_create_el_types_for_legacy_v2_models(mod, matls);
end

solid_el_i = fn_el_types_of_state(el_types, 'solid');
fluid_el_i = fn_el_types_of_state(el_types, 'fluid');
int_el_typ_i = fn_el_types_of_state(el_types, 'fluid_solid_interface');

if ~(any(mod.el_typ_i == solid_el_i, 'all') && any(mod.el_typ_i == fluid_el_i, 'all'))
    %model has no solid or no fluid element types 
    return
end

%Remove existing interface elements - necessary otherwise new elements will
%be added on top and will effectively double the coupling between solid and
%fluid
mod = fn_remove_fluid_solid_interface_els(mod, el_types);

%for legacy v2 calls, need to embed el_types in mod as well
% mod.el_types = el_types;

%New method using find interface function (should work for 2D and 3D up
%to and including this function)
[interface_facets, interface_el_solid, interface_el_fluid] = fn_find_interface(mod, ismember(mod.el_typ_i, solid_el_i), ismember(mod.el_typ_i, fluid_el_i), el_types);

if isempty(interface_facets)
    %No interface found, so do nothing
    return
end

%Nodes in interface_facets need ordering so solid and fluid are 
%on correct sides for all elements - this is why mod.nds data is necessary
%for this function to work. Currently this part is specific to 2D.
no_int_els = size(interface_facets,1);
%Loop through each interface edge in turn and flip node order if necessary
%to they are all same way around
for i = 1:no_int_els
    %work out centre of fluid element adjoining this edge
    % e = els_adjoining_fluid_solid_interface_edges(i, fluid_or_solid(i, :) == 2);
    e = interface_el_fluid(i);
    ec = fn_calc_element_centres(mod.nds, mod.els(e,:));
    if size(mod.nds, 2) == 2
        %line between nodes
        a = mod.nds(interface_facets(i, 2), :) - mod.nds(interface_facets(i, 1), :);
        %line at right angle to line between nodes
        b = [a(2), -a(1)];
    else
        b = cross(...
            mod.nds(interface_facets(i, 2), :) - mod.nds(interface_facets(i, 1), :), ...
            mod.nds(interface_facets(i, 4), :) - mod.nds(interface_facets(i, 1), :));
    end
    %line from first node to ec
    c = ec - mod.nds(interface_facets(i, 1), :);
    %check sign of dot product
    if dot(c, b) < 0
        interface_facets(i, :) = fliplr(interface_facets(i, : ));
    end

end

%Add the new interface elements to the model

%First add their nodal indices to the element indices for the whole model
mod.els = [mod.els; [interface_facets, zeros(no_int_els, size(mod.els, 2) - size(interface_facets, 2))]];

%Get number of nodes per interface face for the different types of
%interface element available as these will need to be selected according to
%number of nodes on faces of elements on interface (only issue in 3d where
%they could be quad- or tri-shaped faces
nds_per_interface_face = zeros(1, numel(int_el_typ_i));
for i = 1:numel(int_el_typ_i)
    tmp = fn_query_el_type_info(el_types{int_el_typ_i(i)});
    nds_per_interface_face(i) = size(tmp.faces, 4 - tmp.dims); %tortuous way to get nodes per face that works in 2D and 3D!
end

%Set interface element types to those with appropriate number of nodes
tmp_int_el_typ_i = zeros(no_int_els, 1);
interface_nds_per_facet = sum(interface_facets > 0, 2);
for i = 1:numel(nds_per_interface_face)
    tmp_int_el_typ_i(interface_nds_per_facet == nds_per_interface_face(i)) = int_el_typ_i(i);
end

%Check none of the interface elements have not been assigned a type
if any(tmp_int_el_typ_i == 0)
    error('Some interface elements could not be defined')
end

%Append the interface element types to the model and zero the associated
%material and absorbing indices (which should be zero for interface
%elements)
mod.el_typ_i = [mod.el_typ_i; tmp_int_el_typ_i];
mod.el_mat_i = [mod.el_mat_i; zeros(no_int_els, 1)];
mod.el_abs_i = [mod.el_abs_i; zeros(no_int_els, 1)];

% mod.el_typ_i = [mod.el_typ_i; repmat(int_el_typ_i, [no_int_els, 1])];

end