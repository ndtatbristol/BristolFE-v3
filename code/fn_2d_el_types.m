function el_types = fn_2d_el_types()
%SUMMARY
%   Returns cell array of available 2D element types. Use in 2D models,
%   usually exactly like this:
%   el_types = fn_2d_el_types()
%INPUTS
%   None
%OUTPUTS
%   el_types - cell array of element names that can be used in 2D models
el_types = fn_query_el_type_info('dims', 2);
end