function el_types = fn_fluid_el_types()
%USAGE
%   el_types = fn_fluid_el_types()
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns cell array of available fluid element types. 
%INPUTS
%   None
%OUTPUTS
%   el_types - cell array of element names that model fluids
%--------------------------------------------------------------------------

el_types = fn_query_el_type_info('state', 'fluid');
end