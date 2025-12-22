function el_types = fn_fluid_el_types()
%SUMMARY
%   Returns cell array of available fluid element types. 
%INPUTS
%   None
%OUTPUTS
%   el_types - cell array of element names that model fluids
%--------------------------------------------------------------------------

el_types = fn_query_el_type_info('state', 'fluid');
end