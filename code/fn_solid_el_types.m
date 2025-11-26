function el_types = fn_solid_el_types()
%SUMMARY
%   Returns cell array of available solid element types. 
%INPUTS
%   None
%OUTPUTS
%   el_types - cell array of element names that model solids
el_types = fn_query_el_type_info('state', 'solid');
end