function res = fn_query_el_type_info(varargin)
%Function for querying element type info table, which is stored in
%fn_el_type_info function

el_type_info = fn_el_type_info();

%No arguments - just return whole lot
if nargin == 0
    res = el_type_info;
    return
end

%One argument - return info for element with that name
if nargin == 1
    name = varargin{1};
    for i = 1:numel(el_type_info)
        if isequal(el_type_info{i}.name, name)
            res = el_type_info{i};
            return
        end
    end
    error('Element name not found')
end

if nargin == 2
    query = varargin{1};
    value = varargin{2};
    %Return cell array of element names that match el_type_info.(query) == value
    res = {};
    for i = 1:numel(el_type_info)
        if isequal(el_type_info{i}.(query), value)
            res{end+1} = el_type_info{i}.name;
        end
    end
else
    error('Expects 0, 1, or 2 inputs')
end
end
