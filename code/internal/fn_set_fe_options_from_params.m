function fe_options = fn_set_fe_options_from_params(params)
%Utility function for extracting fieldfnames prefixed with fe_options_* to
%separate structure. Used in examples.
tmp = fieldnames(params);
fe_options_str = 'fe_options_';
fe_options = [];
for i = 1:numel(tmp)
    if startsWith(tmp{i}, fe_options_str)
        fe_options.(tmp{i}(numel(fe_options_str)+1:end)) = params.(tmp{i});
    end
end
end