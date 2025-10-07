function [mod, el_types] = fn_create_el_types_for_legacy_v2_models(mod, matls)
%el_type used to be part of material description with a cell array, 
%el_typ_i, of size n_els x 1 with string giving element type for each 
%element based on material. In v3, el_typ_i is an n_els x 1 matrix of
%indices into a single cell array el_types that contains strings for each
%element type used in model

%v2 - el_typ_i is cell array or does not exist


if isfield(mod, 'el_typ_i') 
    if iscell(mod.el_typ_i)
        %in this case it just needs converting to indices and a separate lookup
        %of possible el_types
        el_types = unique(mod.el_typ_i);
        tmp = zeros(size(mod.els, 1), 1);
        for i = 1:numel(el_types)
            tmp(strcmp(el_types{i}, a)) = i;
        end
        mod.el_typ_i = tmp;
        mod.el_types = el_types;%this is a bodge - see below
    else
        %in this case, el_typ_i is OK (probably this function has been
        %called already, so just el_types needs setting
        el_types = mod.el_types;
    end
else
    mod.el_typ_i = zeros(size(mod.els, 1), 1);
    un_mats = unique(mod.el_mat_i);
    for m = 1:numel(un_mats)
        el_types{un_mats(m)} = matls{un_mats(m)}.el_typ;
        mod.el_typ_i(mod.el_mat_i == un_mats(m)) = un_mats(m);
    end
    mod.el_types = el_types;
end
if ~isfield(mod, 'el_abs_i')
    mod.el_abs_i = zeros(size(mod.els, 1), 1);
end

end