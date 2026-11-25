function type = fn_matl_type(matl)
%Returns type of material (e.g. solid, fluid)
type = 'unknown';
if isfield(matl, 'D')
    if all(size(matl.D) == 1)
        type = 'fluid';
    elseif all(size(matl.D) == 6)
        type = 'solid';
    end
end
end
