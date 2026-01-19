function nds = fn_nds_on_els(els, els_to_consider)
%Returns list of unique nodes from els(els_to_consider, :)
nds = unique(els(els_to_consider, :));
end