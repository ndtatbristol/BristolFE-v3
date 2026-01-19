function els = fn_els_on_nds(els, nds_to_consider)
%Returns logic vector of which els have one or mode node in nds_to_consder
els = any(ismember(els, nds_to_consider), 2);
end
