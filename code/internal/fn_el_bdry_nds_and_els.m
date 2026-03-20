function [nds_in, nds_on, nds_out, els_in, els_out] = fn_el_bdry_nds_and_els(els, els_in_bdry)
%Returns lists of nodes and elements at boundary between els(els_in_bdry,:) and els(~els_in_bdry,:)
%nds_on are the nodes on the boundary; els_in are elements
%adjacent to and within boundary; nds_in are the other nodes on these
%elements that are not in nds_on. els_out and nds_out are equivalent 
%outside boundary 
els_not_in_bdry =  ~els_in_bdry;
nds_on = intersect(fn_nds_on_els(els, els_in_bdry), fn_nds_on_els(els, els_not_in_bdry));
els_in = fn_els_on_nds(els, nds_on) & els_in_bdry;
els_out = fn_els_on_nds(els, nds_on) & els_not_in_bdry;
nds_in = setdiff(fn_nds_on_els(els, els_in), nds_on);
nds_out = setdiff(fn_nds_on_els(els, els_out), nds_on);
end
