function all_els_in = fn_all_els_inside_bdry(els, nds_in, els_in)
%Returns logic vector of all els inside boundary defined by els(els_in,:)
%and nds_in which are the first set of nodes within boundary

last_nds = nds_in;
next_els = 1;
all_els_in = els_in;
while any(next_els)
    next_els = ~all_els_in & fn_els_on_nds(els, last_nds);
    all_els_in = all_els_in | next_els;
    last_nds = setdiff(fn_nds_on_els(els, next_els), last_nds); 
end

end