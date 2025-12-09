clear

%Test of the various functions to remove nodes from a model and rejig
%element table accordingly
rng(1);

%Initial nodes and elements
nds1 = rand(10,2);
els1 = [1,4,7,5,3,3];

%Table with original node numbers in first column and associated nodes in
%2nd and 3rd columns. All subsequent versions of this should have the same
%data in the same order in 2nd and 3rd columns but different node numbers
%in first column
orig1 = [els1', nds1(els1,:)];

%Remove unused nodes function
[nds2, els2, old_nds, new_nds] = fn_remove_unused_nodes(nds1, els1);
orig2 = [els2', nds2(els2,:)]

isequal(orig1(:,2:end), orig1(:,2:end))

%Can also do the element remapping like this
els2b = fn_remap_matrix(els1, new_nds);
isequal(els2, els2b)

%Sorting nodes (e.g. to get into coordinate order for Pogo)
[nds3, i] = sortrows(nds2, 'descend');
i = fn_inverse_map(i);
els3 = fn_remap_matrix(els2, i);
orig3 = [els3', nds3(els3,:)];

isequal(orig1(:,2:end), orig3(:,2:end))
