clear

%Test of the various functions to remove nodes from a model and rejig
%element table accordingly
rng(1);

%Initial nodes and elements
nds1 = rand(10,2);
els1 = [1,4,7,5,3,3];
%example of something else that references nodes by number
nd_refs1 = [3,5,7];

%Table with original node numbers in first column and associated nodes in
%2nd and 3rd columns. All subsequent versions of this should have the same
%data in the same order in 2nd and 3rd columns but different node numbers
%in first column
el_test1 = [els1', nds1(els1,:)];
ref_test1 = [nd_refs1', nds1(nd_refs1,:)];

%Remove unused nodes function
[nds2, els2, old_nds, new_nds] = fn_remove_unused_nodes(nds1, els1);
el_test2 = [els2', nds2(els2,:)];
nd_refs2 = fn_remap_matrix(nd_refs1, new_nds);
ref_test2 = [nd_refs2', nds2(nd_refs2,:)];

fprintf('Elements after node removal: %i\n',isequal(el_test2(:,2:end), el_test1(:,2:end)));
fprintf('Node references after node removal: %i\n',isequal(ref_test2(:,2:end), ref_test1(:,2:end)));


%Can also do the element remapping like this
els2b = fn_remap_matrix(els1, new_nds);
fprintf('(Alternative elements remapping outside function: %i)\n',isequal(els2, els2b));

%Sorting nodes (e.g. to get into coordinate order for Pogo)
[nds3, i] = sortrows(nds2, 'descend');
i = fn_inverse_map(i);
els3 = fn_remap_matrix(els2, i);
el_test3 = [els3', nds3(els3,:)];
nd_refs3 = fn_remap_matrix(nd_refs2, i);
ref_test3 = [nd_refs3', nds3(nd_refs3,:)];
fprintf('Elements after node sort: %i\n',isequal(el_test3(:,2:end), el_test1(:,2:end)));
fprintf('Node references after node sort: %i\n',isequal(ref_test3(:,2:end), ref_test1(:,2:end)));

%Attempt to recreate what happens around validation model
fprintf('\n\n\nVALIDATION MODELLING\n\n')
mn_mod_nds = rand(10,2);
mn_mod_els = randperm(10);
orig_ref_to_mn_mod_nds = [7,6,9]'; %these need to still refer to same physical nodes at the end
true_refs = [orig_ref_to_mn_mod_nds, mn_mod_nds(orig_ref_to_mn_mod_nds, :)];

val_mod_nds = [mn_mod_nds; rand(3,2)]; %val model
val_mod_els = [[6:10], [13,11,12]]; %some main model nodes are no longer needed, but refs to some new nodes are added
true_els = [val_mod_els', val_mod_nds(val_mod_els, :)];

if 0
    %Drop the unused nodes and update elements
    [val_mod_nds1, val_mod_els1, old_nds1, new_nds1] = fn_remove_unused_nodes(val_mod_nds, val_mod_els);

    %Update references
    ref_to_val_nds1 = fn_remap_matrix(orig_ref_to_mn_mod_nds, new_nds1);

    %Checks at this point
    check_refs1 = [ref_to_val_nds1, val_mod_nds1(ref_to_val_nds1, :)];
    check_els1 = [val_mod_els1', val_mod_nds1(val_mod_els1, :)];
    fprintf('Node references to val model after node drop: %i\n',isequal(check_refs1(:,2:end), true_refs(:,2:end)));
    fprintf('Element node references to val model after node drop: %i\n',isequal(check_els1(:,2:end), true_els(:,2:end)));

    %And the mapping back to original node numbers
    orig_refs_recovered1 = fn_remap_matrix(ref_to_val_nds1, old_nds1);
    fprintf('Equivalent original node references recovered for val model after node drop: %i\n', isequal(orig_refs_recovered1, orig_ref_to_mn_mod_nds));

    %Now add the node sorting step
    [val_mod_nds2, tmp] = sortrows(val_mod_nds1, (size(val_mod_nds1,2):-1:1));
    new_nds2 = fn_inverse_map(tmp);
    val_mod_els2 = fn_remap_matrix(val_mod_els1, new_nds2);
    %Update elements - OK
    check_els2 = [val_mod_els2', val_mod_nds2(val_mod_els2, :)];

    %Update references - OK
    ref_to_val_nds2 = new_nds2(ref_to_val_nds1);
    check_refs2 = [ref_to_val_nds2, val_mod_nds2(ref_to_val_nds2, :)];

    fprintf('Node references to val model after node drop and sort: %i\n',isequal(check_refs2(:,2:end), true_refs(:,2:end)));
    fprintf('Element node references to val model after node drop and sorting: %i\n',isequal(check_els2(:,2:end), true_els(:,2:end)));

    %And the mapping back to original node numbers - OK
    old_nds2 = fn_inverse_map(new_nds2);
    orig_refs_recovered2 = fn_remap_matrix(fn_remap_matrix(ref_to_val_nds2, old_nds2), old_nds1); %undo both mappings!
    fprintf('Equivalent original node references recovered for val model after node drop: %i\n', isequal(orig_refs_recovered2, orig_ref_to_mn_mod_nds));
end
%KEY STEPS IN FUNCTION

%Drop the unused nodes and update elements
[val_mod_nds1, val_mod_els1, old_nds1, new_nds1] = fn_remove_unused_nodes(val_mod_nds, val_mod_els);

%Sort nodes
[val_mod_nds2, tmp] = sortrows(val_mod_nds1, (size(val_mod_nds1,2):-1:1));
new_nds2 = fn_inverse_map(tmp);
%Update elements
val_mod_els2 = fn_remap_matrix(val_mod_els1, new_nds2);
%Check elements
check_els2 = [val_mod_els2', val_mod_nds2(val_mod_els2, :)];

%Work out the old_nds and new_nds return values
j = new_nds1 == 0;
new_nds1(j) = 1;
new_nds = new_nds2(new_nds1);
old_nds = fn_inverse_map(new_nds);
new_nds(j) = 0;

fprintf('Element node references to val model after node drop and sorting: %i\n',isequal(check_els2(:,2:end), true_els(:,2:end)));

%END OF FUNCTION returns old_nds & new_nds

%IN MAIN CODE

%Update references
ref_to_val_nds = fn_remap_matrix(orig_ref_to_mn_mod_nds, new_nds);
check_refs = [ref_to_val_nds, val_mod_nds2(ref_to_val_nds, :)];
fprintf('Node references to val model after node drop and sort: %i\n',isequal(check_refs(:,2:end), true_refs(:,2:end)));

%Check original refs can be recovered
orig_refs_recovered = fn_remap_matrix(ref_to_val_nds, old_nds); %undo both mappings!
fprintf('Equivalent original node references recovered for val model after node drop and sort: %i\n', isequal(orig_refs_recovered, orig_ref_to_mn_mod_nds));
