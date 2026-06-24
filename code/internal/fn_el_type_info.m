function el_type_info = fn_el_type_info()

%2D modelling - solids - plane strain
i = 1;
el_type_info{i}.name = 'CPE3';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'triangular';

i = i + 1;
el_type_info{i}.name = 'CPE4';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'quadrilateral';

i = i + 1;
el_type_info{i}.name = 'CPE4R'; 
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'quadrilateral';

%2D modelling - solids - plane stress
i = i + 1;
el_type_info{i}.name = 'CPS3'; 
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'triangular';

i = i + 1;
el_type_info{i}.name = 'CPS4';  
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'quadrilateral';

i = i + 1;
el_type_info{i}.name = 'CPS4R';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'quadrilateral';

%test elements
i = i + 1;
el_type_info{i}.name = 'CPE3_new';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'triangular';

%2D modelling - fluids
i = i + 1;
el_type_info{i}.name = 'AC2D3';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'triangular';

i = i + 1;
el_type_info{i}.name = 'AC2D4';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'quadrilateral';

i = i + 1;
el_type_info{i}.name = 'AC2D4R';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'quadrilateral';

%test elements
% i = i + 1;
% el_type_info{i}.name = 'CPE4_new';
% el_type_info{i}.dims = 2;
% el_type_info{i}.state = 'solid';
% el_type_info{i}.shape = 'quadrilateral';
% 
% i = i + 1;
% el_type_info{i}.name = 'AC2D3_new';
% el_type_info{i}.dims = 2;
% el_type_info{i}.state = 'fluid';
% el_type_info{i}.shape = 'triangular';
% 
% i = i + 1;
% el_type_info{i}.name = 'AC2D4_new';
% el_type_info{i}.dims = 2;
% el_type_info{i}.state = 'fluid';
% el_type_info{i}.shape = 'quadrilateral';

%2D modelling - fluid-solid interface
i = i + 1;
el_type_info{i}.name = 'ASI2D2';
el_type_info{i}.dims = 2;
el_type_info{i}.state = 'fluid_solid_interface';
el_type_info{i}.shape = 'line';

%3D modelling - solids
i = i + 1;
el_type_info{i}.name = 'C3D4';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'tetrahedral';

i = i + 1;
el_type_info{i}.name = 'C3D6';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'triangular prism';

i = i + 1;
el_type_info{i}.name = 'C3D8';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'hexahedral';

i = i + 1;
el_type_info{i}.name = 'C3D8R';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'solid';
el_type_info{i}.shape = 'hexahedral';

%3D modelling - fluids
i = i + 1;
el_type_info{i}.name = 'AC3D4';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'tetrahedral';

i = i + 1;
el_type_info{i}.name = 'AC3D8';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'hexahedral';

i = i + 1;
el_type_info{i}.name = 'AC3D8R';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'fluid';
el_type_info{i}.shape = 'hexahedral';

%3D modelling - fluid-solid interface
i = i + 1;
el_type_info{i}.name = 'ASI3D3';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'fluid_solid_interface';
el_type_info{i}.shape = 'triangular';

i = i + 1;
el_type_info{i}.name = 'ASI3D4';
el_type_info{i}.dims = 3;
el_type_info{i}.state = 'fluid_solid_interface';
el_type_info{i}.shape = 'quadrilateral';

%Add the face indices based on underlying shape
el_type_info = fn_faces(el_type_info);
end


function el_type_info = fn_faces(el_type_info)
for i = 1:numel(el_type_info)
    switch el_type_info{i}.shape
        case 'triangular' %2D triangular
            fc_i = [
                1,2
                2,3
                3,1];
        case 'quadrilateral' %2d quadrilaterals
            fc_i = [
                1,2
                2,3
                3,4
                4,1];
        case 'line' %2D line
            fc_i = [
                1,2];
        case 'tetrahedral' %3D tetrahedral
            fc_i = [
                1,2,3
                1,2,4
                2,3,4
                1,3,4];
        case 'hexahedral' %3D hexahedral
            fc_i = [
                1,2,3,4
                1,2,6,5
                2,3,7,6
                3,4,8,7
                4,1,5,8
                5,6,7,8
                ];
    end
    el_type_info{i}.faces = fc_i;
end
end