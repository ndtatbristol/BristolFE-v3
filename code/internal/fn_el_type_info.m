function el_type_info = fn_el_type_info()

i = 1;
el_type_info(i).name = 'CPE3';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = i + 1;
el_type_info(i).name = 'CPS3';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = i + 1;
el_type_info(i).name = 'AC2D3';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid';

i = i + 1;
el_type_info(i).name = 'AC2D4';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid';

i = i + 1;
el_type_info(i).name = 'AC2D4R';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid';

i = i + 1;
el_type_info(i).name = 'ASI2D2';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid_solid_interface';

i = i + 1;
el_type_info(i).name = 'CPE4';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = i + 1;
el_type_info(i).name = 'CPE4R';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = i + 1;
el_type_info(i).name = 'CPS4';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = i + 1;
el_type_info(i).name = 'CPS4R';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';
end