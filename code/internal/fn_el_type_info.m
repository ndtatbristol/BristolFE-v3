function el_type_info = fn_el_type_info()

i = 1;
el_type_info(i).name = 'CPE3';
el_type_info(i).dims = 2;
el_type_info(i).state = 'solid';

i = 2;
el_type_info(i).name = 'AC2D3';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid';

i = 3;
el_type_info(i).name = 'ASI2D2';
el_type_info(i).dims = 2;
el_type_info(i).state = 'fluid_solid_interface';

end