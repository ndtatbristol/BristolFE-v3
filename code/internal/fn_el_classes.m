function el_class = fn_el_classes()
%Utility function to keep track of types of elements to save having to pass
%them around all the time


el_class.solid = {'CPE3'}; 
el_class.fluid = {'AC2D3'};
el_class.fluid_solid_interface = {'ASI2D2'};
end