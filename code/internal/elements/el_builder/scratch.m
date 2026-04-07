%Current AC2D3 terms

%Shape functions
n(1) = q(1);
n(2) = q(2);
n(3) = 1 - q(1) - q(2);

%Determinant of Jacobian, OK
detJ = nds_1_1*nds_2_2 - nds_1_2*nds_2_1 - nds_1_1*nds_3_2 + nds_1_2*nds_3_1 + nds_2_1*nds_3_2 - nds_2_2*nds_3_1
XXX =  nds_1_1*nds_2_2 - nds_1_2*nds_2_1 - nds_1_1*nds_3_2 + nds_1_2*nds_3_1 + nds_2_1*nds_3_2 - nds_2_2*nds_3_1

%Inverse Jacobian (ignoring 1/detJ factor)
invJ_times_detJ = [
    nds_2_2 - nds_3_2, nds_3_1 - nds_2_1
    nds_3_2 - nds_1_2, nds_1_1 - nds_3_1]

%Final B-matrix
B = [ 
     (nds_2_2 - nds_3_2)/J, -(nds_1_2 - nds_3_2)/J, (nds_1_2 - nds_3_2)/J - (nds_2_2 - nds_3_2)/J
    -(nds_2_1 - nds_3_1)/J,  (nds_1_1 - nds_3_1)/J, (nds_2_1 - nds_3_1)/J - (nds_1_1 - nds_3_1)/J]
                         0,                      0,                                             0]

XXX=[
     (nds_2_2 - nds_3_2)/detJ, -(nds_2_1 - nds_3_1)/detJ, (nds_2_1 - nds_3_1)/detJ - (nds_2_2 - nds_3_2)/detJ
    -(nds_1_2 - nds_3_2)/detJ,  (nds_1_1 - nds_3_1)/detJ, (nds_1_2 - nds_3_2)/detJ - (nds_1_1 - nds_3_1)/detJ
                        0,                         0,                                                   0]