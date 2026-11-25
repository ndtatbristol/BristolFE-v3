function [cL,cS] = fn_vels_from_matl_def(matl)
%This recovers velocities from fluid and isotropic solid material
%descriptions. Useful when injecting fields based on analytical plane wave models

if all(size(matl.D) == 6)

    tol = 1e-6;

    if abs(matl.D(1,1) - (matl.D(1,2) + 2*matl.D(4,4))) > tol*abs(matl.D(1,1))
        warning('Stiffness matrix is not isotropic.');
    end

    lambda = mean([C(1,2), C(1,3), C(2,3)]);
    mu     = mean([C(4,4), C(5,5), C(6,6)]);

    cL = sqrt((lambda + 2*mu)/matl.rho);
    cS = sqrt(mu/rho);

elseif all(size(matl.D) == 1)
    cL = sqrt(matl.D/matl.rho);
    cS = [];
end

end