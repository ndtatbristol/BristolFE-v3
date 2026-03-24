function [K, M] = fn_get_pogo_matrices(el_type, nds, els, D, rho)
%returns K and M matrices for element of type specified for checking
%against BristolFE equivalent
mod.nds = nds;
mod.els = els;
n_els = size(els, 1);
mod.el_mat_i = ones(n_els, 1);
mod.el_typ_i = ones(n_els, 1);;
el_types = {el_type};

matls{1}.D = D;
matls{1}.rho = rho;


steps{1}.load.time = linspace(0,1e-3,100);%dummy
steps{1}.load.frc_nds = [1];%dummy
steps{1}.load.frc_dfs = [1];%dummy
steps{1}.load.frcs = zeros(size(steps{1}.load.time));%dummy
steps{1}.mon.nds = [1];
steps{1}.mon.dfs = [1];

fe_options.dof_to_use = [1,2,3];
fe_options.field_output_every_n_frames = inf;
fe_options.solver_precision = 'double';
fe_options.pogo_verbosity = 2;
[~, mats] = fn_pogoFE(mod, matls, el_types, steps, fe_options);
K = mats.K;
M = mats.M;
end