function [force_inp, force_set, force_in_set, force_out_set] = fn_convert_disps_to_forces_v2(K_sub, C_sub, M_sub, time_step, disp_inp, lyrs, in_or_out, solver_mode)
% C_sub(:) = 0;%this is to fix a problem with C export from Pogo only!
force_in_set = lyrs == 2;
force_out_set = lyrs == 3;
force_set = force_in_set | force_out_set;
switch in_or_out
    case 'in'
        freeze_set = [lyrs == 3 | lyrs == 4];
    case 'out'
        freeze_set = [lyrs == 1 | lyrs == 2];
end
disp_inp(freeze_set, :) = 0;
tmp = [zeros(size(disp_inp, 1), 2), disp_inp];
accn = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / time_step ^ 2;
switch lower(solver_mode)
    case {'vel at last half time step', 'explicit', 'exp'}
        vel = (tmp(:, 2:end - 1) - tmp(:, 1:end - 2)) / time_step;
    case {'vel at curent time step', 'implicit', 'imp'}
        vel = (tmp(:, 3:end) - tmp(:, 1:end - 2)) / (2 * time_step);
    case {'predictor corrector', 'pc'}
        I_sub = speye(size(K_sub));
        M_sub = spdiags(sum(M_sub).', 0, size(M_sub,1), size(M_sub,2));
        inv_M_sub = M_sub \ speye(size(M_sub));
        vel = (I_sub - time_step / 2 * C_sub * inv_M_sub) \ ((tmp(:, 2:end - 1) - tmp(:, 1:end - 2)) / time_step - C_sub * inv_M_sub * K_sub * tmp(:, 2:end - 1) * time_step / 2);
        % error('Subdomain disp to force mapping not implemented for predictor-corrector solver yet');
        % vel = (tmp(:, 2:end - 1) - tmp(:, 1:end - 2)) / time_step + time_step * M_sub \ (f - K * tmp(:, 2:end - 1));
end
disp = [zeros(size(disp_inp, 1), 1), disp_inp(:, 1:end - 1)];
% disp = [zeros(size(disp_inp, 1), 2), disp_inp(:, 1:end - 2)];
% disp = disp_inp;
force_inp  = ...
        M_sub(force_set, :) * accn + ...
        C_sub(force_set, :) * vel + ...
        K_sub(force_set, :) * disp;


end