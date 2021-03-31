function [converged, scale] = check_convergence_adaptive(prim1, prim2, dual1, dual2, ...
                                                         params, tau_i, tau_d, mu)
% Checks whether the following is satisfied, where ||  ||_F is frob norm
% || prim1 - prim2 ||_F <= params.eps_p_
% || dual1 - dual2 ||_F <= params.eps_d_

% Return scaling factor for adaptive ADMM rho
% mu, tau_i, and tau_d relate to adaptive ADMM only

primRes = norm(prim1 - prim2, 'fro');
dualRes = norm(dual1 - dual2, 'fro');
converged = primRes <= params.eps_p_ && dualRes <= params.eps_d_;

scale = 1;
if primRes > mu * dualRes
    scale = tau_i;
elseif dualRes > mu * primRes
    scale = 1 / tau_d;
end

end
