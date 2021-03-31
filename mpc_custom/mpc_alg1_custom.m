function [x, u, time, iters] = mpc_alg1_custom(sys, x0, params, tau_i, tau_d, mu)
% This function differs from the toolbox mpc_distributed in the following
% ways:
%   It returns the entire future trajectory + inputs
%   It uses adaptive ADMM tuning for rho (params are tau_i, tau_d, mu)
%   It includes only mechanisms for algorithm 1
%   It ignores solverMode in the params and defaults to fastest solver

%% Setup
sanity_check_actuation(sys)
params.sanity_check_alg_1();

% For ease of notation
Nx = sys.Nx; Nu = sys.Nu;
tFIR = params.tFIR_;

maxIters     = params.maxIters_;

nVals = Nx*tFIR + Nu*(tFIR-1);

% Track runtime and iterations
times        = zeros(Nx, 1); % per state

% Get indices corresponding to rows / columns / localities
PsiSupp  = get_psi_sparsity(sys, params); % Toeplitz matrix
PhiSupp  = PsiSupp(:, 1:Nx);              % First block column
r   = assign_rows_phi(sys, tFIR);
c   = assign_cols_phi(sys);
s_r = get_row_locality(PhiSupp);
s_c = get_col_locality(PhiSupp);

% Constraints and costs
C     = build_cost_mtx(params);
K     = build_constr_mtx(sys, params);

% ADMM variables
Phi    = zeros(nVals, Nx);
Psi    = zeros(nVals, Nx);
Lambda = zeros(nVals, Nx);

% Precalculate items for column-wise update
[zabs, eyes, zabis] = precalculate_col(sys, tFIR, s_c);

%% MPC
for iters=1:maxIters % ADMM (outer loop)
    rho      = params.rho_; % in case of adaptive ADMM
    Psi_prev = Psi;

    % Step 4: Solve for Phi
    Phi_rows = cell(nVals, 1);
    
    % Solve for Phi for uncoupled rows
    for i = 1:Nx
        for j = 1:length(r{i})
            row   = r{i}{j};
            x_loc = x0(s_r{row}); % observe local state
            cost_ = C(row, row);

            if K(row, row) % has constraint
                if row <= tFIR*Nx % state constraint
                    b1_ = params.stateUB_(i) / K(row, row);
                    b2_ = params.stateLB_(i) / K(row, row);
                else % input constraint
                    inputIdx = find(sys.B2(i,:));
                    b1_ = params.inputUB_(inputIdx) / K(row, row);
                    b2_ = params.inputLB_(inputIdx) / K(row, row);
                end
                b1  = max(b1_,b2_); b2 = min(b1_,b2_); % in case of negative signs
                solverMode = MPCSolverMode.Explicit;
            else % no constraint, use closed form
                solverMode = MPCSolverMode.ClosedForm;
            end

            if solverMode == MPCSolverMode.ClosedForm
                tic;
                Phi_rows{row} = eqn_16a_closed(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), cost_, rho);
                times(i) = times(i) + toc;
            else % use explicit form
                tic;
                Phi_rows{row} = eqn_16a_explicit(x_loc, Psi(row, s_r{row}), Lambda(row, s_r{row}), b1, b2, cost_, rho);
                times(i) = times(i) + toc;
            end
        end
    end
    
    % Note: skipped a bunch of steps (omitted since no coupling here)
    
    % Step 10: Build entire Phi matrix
    Phi = build_from_rows(r, s_r, Phi_rows, size(Phi));
        
    % Step 11: Solve (16b) to get local Psi
    Psi_cols = cell(Nx, 1);
    for i = 1:Nx
        tic;
        Psi_cols{i} = eqn_16b(Phi(s_c{i}, c{i}{1}), Lambda(s_c{i}, c{i}{1}), zabs{i}, eyes{i}, zabis{i});
        times(i) = times(i) + toc;        
    end
    
    % Step 12: Build entire Psi matrix
    Psi = build_from_cols(c, s_c, Psi_cols, size(Psi));
    
    % Step 13: Update Lambda
    Lambda = Lambda + Phi - Psi;
    
    % Step 14: Check convergence of ADMM (outer loop)
    converged = true;
    for i = 1:Nx
        phi_      = [];
        psi_      = [];
        psi_prev_ = [];
        for j = 1:length(r{i})
            % Due to dimensionality issues, not stacking rows
            % Instead, just make one huge row
            % (since we're checking Frob norm, doesn't matter)
            row = r{i}{j};
            phi_      = [phi_, Phi(row, s_r{row})];
            psi_      = [psi_, Psi(row, s_r{row})];
            psi_prev_ = [psi_prev_, Psi_prev(row, s_r{row})];
        end
        
        [conv, scale] = check_convergence_adaptive(phi_, psi_, psi_, psi_prev_, ...
                                                   params, tau_i, tau_d, mu);
        if ~conv
            converged   = false;
            params.rho_ = params.rho_ * scale;
            break; % if one fails, can stop checking the rest
        end
    end
    
    if converged
        break; % exit ADMM iterations
    end
end

if ~converged
    fprintf('ADMM reached %d iters and did not converge\n', maxIters);
end

time = mean(times);

% Compute control + state
u = zeros(Nu, tFIR);
x = zeros(Nx, tFIR);

% Split Phi for ease of notation
Phix = Phi(1:Nx*tFIR, :);
Phiu = Phi(Nx*tFIR+1:end, :);

x(:,1) = x0;

for k=1:tFIR
    if k>1
        x_start = (k-1)*Nx+1;
        x_end   = k*Nx;
        x(:,k)  = Phix(x_start:x_end,:)*x0;
    end
    
    if k<tFIR
        u_start = (k-1)*Nu+1;
        u_end   = k*Nu;
        u(:,k)  = Phiu(u_start:u_end,:)*x0;
    end
end

end