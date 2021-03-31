function [x, u] = simulate_layered(sys, ctrller, simParams, ...
                                   mpcParams, tau_i, tau_d, mu, ...                                            
                                   wsDist, tSim, tMPC, tOPFStart, uMax)          
x      = zeros(sys.Nx, tSim);
u      = zeros(sys.Nu, tSim);

x_err  = zeros(sys.Nx, tSim); % error from trajectory
u_err  = zeros(sys.Nu, tSim);
xTraj  = zeros(sys.Nx, tSim); % suggested trajectories by MPC
uTraj  = zeros(sys.Nu, tSim); % not used in simulation
w_hat  = zeros(sys.Nw, tMPC); 

for t=1:tSim-1
    fprintf('Layered: Calculating time %d of %d\n', t+1, tSim);
        
    % Time to run new mid-tier MPC
    if ~mod(t-tOPFStart-1, tMPC)
        mpcParams.rho_ = 1e1;
            
        xTraj(:,t) = x(:,t);
        [xs, us, ~] = mpc_alg1_custom(sys, xTraj(:,t), mpcParams, tau_i, tau_d, mu);
        for tm=0:min(mpcParams.tFIR_-2, tSim-t-1)
            t_ = t + tm;
               
            xTraj(:, t_) = xs(:, tm+1);
            uTraj(:, t_) = us(:, tm+1);
        end            
    end

    if ~mod(t-tOPFStart-1, tMPC)
        % Reset error + offline simulation (since traj changed)
        x_err(:,t) = zeros(sys.Nx, 1);
        w_hat      = zeros(sys.Nw, tMPC);
    end
        
    err = x(:, t) - xTraj(:, t);
    wt  = err - x_err(:, t) + wsDist(:,t);        
    to  = mod(t-tOPFStart-1, tMPC) + 1; % offline time idx
        
    [x_err(:,t+1), u_err(:,t), w_hat] = simulate_state_fdbk_per_timestep(sys, ctrller, simParams, w_hat, wt, x_err(:, t), to);
       
    u(:, t)   = uTraj(:, t) + u_err(:, t);        
    % Saturate the output
    u(:,t)    = sat(u(:,t), uMax);
    x(:, t+1) = sys.A*x(:, t) + sys.B2*u(:, t) + sys.B1*simParams.w_(:,t);
end