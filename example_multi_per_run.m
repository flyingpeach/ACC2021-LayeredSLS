function [unsatCenLinCost, cenLinCost, cenMPCCost, locLayeredCost] = example_multi_per_run(seed)
% Differences compared to example_single:
% tSim = 40 instead of tSim = 60 (the last 20 seconds are low-norm anyway)
% Lack of extra plot-specific disturbance (since we are not plotting)
% No plotting

rng(seed);

tSim = 40;
setup_common;
setup_ctrllers;
setup_opf;
setup_dists;

%% Simulate offline controllers
simParams       = SimParams();
simParams.tSim_ = tSim; 
simParams.w_    = wsOPF + wsDist;

xsUnsatCenLin = zeros(sys.Nx, tSim);
usUnsatCenLin = zeros(sys.Nu, tSim);
xsCenLin      = zeros(sys.Nx, tSim);
usCenLin      = zeros(sys.Nu, tSim);

for t=1:tSim-1
    usUnsatCenLin(:, t)   = -KOpt*xsUnsatCenLin(:, t); 
    xsUnsatCenLin(:, t+1) = sys.A*xsUnsatCenLin(:,t) + sys.B2*usUnsatCenLin(:,t) + sys.B1*simParams.w_(:,t);

    usCenLin(:, t)   = -KOpt*xsCenLin(:,t);
    usCenLin(:, t)   = sat(usCenLin(:,t), uMax);
    xsCenLin(:, t+1) = sys.A*xsCenLin(:,t) + sys.B2*usCenLin(:,t) + sys.B1*simParams.w_(:,t);  
end

%% Simulate global MPC
xsCenMPC = zeros(sys.Nx, tSim);
usCenMPC = zeros(sys.Nu, tSim);

for t=1:tSim-1       
    fprintf('CenMPC: Calculating time %d of %d\n', t+1, tSim);
    [x, u] = mpc_centralized(sys, xsCenMPC(:,t), mpcParamsGlob);
        
    usCenMPC(:,t)   = sat(u, uMax);
    xsCenMPC(:,t+1) = x + sys.B1*simParams.w_(:,t);
end

%% Simulate layered case
[xsLocLayered, usLocLayered] = simulate_layered(sys, ctrllerLr, simParams, ...
                                                mpcParamsLr, tau_i, tau_d, muAdapt, ...
                                                wsDist, tSim, tMPC, tOPFStart, uMax);

%% Calculate performance
unsatCenLinCost = get_cost_fn(mpcParamsGlob, xsUnsatCenLin, usUnsatCenLin);
cenLinCost      = get_cost_fn(mpcParamsGlob, xsCenLin, usCenLin);
cenMPCCost      = get_cost_fn(mpcParamsGlob, xsCenMPC, usCenMPC);
locLayeredCost  = get_cost_fn(mpcParamsGlob, xsLocLayered, usLocLayered);

baseline = unsatCenLinCost;

fprintf('\nNormalized costs this run:\n')
fprintf('UnsatCenLin : %.3f\n', unsatCenLinCost / baseline);
fprintf('CenLin      : %.3f\n', cenLinCost / baseline);
fprintf('CenMPC      : %.3f\n', cenMPCCost / baseline);
fprintf('LocLayered  : %.3f\n', locLayeredCost / baseline);
