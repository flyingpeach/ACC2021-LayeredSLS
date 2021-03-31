%% Controller for UnsatCentLin and CentLin
A = full(sys.A);
B = full(sys.B2);
QSqrt = full(sys.C1(1:sys.Nx,:));
RSqrt = full(sys.D12(sys.Nx+1:end, :));

Q = QSqrt * QSqrt;
R = RSqrt * RSqrt;

[S, ~] = idare(A, B, Q, R);
KOpt = inv(B'*S*B + R)*B'*S*A;

%% Controller for Global MPC
mpcParamsGlob           = MPCParams();
mpcParamsGlob.locality_ = 1000; % Effectively global
mpcParamsGlob.tFIR_     = tMPC; % Otherwise worse than 3-layer

mpcParamsGlob.QSqrt_ = sys.C1(1:sys.Nx,:);
mpcParamsGlob.RSqrt_ = sys.D12(sys.Nx+1:end, :);

mpcParamsGlob.maxIters_ = 10000;
mpcParamsGlob.eps_p_    = 1e-3;
mpcParamsGlob.eps_d_    = 1e-3;

mpcParamsGlob.inputConsMtx_ = eye(sys.Nu);
mpcParamsGlob.inputUB_      = uMax * ones(sys.Nu, 1);
mpcParamsGlob.inputLB_      = -mpcParamsGlob.inputUB_;

%% Bottom layer controller for LocLayered
slsParamsLr       = SLSParams();
slsParamsLr.T_    = 5;

slsParamsLr.add_objective(SLSObjective.H2, 1);
slsParamsLr.add_constraint(SLSConstraint.CommSpeed, commspeed);
slsParamsLr.add_constraint(SLSConstraint.Locality, locality);

clMapsLr  = state_fdbk_sls(sys, slsParamsLr);
ctrllerLr = Ctrller.ctrller_from_cl_maps(clMapsLr);

%% Top layer controller for LocLayered
% Adaptive parameters
tau_i   = 2;
tau_d   = 2;
muAdapt = 10;

mpcParamsLr           = copy(mpcParamsGlob);
mpcParamsLr.locality_ = locality; % Make it local

% Tuning for our particular example to converge
mpcParamsLr.maxIters_ = 10000;
mpcParamsLr.eps_p_    = 3e-4;
mpcParamsLr.eps_d_    = 3e-4;
