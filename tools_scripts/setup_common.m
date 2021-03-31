% System
gridSize       = 5;
numNodes       = gridSize * gridSize;
actDens        = 1;
numActs        = round(actDens*numNodes);
actuatedNodes  = randsample(numNodes, numActs);

connectThresh  = 0.65; % lower = more connected grid
locality       = 3;
commspeed      = 2;

% Bounds for DC-OPF problem
theta_min = -pi;
theta_max = pi;
PG_min    = 0;
PG_max    = 25;

% Saturation limits
uMax      = 9;

% Simulation parameters
tOPFStart = 10;  % when to solve the first OPF
tOPF      = 100; % how often to solve OPF
tMPC      = 20; % how often to solve MPC, must divide evenly into tOPF

Ts = 0.2; % sampling time

% Sanity check
if mod(tOPF, tMPC)
    error('tMPC must divide evenly into tOPF');
end

%% Setup power grid and systems
[adjMtx, nodeCoords, susceptMtx, inertiasInv, dampings] = setup_topology(gridSize, connectThresh);
sys = setup_system(actuatedNodes, adjMtx, susceptMtx, inertiasInv, dampings, Ts);
