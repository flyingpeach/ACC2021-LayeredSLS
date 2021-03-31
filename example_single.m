rng(215);

tSim = 60;
setup_common;
setup_ctrllers;
setup_opf;
setup_dists;

% Add a big disturbance for kicks
bigDistNode = 17;
bigDistFreq = bigDistNode*2;
wsDist(bigDistFreq, 40) = rand() + faultSizeMean;

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

fprintf('\nNormalized costs:\n')
fprintf('UnsatCenLin : %.3f\n', unsatCenLinCost / baseline);
fprintf('CenLin      : %.3f\n', cenLinCost / baseline);
fprintf('CenMPC      : %.3f\n', cenMPCCost / baseline);
fprintf('LocLayered  : %.3f\n', locLayeredCost / baseline);

%% Plotting
% Note: this is distinct from figure plotting for the paper

% These are the nodes where saturation makes a difference
% Only for plotting purposes
badStates = find(max(xsUnsatCenLin - xsCenLin, [] ,2) > 1);
badNodes  = unique(ceil(badStates/2));

plotNode  = 16;
phaseIdx  = plotNode*2-1;
freqIdx   = plotNode*2;
timeVals  = (1:tSim) * Ts;

figure();
subplot(4,1,1); title('Phase'); hold on;
plot(timeVals, xsCenLin(phaseIdx, :) + thetas(plotNode, :));
plot(timeVals, xsCenMPC(phaseIdx, :) + thetas(plotNode, :));
plot(timeVals, xsLocLayered(phaseIdx, :) + thetas(plotNode, :));
plot(timeVals, thetas(plotNode, :));
legend('CenLin', 'CenMPC', 'LocLayered', 'Setpoint');
ylim([-10 20]);

subplot(4,1,2); title('Frequency'); hold on;
plot(timeVals, xsCenLin(freqIdx, :));
plot(timeVals, xsCenMPC(freqIdx, :));
plot(timeVals, xsLocLayered(freqIdx, :));
ylim([-40 30]);

% Be careful: assumes there is actuation at this node
subplot(4,1,3); title('Actuation'); hold on;
plot(timeVals, usCenLin(plotNode, :));
plot(timeVals, usCenMPC(plotNode, :));
plot(timeVals, usLocLayered(plotNode, :));
xlabel('Time');

subplot(4,1,4); title('Disturbances'); hold on;
plot(timeVals, wsDist(freqIdx, :));
