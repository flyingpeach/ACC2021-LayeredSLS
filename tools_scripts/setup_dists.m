%% Randomly generated disturbances
wsDist = zeros(sys.Nw, tSim);

% smaller disturbances
distSize     = 0.5;                % max disturbance size
numDistNodes = round(numNodes/1);  % how many nodes to disturb
numDistTimes = round(tSim/4);      % how many times to disturb at

% bigger disturbances
faultSizeMean = 1.5;
numFaultNodes = round(numNodes/4);
numFaultTimes = 1;

distNodes = randsample(numNodes, numDistNodes);
faultNodes = randsample(numNodes, numFaultNodes);

for i=1:length(distNodes)
    freqDist = distNodes(i)*2; % convert index    
    
    distTimes = randsample(tSim-10, numDistTimes);
    wsDist(freqDist, distTimes) = rand(numDistTimes, 1) * distSize * 2 - distSize;    
end

for i=1:length(faultNodes)
    freqDist = faultNodes(i)*2; % convert index
    
    faultTimes = randsample(tSim-10, numFaultTimes);
    wsDist(freqDist, faultTimes) = rand(numFaultTimes, 1) + faultSizeMean;
end