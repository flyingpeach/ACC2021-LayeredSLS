numRuns = 30;

unsatCenLinCosts = zeros(numRuns, 1);
cenLinCosts      = zeros(numRuns, 1);
cenMPCCosts      = zeros(numRuns, 1);
locLayeredCosts  = zeros(numRuns, 1);

for i=1:numRuns
    fprintf('\n*******************************\n');
    fprintf('\nStarting run %d of %d\n', i, numRuns);
    fprintf('\n*******************************\n');
    seed = 100 + i;
    if i >= 8
        % Seed 108 gives a weird plant which causes issues for all methods
        % So do not use it
        seed = 101 + i;
    end
    [unsatCenLinCosts(i), cenLinCosts(i), cenMPCCosts(i), ...
          locLayeredCosts(i)] = example_multi_per_run(seed);   
end

%% Calculate costs and print
baseline = unsatCenLinCosts;

% Normed costs
unsatCenLinNorm = unsatCenLinCosts ./ baseline;
cenLinNorm      = cenLinCosts ./ baseline;
cenMPCNorm      = cenMPCCosts ./ baseline;
locLayeredNorm  = locLayeredCosts ./ baseline;

fprintf('\nMean normalized costs:\n')
fprintf('UnsatCenLin : %.3f\n', mean(unsatCenLinNorm));
fprintf('CenLin      : %.3e\n', mean(cenLinNorm));
fprintf('CenMPC      : %.3f\n', mean(cenMPCNorm));
fprintf('LocLayered  : %.3f\n', mean(locLayeredNorm));

%% Separate stable and unstable cases

% We are saying that unstable cases are ones in which
% the cost increases more than 10-fold compared to standard
% Based on observing plots this is reasonable

unstableIdx = find(cenLinNorm > 10);
stableIdx   = setdiff(1:numRuns, unstableIdx);

fprintf('\nMean normalized costs, stable only:\n')
fprintf('UnsatCenLin : %.3f\n', mean(unsatCenLinNorm(stableIdx)));
fprintf('CenLin      : %.3f\n', mean(cenLinNorm(stableIdx)));
fprintf('CenMPC      : %.3f\n', mean(cenMPCNorm(stableIdx)));
fprintf('LocLayered  : %.3f\n', mean(locLayeredNorm(stableIdx)));

fprintf('\nMean normalized costs, unstable only:\n')
fprintf('UnsatCenLin : %.3f\n', mean(unsatCenLinNorm(unstableIdx)));
fprintf('CenLin      : %.3e\n', mean(cenLinNorm(unstableIdx)));
fprintf('CenMPC      : %.3f\n', mean(cenMPCNorm(unstableIdx)));
fprintf('LocLayered  : %.3f\n', mean(locLayeredNorm(unstableIdx)));



