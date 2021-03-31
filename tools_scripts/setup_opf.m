thetas = zeros(numNodes, tSim); % OPF set points
wsOPF  = zeros(sys.Nw, tSim);

theta_mins = theta_min * ones(numNodes, 1);
theta_maxs = theta_max * ones(numNodes, 1);
PG_mins    = PG_min    * ones(numNodes, 1);
PG_maxs    = PG_max    * ones(numNodes, 1);

% Pre-generate OPF values
for t=1:tSim-1
    if ~mod(t-tOPFStart, tOPF)  
        % generate random load profile
        PL                  = rand(numNodes, 1) * PG_max;
        [thetas(:,t+1), PG] = solve_opf(susceptMtx, PL, theta_mins, theta_maxs, PG_mins, PG_maxs);

        w_theta          = thetas(:,t) - thetas(:,t+1);
        wsOPF(1:2:end,t) = w_theta;
    else
        thetas(:,t+1) = thetas(:,t);
    end
end