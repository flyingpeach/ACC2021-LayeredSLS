function [delta, PG] = solve_opf(B, PL, delta_min, delta_max, PG_min, PG_max)
% This function solves the OPF problem, as stated in Frank et al. In
% particular, we minimize (5) subject to (8), (11) and (72) (in replacement
% of (6)). 

dim = max(size(PG_min)); 

%% Define the cost function

Q = eye(dim+1); p = ones(dim+1,1); % We can change this. According to Frank's, Q needs to be diagonal.

%% Solve the OPF

cvx_begin

    variable PG(dim+1, 1)
    variable delta(dim+1, 1)

    minimize (PG'*Q*PG + p'*PG) % Eq (5) with quadratic cost

    subject to 

    PG_min(1:dim) <= PG(1:dim); PG(1:dim) <= PG_max(1:dim); % Eq (8) (this bound is only applicable to the components that sit in set G, which is a subset of the entire network (G\subset N as in the paper). So we need to include that here as an input from the plant tolopology, I think)
    delta_min(1:dim) <= delta(1:dim); delta(1:dim) <= delta_max(1:dim); % Eq (11)
    for i = 1:dim
        PG(i) == B(i,:)*(delta(i)*ones-delta(1:dim)) + PL(i); % Eq (72)
    end
    delta(dim+1) == 0; % Slack bus

cvx_end

%% Don't return slack bus
delta = delta(1:dim);
PG    = PG(1:dim);