function [x, u, w_hat] = simulate_state_fdbk_per_timestep(sys, ctrller, params, w_hat, wt, x0, t)
params.sanity_check();

Rc = ctrller.Rc_;
Mc = ctrller.Mc_;
T  = length(Rc);

u     = zeros(sys.Nu, 1);
x_hat = zeros(sys.Nx, 1);

for k=1:1:min(t-1, T)
   u = u + Mc{k}*w_hat(:,t-k);
end

for k=1:1:min(t-1, T-1)
   x_hat = x_hat + Rc{k+1}*w_hat(:,t-k);       
end 
    
x = sys.A*x0 + sys.B1*wt+ sys.B2*u;
w_hat(:,t) = x - x_hat;

end
