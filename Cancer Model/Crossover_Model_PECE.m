function [t_all, Y_all] = Crossover_Model_PECE(t1, t2, T, h, Y0, C, alpha_const, alphaFun)
% piecewise integrate:
% 0 : t1   -> integer (ode45)
% t1 : t2  -> fractional Caputo (fracCaputoPECE) with order alpha_const
% t2 : T   -> variable-order (voFracPECE) using alphaFun

% 1) integer part
tspan1 = [0 t1];
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t1_vec, Y1] = ode45(@(t,y) Cancer_Model(t,y,C), tspan1, Y0, opts);

% prepare initial cond at t1
Y_t1 = Y1(end,:)';

% 2) fractional part (we integrate fractional with initial condition at t1)
% we will shift time so frac solver starts at 0 with Y0 = Y_t1
T2 = t2 - t1;
[t_f, Yf] = Frac_Caputo_PECE(alpha_const, T2, h, Y_t1, C, @(tt,yy,prms) Cancer_Model(tt + t1, yy, prms));

% adjust times to absolute
t_f = t_f + t1;

% 3) variable-order part: start from last point of fractional at t2
Y_t2 = Yf(end,:)';
T3 = T - t2;
% define shifted alphaFun for absolute time
alphaFunShifted = @(tau) alphaFun(tau + t2); % or just use alphaFun directly with absolute time
[t_v, Yv] = VO_Frac_PECE(alphaFunShifted, T3, h, Y_t2, C, @(tt,yy,prms) Cancer_Model(tt + t2, yy, prms));

t_v = t_v + t2;

% assemble
t_all = [t1_vec; t_f(2:end); t_v(2:end)]; % avoid duplicate endpoints
Y_all = [Y1; Yf(2:end,:); Yv(2:end,:)];
end
