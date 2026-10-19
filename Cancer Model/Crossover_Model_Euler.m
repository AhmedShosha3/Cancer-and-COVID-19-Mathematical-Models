function [t_all, Y_all] = Crossover_Model_Euler(t1, t2, T, h, Y0, C, alpha_const, alphaFun)

% Piecewise integrate:
% 0  : t1   -> integer (cancerSystemEuler)
% t1 : t2  -> fractional Caputo (fracCaputoEuler) with order alpha_const
% t2 : T   -> variable-order (voFracEuler) using alphaFun
%
% Inputs same as original crossoverModel2. Uses step h as dt for Euler solvers.

% 1) integer part using your cancerSystemEuler (explicit Euler)
tspan1 = [0, t1];
% cancerSystemEuler signature: [Tvec, X] = cancerSystemEuler(X0, C, tspan, dt)
[t1_vec, Y1] = Cancer_Model_Euler(Y0', C, tspan1, h); % note Y0' because cancerSystemEuler expects row X0

% ensure column convention
Y_t1 = Y1(end, :)'; % column

% 2) fractional part (constant-order) using your fracCaputoEuler
% We shift time so fractional solver starts at 0 with Y0 = Y_t1
T2 = t2 - t1;
% fracCaputoEuler(alpha, T, h, Y0, params, rhsFun)
rhs_shift_f = @(tt, yy, prms) Cancer_Model(tt + t1, yy, prms); % tt is local (0..T2)
[t_f, Yf] = Frac_Caputo_Euler(alpha_const, T2, h, Y_t1, C, rhs_shift_f);

% adjust times to absolute
t_f = t_f + t1;

% 3) variable-order part: start from last point of fractional at t2
Y_t2 = Yf(end, :)';
T3 = T - t2;

% define shifted alphaFun for local time within vo solver (voFracEuler expects alphaFun(local_time))
% We'll pass alphaFunShifted which maps local tau (0..T3) to absolute time alpha
alphaFunShifted = @(tau) alphaFun(tau + t2);

% rhs for variable-order part: local time tt maps to absolute time tt+t2 inside cancerSystem
rhs_shift_v = @(tt, yy, prms) Cancer_Model(tt + t2, yy, prms);

% call variable-order Euler solver
[t_v, Yv] = VO_Frac_Euler(alphaFunShifted, T3, h, Y_t2, C, rhs_shift_v);

% shift times to absolute
t_v = t_v + t2;

% assemble full solution (avoid duplicate endpoints)
t_all = [t1_vec(:); t_f(2:end); t_v(2:end)];
Y_all = [Y1; Yf(2:end, :); Yv(2:end, :)];

end
