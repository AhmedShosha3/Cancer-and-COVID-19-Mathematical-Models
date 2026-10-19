function [t, Y] = VO_Frac_Euler(alphaFun, T, h, Y0, params, rhsFun)
% voFracEuler  Variable-order Caputo fractional Euler-like solver
%
% Solves D^{alpha(t)} y(t) = f(t,y), 0 < alpha(t) <= 1, with y(0) = Y0
% using an explicit convolutional scheme where coefficients at step n use
% the local order alpha_n = alphaFun(t_n).
%
% Inputs:
%   alphaFun - function handle alphaFun(t) giving alpha at absolute time t (here t is local 0..T)
%   T        - final local time (solver runs from 0..T)
%   h        - time step
%   Y0       - initial state vector (column or row)
%   params   - additional parameters passed to rhsFun
%   rhsFun   - function handle f = rhsFun(t, y, params), t is local time
%
% Outputs:
%   t - (N+1)x1 time vector (local 0..T)
%   Y - (N+1) x dim solution matrix (rows correspond to times)
%
% NOTE: This is a straightforward variable-order generalization of the
% fractional-Euler convolution form. It's a practical/approximate scheme.

% input checks
if h <= 0 || T < 0
    error('h must be > 0 and T >= 0');
end

% ensure Y0 column
Y0 = Y0(:);
dim = numel(Y0);

N = floor(T / h);
t = (0:N)' * h;
Y = zeros(N+1, dim);
Y(1, :) = Y0.';  % row

% store f evaluations
F = zeros(N+1, dim);
F(1, :) = rhsFun(t(1), Y(1,:)', params)';

% time-stepping
for n = 1:N
    % local time at step n (we will use alpha at t_n)
    tn = t(n+1); % note: using t_{n} or t_{n+1} is a choice; we pick t_n+1 for alpha evaluation so it's causal to the new step
    alpha_n = alphaFun(tn); 
    if ~(alpha_n > 0 && alpha_n <= 1)
        error('alphaFun must return values in (0,1]. Received %g at time %g', alpha_n, tn);
    end

    % precompute constants
    gammaTerm = gamma(alpha_n + 1);
    hAlphaOverGamma = h^alpha_n / gammaTerm;

    % convolutional coefficients for step n: b_m^{(n)} = ( (m+1)^{alpha_n} - m^{alpha_n} ), m = 0..(n-1)
    sumTerm = zeros(1, dim);
    for j = 0:(n-1)
        m = n - j - 1; % m >= 0
        bj = ( (m+1)^alpha_n - (m)^alpha_n );
        sumTerm = sumTerm + bj * F(j+1, :);
    end

    % explicit update: use initial Y(1,:) which corresponds to y0
    Y(n+1, :) = Y(1, :) + hAlphaOverGamma * sumTerm;

    % compute f at new point for future steps
    F(n+1, :) = rhsFun(t(n+1), Y(n+1, :)', params)';

end

end
