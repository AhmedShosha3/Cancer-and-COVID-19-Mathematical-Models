function [t, Y] = Frac_Caputo_Euler(alpha, T, h, Y0, params, rhsFun)
% fracCaputoEuler  Explicit Euler-type solver for Caputo fractional ODE
%
% Solves D^alpha y(t) = f(t,y), 0 < alpha <= 1, with y(0) = Y0
% using an explicit fractional-Euler (convolution) scheme:
%
% % %   y_n = y_0 + (h^alpha / Gamma(alpha+1)) * sum_{j=0}^{n-1} b_{n-j-1} f(t_j,y_j)
%   where b_m = (m+1)^alpha - m^alpha
%
% Inputs:
%   alpha   - fractional order (0 < alpha <= 1)
%   T       - final time
%   h       - time step
%   Y0      - initial state vector (column or row)
%   params  - additional parameters passed to rhsFun (can be [])
%   rhsFun  - function handle f = rhsFun(t, y, params), returns column vector
%
% Outputs:
%   t       - (N+1)x1 time vector
%   Y       - (N+1) x dim matrix of solution (rows correspond to times)
%
% Example usage:
%   rhs = @(t,y,p) -p.k * y;              % scalar linear example
%   [t,Y] = fracCaputoEuler(0.8, 1, 0.01, 1, struct('k',1), rhs);
%

% --- input checks -------------------------------------------------------
if ~(alpha > 0 && alpha <= 1)
    error('alpha must be in (0,1].');
end
if h <= 0 || T <= 0
    error('h and T must be positive.');
end

% ensure Y0 is a column for internal operations
Y0 = Y0(:);
dim = numel(Y0);

N = floor(T / h);
t = (0:N)' * h;         % column vector (N+1 x 1)
Y = zeros(N+1, dim);    % solution storage (rows are time steps)
Y(1, :) = Y0.';         % initial condition (row)

% Precompute constant factor
gammaTerm = gamma(alpha + 1);
hAlphaOverGamma = h^alpha / gammaTerm;

% store f evaluations at previous grid points for efficiency
F = zeros(N+1, dim);
F(1, :) = rhsFun(t(1), Y(1,:)', params)';  % f(t0, y0) as row

% time-stepping loop (explicit convolution form)
for n = 1:N
    % compute convolution sum: sum_{j=0}^{n-1} b_{n-j-1} * F(j+1,:)
    sumTerm = zeros(1, dim);
    % we compute for j = 0 .. n-1, index in F is j+1
    for j = 0:(n-1)
        m = n - j - 1;  % m >= 0
        bj = ( (m+1)^alpha - (m)^alpha );   % b_m
        sumTerm = sumTerm + bj * F(j+1, :);
    end

    % explicit update (uses only previous f-values)
    Y(n+1, :) = Y(1, :) + hAlphaOverGamma * sumTerm;

    % compute f at new point for future steps
    F(n+1, :) = rhsFun(t(n+1), Y(n+1, :)', params)';
end

end
