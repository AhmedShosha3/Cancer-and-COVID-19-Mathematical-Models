function [t, y] = Frac_Caputo_PECE(alpha, Tmax, h, y0, params, rhsFun)
% fracCaputoPECE  Solve system D^alpha y = f(t,y) (Caputo) using PECE
%
% alpha  : fractional order (0<alpha<=1)
% T      : final time
% h      : time step
% Y0     : initial vector (column) y(0)
% params : additional params passed to rhsFun
% rhsFun : function handle f = rhsFun(t, y, params) returning column vector
%
% Returns: t (column), Y (N+1 x dim)

if alpha <= 0 || alpha > 1
    error('alpha must be in (0,1].');
end

N = floor(Tmax/h);
t = (0:h:N*h)';    % column
dim = numel(y0);

y = zeros(N+1, dim);
y(1,:) = y0(:)';

% Precompute Gamma factor
gammaTerm = gamma(alpha + 1);

% Precompute coefficient kernel b_{n-j-1} = (n-j)^{alpha} - (n-j-1)^{alpha}
% We'll compute on the fly inside loop (efficient enough for moderate N)

% store f evaluations at grid points
F = zeros(N+1, dim);
F(1,:) = rhsFun(t(1), y(1,:)', params)';

for n = 1:N
    % predictor: explicit formula
    sumTerm = zeros(1,dim);
    for j = 0:n-1
        bj = ( (n-j)^(alpha) - (n-j-1)^(alpha) );
        sumTerm = sumTerm + bj * F(j+1,:); % F index shift: j=0 -> 1
    end
    Yp = y(1,:) + (h^alpha / gammaTerm) * sumTerm; % predictor row

    % predictor f
    fp = rhsFun(t(n+1), Yp', params)'; % row

    % corrector: include fp with weight 1 (PECE style)
    % some PECE variants put different weight; this version is a standard simple correction
    sumTerm2 = sumTerm + 1 * fp;
    yn = y(1,:) + (h^alpha / gammaTerm) * sumTerm2;

    % accept Yn, compute f_n
    y(n+1,:) = yn;
    F(n+1,:) = rhsFun(t(n+1), yn', params)';
end

% transpose to column-based output already set
end
