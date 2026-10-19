function [t, Y] = VO_Frac_PECE(alphaFun, T, h, Y0, params, rhsFun)
% voFracPECE - variable-order fractional PECE solver
% alphaFun : function handle @(t) returning alpha in (0,1]
% other arguments like fracCaputoPECE

N = floor(T/h);
t = (0:h:N*h)';
dim = numel(Y0);
Y = zeros(N+1, dim);
Y(1,:) = Y0(:)';

F = zeros(N+1, dim);
F(1,:) = rhsFun(t(1), Y(1,:)', params)';

for n = 1:N
    a = alphaFun(t(n+1));  % use alpha at the new time (can also use t(n))
    if a <= 0 || a > 1
        error('alphaFun must return values in (0,1].');
    end
    gammaTerm = gamma(a + 1);
    sumTerm = zeros(1,dim);
    for j = 0:n-1
        bj = ( (n-j)^(a) - (n-j-1)^(a) );
        sumTerm = sumTerm + bj * F(j+1,:);
    end
    Yp = Y(1,:) + (h^a / gammaTerm) * sumTerm;
    fp = rhsFun(t(n+1), Yp', params)';
    sumTerm2 = sumTerm + 1 * fp;
    Yn = Y(1,:) + (h^a / gammaTerm) * sumTerm2;

    Y(n+1,:) = Yn;
    F(n+1,:) = rhsFun(t(n+1), Yn', params)';
end
end
