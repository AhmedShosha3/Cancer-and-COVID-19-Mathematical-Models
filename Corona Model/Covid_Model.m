function dydt = Covid_Model(t, y, params)

%x = X(1); y = X(2); z = X(3);
%C1 = C(1); C2 = C(2); C3 = C(3); C4 = C(4);
%C5 = C(5); C6 = C(6); C7 = C(7); C8 = C(8);
%dx = x*(1 - x) - C1*x*y - C2*x*z;
%dy = C3*y*(1 - y) - C4*x*y;
%dz = C5*(x*z/(x + C6)) - C7*x*z - C8*z;
%dX = [dx; dy; dz];

% State variables
S = y(1);
E = y(2);
I = y(3);
P = y(4);
A = y(5);
H = y(6);
R = y(7);
F = y(8);

N = S + E + I + P + A + H + R + F;

% Parameters
beta   = params(1);
beta_p = params(2);
l      = params(3);
K      = params(4);
rho1   = params(5);
rho2   = params(6);
ga     = params(7);
gi     = params(8);
gr     = params(9);
di     = params(10);
dp     = params(11);
dh     = params(12);


dS = -beta*I*S/N - l*beta*H*S/N - beta_p*P*S/N;
dE =  beta*I*S/N + l*beta*H*S/N + beta_p*P*S/N - K*E;
dI =  K*rho1*E - (ga + gi + di)*I;
dP =  K*rho2*E - (ga + gi + dp)*P;
dA =  K*(1 - rho1 - rho2)*E;
dH =  ga*(I + P) - (gr + dh)*H;
dR =  gi*(I + P) + gr*H;
dF =  di*I + dp*P + dh*H;

dydt = [dS; dE; dI; dP; dA; dH; dR; dF];



end
