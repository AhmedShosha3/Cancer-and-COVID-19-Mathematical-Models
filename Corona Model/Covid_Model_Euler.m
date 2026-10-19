function [T,Y] = Covid_Model_Euler(y, params,t,dt)

%x = X(1); y = X(2); z = X(3);
%C1 = C(1); C2 = C(2); C3 = C(3); C4 = C(4);
%C5 = C(5); C6 = C(6); C7 = C(7); C8 = C(8);
%dx = x*(1 - x) - C1*x*y - C2*x*z;
%dy = C3*y*(1 - y) - C4*x*y;
%dz = C5*(x*z/(x + C6)) - C7*x*z - C8*z;
%dX = [dx; dy; dz];
T = 0:dt:t;
N = length(T);

% Initialize solution
Y = zeros(N, 8);
Y(1, :) = y;

% Euler integration loop
for k = 1:N-1

% State variables
S = Y(k,1);
E = Y(k,2);
I = Y(k,3);
P = Y(k,4);
A = Y(k,5);
H = Y(k,6);
R = Y(k,7);
F = Y(k,8);

Npop = S + E + I + P + A + H + R + F;

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


dS = -beta*I*S/Npop - l*beta*H*S/Npop - beta_p*P*S/Npop;
dE =  beta*I*S/Npop + l*beta*H*S/Npop + beta_p*P*S/Npop - K*E;
dI =  K*rho1*E - (ga + gi + di)*I;
dP =  K*rho2*E - (ga + gi + dp)*P;
dA =  K*(1 - rho1 - rho2)*E;
dH =  ga*(I + P) - (gr + dh)*H;
dR =  gi*(I + P) + gr*H;
dF =  di*I + dp*P + dh*H;

Y(k+1,1)= S +dt*dS;
Y(k+1,2)= E +dt*dE;
Y(k+1,3)= I +dt*dI;
Y(k+1,4)= P +dt*dP;
Y(k+1,5)= A +dt*dA;
Y(k+1,6)= H +dt*dH;
Y(k+1,7)= R +dt*dR;
Y(k+1,8)= F +dt*dF;
end




end