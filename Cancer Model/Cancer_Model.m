function dX = Cancer_Model(t, X, C)
% cancerSystem - RHS of integer-order cancer model

% X = [x; y; z]
% C = [C1 C2 C3 C4 C5 C6 C7 C8]

x = X(1); y = X(2); z = X(3);

C1 = C(1); C2 = C(2); C3 = C(3); C4 = C(4);
C5 = C(5); C6 = C(6); C7 = C(7); C8 = C(8);

dx = x*(1 - x) - C1*x*y - C2*x*z;
dy = C3*y*(1 - y) - C4*x*y;
dz = C5*(x*z/(x + C6)) - C7*x*z - C8*z;

dX = [dx; dy; dz];
end
