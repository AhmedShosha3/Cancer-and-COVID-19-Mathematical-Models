function [T, X] = Cancer_Model_Euler(X0, C, tspan, dt)
% cancerSystemEuler - Solves the cancer system using Euler's method
%
% Inputs:
%   X0    - Initial conditions [x0, y0, z0]
%   C     - Parameter vector [C1 C2 C3 C4 C5 C6 C7 C8]
%   tspan - [t0, tf] time interval
%   dt    - Time step
%
% Outputs:
%   T - Time vector
%   X - Solution matrix [x, y, z] at each time step

    % Time vector
    T = tspan(1):dt:tspan(2);
    N = length(T);

    % Initialize solution
    X = zeros(N, 3);
    X(1, :) = X0;

    % Euler integration loop
    for k = 1:N-1
        x = X(k, 1);
        y = X(k, 2);
        z = X(k, 3);

        % Parameters
        C1 = C(1); C2 = C(2); C3 = C(3); C4 = C(4);
        C5 = C(5); C6 = C(6); C7 = C(7); C8 = C(8);

        % Compute derivatives
        dx = x*(1 - x) - C1*x*y - C2*x*z;
        dy = C3*y*(1 - y) - C4*x*y;
        dz = C5*(x*z/(x + C6)) - C7*x*z - C8*z;

        % Euler update
        X(k+1, 1) = x + dt*dx;
        X(k+1, 2) = y + dt*dy;
        X(k+1, 3) = z + dt*dz;
    end
end
