clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5]; % parameters
X0 = [0.1; 0.1; 0.1]; % initial conditions
tspan = [0 300]; % time interval
dt = 0.01; % time step

% Preallocate time vector
T = tspan(1):dt:tspan(2);
N = length(T);

% Preallocate solution
X = zeros(N,3);
X(1,:) = X0';

% 2D REAL-TIME PLOT
figure('Color','white');
h1 = plot(NaN,NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN,NaN,'g','LineWidth',2);
h3 = plot(NaN,NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title('Integer-order Cancer Model (Euler)');

% 3D REAL-TIME PHASE PLOT
figure('Color','white');
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
hold on;
xlabel('x'); 
ylabel('y'); 
zlabel('z');
title('Integer-order Cancer Model (Euler)');

% Set axis limits
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

view(-45, 20); % Change viewpoint

set(gca, 'Color', 'white'); % Set background to white


% EULER SOLVER
disp('Running Integer Euler Solver...');
tic;

X_temp = X0;

updateStep = 50;

for k = 1:N-1
    
    % Use Euler method function for single step
    tspan_step = [T(k) T(k+1)];
    [~, X_step] = Cancer_Model_Euler(X_temp, C, tspan_step, dt);
    
    X(k+1,:) = X_step(end,:);
    X_temp = X(k+1,:).';
    
    if mod(k, updateStep) == 0
        
        % 2D update
        set(h1, 'XData', T(1:k), 'YData', X(1:k,1));
        set(h2, 'XData', T(1:k), 'YData', X(1:k,2));
        set(h3, 'XData', T(1:k), 'YData', X(1:k,3));
        
        % 3D update 
        set(h3d, 'XData', X(1:k,1), 'YData', X(1:k,2), 'ZData', X(1:k,3));
        
        drawnow;
    end
    
end

cpuEuler = toc;
fprintf('CPU Time (Euler function) = %.6f seconds\n', cpuEuler);

drawnow;