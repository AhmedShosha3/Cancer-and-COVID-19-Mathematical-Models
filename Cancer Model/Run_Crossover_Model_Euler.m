clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5]; % Parameters

Y0 = [0.1; 0.1; 0.1]; % Initial conditions

% Time parameters
t1 = 100; t2 = 200; T = 300; h = 0.02;
alpha_const = 0.85;
alphaFun = @(tt) 0.7 + 0.25*sin(0.03*tt);

% Measure CPU time
disp('Running Crossover Model (Euler)...');
tic;
[t_all, Y_all] = Crossover_Model_Euler(t1, t2, T, h, Y0, C, alpha_const, alphaFun);
cpuTime = toc;
fprintf('CPU Time (Crossover Model) = %.6f seconds\n', cpuTime);

% 2D REAL-TIME PLOT
figure('Color','white');
h1 = plot(NaN,NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN,NaN,'g','LineWidth',2);
h3 = plot(NaN,NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title('Crossover model Euler (int to frac to varfrac)');

% 3D REAL-TIME PHASE PLOT
figure('Color','white');
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
hold on;
xlabel('x'); ylabel('y'); zlabel('z');
title(' Crossover model Euler ');

% Set axis limits
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

view(-45, 20); % Change viewpoint

set(gca, 'Color', 'white'); % Set background to white

% REAL-TIME UPDATES
updateStep = 25;  

for k = 1:length(t_all)
    
    % 2D update
    set(h1,'XData',t_all(1:k),'YData',Y_all(1:k,1));
    set(h2,'XData',t_all(1:k),'YData',Y_all(1:k,2));
    set(h3,'XData',t_all(1:k),'YData',Y_all(1:k,3));
    
    % 3D update
    if mod(k, updateStep) == 0
        set(h3d,'XData',Y_all(1:k,1),'YData',Y_all(1:k,2),'ZData',Y_all(1:k,3));
        drawnow;
    end
end

drawnow;