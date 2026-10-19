clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5];    % model parameters
Y0 = [0.1; 0.1; 0.1];                 % initial populations
alpha = 0.9;      % fractional order
Tmax = 300;       % simulation time
h = 0.02;         % step size
T = 0:h:Tmax;     % time vector
N = length(T);
Y = zeros(N,3);
Y(1,:) = Y0';
rhsFun = @(t,y,params) Cancer_RHS_Wrapper(t,y,params);

% 2D REAL-TIME PLOT
figure('Color','white');
h1 = plot(NaN, NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN, NaN,'g','LineWidth',2);
h3 = plot(NaN, NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title(sprintf('Fractional Euler Method (\\alpha = %.2f)', alpha));

% 3D REAL-TIME PHASE PLOT
figure('Color','white');
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
hold on;
xlabel('x'); 
ylabel('y'); 
zlabel('z');
title(sprintf('Fractional Euler Method (\\alpha = %.2f)', alpha));

% Set axis limits
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);

view(-45, 20); % Change viewpoint

set(gca, 'Color', 'white'); % Set background to white

% FRACTIONAL EULER METHOD (Caputo)
disp('Running Fractional Euler Solver...');
tic;

Ytemp = Y0;
updateStep = 30;

for k = 1:N-1
    % one fractional Euler step from T(k) to T(k+1)
    [~, Y_step] = Frac_Caputo_Euler(alpha, h, h, Ytemp, C, rhsFun);
    Y(k+1,:) = Y_step(end,:);
    Ytemp = Y(k+1,:)';
    
    if mod(k, updateStep) == 0
        
        % 2D update
        set(h1,'XData',T(1:k),'YData',Y(1:k,1));
        set(h2,'XData',T(1:k),'YData',Y(1:k,2));
        set(h3,'XData',T(1:k),'YData',Y(1:k,3));
        
        % 3D update - with swapped x and y
        set(h3d, 'XData', Y(1:k,1), 'YData', Y(1:k,2) , 'ZData', Y(1:k,3));
        
        drawnow;
    end
end

cpu_time = toc;
fprintf('CPU Time (Euler Fractional) = %.6f seconds\n', cpu_time);

drawnow;