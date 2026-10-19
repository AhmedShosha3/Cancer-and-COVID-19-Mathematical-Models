clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5]; % PARAMETERS

Y0 = [0.1; 0.1; 0.1];% Initial conditions

% Solver settings
alpha = 0.9;    % fractional order
Tmax = 300;
h = 0.02;

% Preallocate time vector
T = 0:h:Tmax;
N = length(T);

% Preallocate solution
Y = zeros(N,3);
Y(1,:) = Y0';

% rhs wrapper
rhsFun = @(t,y,params) Cancer_RHS_Wrapper(t,y,params);

% 2D REAL-TIME PLOT
figure('Color','white');
h1 = plot(NaN, NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN, NaN,'g','LineWidth',2);
h3 = plot(NaN, NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title(sprintf('Fractional Caputo model, alpha = %.3g', alpha));

% 3D REAL-TIME PHASE PLOT
figure('Color','white');
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
hold on;
xlabel('x'); 
ylabel('y'); 
zlabel('z');
title(sprintf('Phase portrait (fractional), alpha = %.3g', alpha));

% Set axis limits
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);


view(-45, 20); % Change viewpoint

set(gca, 'Color', 'white'); % Set background to white

% Fractional PECE method
disp('Running Fractional PECE Solver...');
tic;  % Start CPU time

Y_temp = Y0;
updateStep = 30;

for k = 1:N-1
    tspan_step = [T(k) T(k+1)];
    [t_step, Y_step] = Frac_Caputo_PECE(alpha, T(k+1)-T(k), h, Y_temp, C, rhsFun);
    Y(k+1,:) = Y_step(end,:);   % store last value
    Y_temp = Y(k+1,:)';
    
    if mod(k, updateStep) == 0
        
        % 2D update
        set(h1,'XData',T(1:k),'YData',Y(1:k,1));
        set(h2,'XData',T(1:k),'YData',Y(1:k,2));
        set(h3,'XData',T(1:k),'YData',Y(1:k,3));
        
        % 3D update 
        set(h3d, 'XData', Y(1:k,1), 'YData', Y(1:k,2) , 'ZData', Y(1:k,3));
        
        drawnow;
    end
end

cpuFrac = toc;  % End CPU time
fprintf('CPU Time (fractional PECE) = %.6f seconds\n', cpuFrac);

drawnow;