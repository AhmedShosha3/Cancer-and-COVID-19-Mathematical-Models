clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5]; % PARAMETERS
Y0 = [0.1; 0.1; 0.1]; 
Tmax = 300; 
h = 0.02; 

% define alpha(t) 
alphaFun = @(tt) 0.75 + 0.2*sin(0.05*tt); 
rhsFun = @(t,y,params) Cancer_RHS_Wrapper(t,y,params); 

% Preallocate time vector
t = 0:h:Tmax;
N = length(t);

% Preallocate solution
Y = zeros(N,3);
Y(1,:) = Y0';

% 2D REAL-TIME PLOT
figure('Color','white');
h1 = plot(NaN, NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN, NaN,'g','LineWidth',2);
h3 = plot(NaN, NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title('Variable-Order Fractional Model');

% 3D REAL-TIME PHASE PLOT
figure('Color','white');
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
hold on;
xlabel('X'); 
ylabel('Y'); 
zlabel('z');
title('Variable-Order Fractional Model');

% Set axis limits
xlim([0 1]);
ylim([0 1]);
zlim([0 1]);


view(-45, 20); % Change viewpoint

set(gca, 'Color', 'white'); % Set background to white

% Variable-order PECE method
disp('Running Variable-Order Fractional (PECE)...');
tic;

Y_temp = Y0;

updateStep = 30;

for k = 1:N-1
    tspan_step = [t(k) t(k+1)];
    [t_step, Y_step] = VO_Frac_PECE(alphaFun, t(k+1)-t(k), h, Y_temp, C, rhsFun); 
    Y(k+1,:) = Y_step(end,:);
    Y_temp = Y(k+1,:)';
    
    if mod(k, updateStep) == 0
        
        % 2D update
        set(h1,'XData',t(1:k),'YData',Y(1:k,1));
        set(h2,'XData',t(1:k),'YData',Y(1:k,2));
        set(h3,'XData',t(1:k),'YData',Y(1:k,3));
        
        % 3D update - with swapped x and y
        set(h3d, 'XData', Y(1:k,1), 'YData', Y(1:k,2), 'ZData', Y(1:k,3));
        
        drawnow;
    end
end

cpuTime = toc;
fprintf('CPU Time (variable-order PECE) = %.6f seconds\n', cpuTime);

drawnow;

% Plot alpha(t)
figure;
plot(t, arrayfun(alphaFun,t), 'k','LineWidth',2);
xlabel('t'); ylabel('\alpha(t)');
title('Variable-order function \alpha(t)');
drawnow;