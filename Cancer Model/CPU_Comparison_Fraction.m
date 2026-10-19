clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5];    % model parameters
Y0 = [0.1; 0.1; 0.1];                 % initial populations
alpha = 0.9;      % fractional order
Tmax = 300;       % simulation time
h = 0.1;         % step size
T = 0:h:Tmax;     % time vector
N = length(T);

rhsFun = @(t,y,params) Cancer_RHS_Wrapper(t,y,params);

% FRACTIONAL PECE METHOD
fprintf('Running Fractional PECE Solver...\n');
Y_PECE = zeros(N,3);
Y_PECE(1,:) = Y0';
Y_temp = Y0;

tic;
for k = 1:N-1
    tspan_step = [T(k) T(k+1)];
    [t_step, Y_step] = Frac_Caputo_PECE(alpha, T(k+1)-T(k), h, Y_temp, C, rhsFun);
    Y_PECE(k+1,:) = Y_step(end,:);
    Y_temp = Y_PECE(k+1,:)';
end
cpuPECE = toc;
fprintf('CPU Time (PECE): %.6f seconds\n', cpuPECE);

% FRACTIONAL EULER METHOD

fprintf('Running Fractional Euler Solver...\n');
Y_Euler = zeros(N,3);
Y_Euler(1,:) = Y0';
Y_temp = Y0;

tic;
for k = 1:N-1
    [~, Y_step] = Frac_Caputo_Euler(alpha, h, h, Y_temp, C, rhsFun);
    Y_Euler(k+1,:) = Y_step(end,:);
    Y_temp = Y_Euler(k+1,:)';
end
cpuEuler = toc;
fprintf('CPU Time (Euler): %.6f seconds\n', cpuEuler);

% TIME SERIES COMPARISON
figure('Color','white','Position',[100 100 1200 800]);

% ------- Tumor x -------
subplot(3,1,1)
plot(T, Y_PECE(:,1),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,1),'b--','LineWidth',2);
ylabel('Tumor Population (x)','FontSize',11);
legend('PECE','Euler','Location','best');

title(sprintf('Fractional Methods Comparison (\\alpha = %.2f)', alpha));

% ------- Healthy y -------
subplot(3,1,2)
plot(T, Y_PECE(:,2),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,2),'b--','LineWidth',2);
ylabel('Healthy Cells (y)','FontSize',11);
legend('PECE','Euler','Location','best');

% ------- Immune z -------
subplot(3,1,3)
plot(T, Y_PECE(:,3),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,3),'b--','LineWidth',2);
ylabel('Immune Cells (z)','FontSize',11);
xlabel('Time','FontSize',11);
legend('PECE','Euler','Location','best');

% 3D PHASE PORTRAIT COMPARISON
figure('Color','white','Position',[150 150 1200 500]);

% PECE Phase Portrait
subplot(1,2,1)
plot3(Y_PECE(:,1), Y_PECE(:,2), Y_PECE(:,3),'r','LineWidth',2);
hold on;
plot3(Y_PECE(1,1), Y_PECE(1,2), Y_PECE(1,3),'go','MarkerSize',10,'LineWidth',2);
plot3(Y_PECE(end,1), Y_PECE(end,2), Y_PECE(end,3),'rs','MarkerSize',10,'LineWidth',2);
xlabel('Tumor (x)','FontSize',11); 
ylabel('Healthy (y)','FontSize',11); 
zlabel('Immune (z)','FontSize',11);
title(sprintf('PECE Method (CPU: %.4f s)', cpuPECE));
xlim([0 1]); ylim([0 1]); zlim([0 1]);
view(-45, 20);
legend('Trajectory','Start','End','Location','best');

% Euler Phase Portrait
subplot(1,2,2)
plot3(Y_Euler(:,1), Y_Euler(:,2), Y_Euler(:,3),'b','LineWidth',2);
hold on;
plot3(Y_Euler(1,1), Y_Euler(1,2), Y_Euler(1,3),'go','MarkerSize',10,'LineWidth',2);
plot3(Y_Euler(end,1), Y_Euler(end,2), Y_Euler(end,3),'rs','MarkerSize',10,'LineWidth',2);
xlabel('Tumor (x)','FontSize',11); 
ylabel('Healthy (y)','FontSize',11); 
zlabel('Immune (z)','FontSize',11);
title(sprintf('Euler Method (CPU: %.4f s)', cpuEuler));
xlim([0 1]); ylim([0 1]); zlim([0 1]);
view(-45, 20);
legend('Trajectory','Start','End','Location','best');

% DIFFERENCE ANALYSIS

figure('Color','white','Position',[200 200 1000 600]);

% Calculate absolute differences
diff_x = abs(Y_PECE(:,1) - Y_Euler(:,1));
diff_y = abs(Y_PECE(:,2) - Y_Euler(:,2));
diff_z = abs(Y_PECE(:,3) - Y_Euler(:,3));

subplot(3,1,1)
plot(T, diff_x,'m','LineWidth',2);
ylabel('|x_{PECE} - x_{Euler}|','FontSize',11);
title('Absolute Difference Between Methods');

subplot(3,1,2)
plot(T, diff_y,'m','LineWidth',2);
ylabel('|y_{PECE} - y_{Euler}|','FontSize',11);

subplot(3,1,3)
plot(T, diff_z,'m','LineWidth',2);
ylabel('|z_{PECE} - z_{Euler}|','FontSize',11);
xlabel('Time','FontSize',11);

% CPU TIME COMPARISON SUMMARY
fprintf('\n========================================\n');
fprintf('       CPU TIME COMPARISON\n');
fprintf('========================================\n');
fprintf('Fractional Order (alpha): %.2f\n', alpha);
fprintf('Time Step (h):            %.4f\n', h);
fprintf('Simulation Time (Tmax):   %.1f\n', Tmax);
fprintf('Number of Steps:          %d\n', N);
fprintf('----------------------------------------\n');
fprintf('PECE Method:              %.6f seconds\n', cpuPECE);
fprintf('Euler Method:             %.6f seconds\n', cpuEuler);
fprintf('----------------------------------------\n');
fprintf('Speed Ratio (PECE/Euler): %.3f\n', cpuPECE / cpuEuler);
fprintf('========================================\n');

% Maximum absolute errors
fprintf('\n========================================\n');
fprintf('     MAXIMUM ABSOLUTE DIFFERENCES\n');
fprintf('========================================\n');
fprintf('Max |x_PECE - x_Euler|:   %.6e\n', max(diff_x));
fprintf('Max |y_PECE - y_Euler|:   %.6e\n', max(diff_y));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_z));
fprintf('========================================\n');