clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5]; % Parameters

X0 = [0.1; 0.1; 0.1]; % Initial conditions

% Time span
tspan = [0 100];
dt = 0.01;   % Euler time step

% Euler Method
tic;
[T_euler, X_euler] = Cancer_Model_Euler(X0, C, tspan, dt);
CPU_Euler = toc;
fprintf('CPU time (Euler): %.6f seconds\n', CPU_Euler);

% ODE Method

tic;
[t_ode45, X_ode45] = ode45(@(t,X) Cancer_Model(t,X,C), tspan, X0);
CPU_ODE = toc;
fprintf('CPU time (ODE): %.6f seconds\n', CPU_ODE);


% Plot Comparison
figure;

% ------- Tumor x -------
subplot(3,1,1)
plot(T_euler, X_euler(:,1),'r','LineWidth',1.4); hold on;
plot(t_ode45, X_ode45(:,1),'k.','LineWidth',1.4);
ylabel('Tumor x');
legend('Euler','ODE','Location','best');
grid on;

% ------- Healthy y -------
subplot(3,1,2)
plot(T_euler, X_euler(:,2),'g','LineWidth',1.4); hold on;
plot(t_ode45, X_ode45(:,2),'k.','LineWidth',1.4);
ylabel('Healthy y');
legend('Euler','ODE','Location','best');
grid on;

% ------- Immune z -------
subplot(3,1,3)
plot(T_euler, X_euler(:,3),'b','LineWidth',1.4); hold on;
plot(t_ode45, X_ode45(:,3),'k.','LineWidth',1.4);
ylabel('Immune z'); xlabel('Time');
legend('Euler','ODE','Location','best');
grid on;

% ------- Title for entire figure -------
annotation('textbox', [0, 0.95, 1, 0.04], ...
    'String', 'Comparison of Euler Method vs ODE (Reference)', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 12);

%   Print CPU Comparison
fprintf('\n=== CPU Time Comparison ===\n');
fprintf('Euler (dt = %.4f):%.6f seconds\n', dt, CPU_Euler);
fprintf('ode45:              %.6f seconds\n', CPU_ODE);
fprintf('Speed Ratio (Euler / ODE) = %.3f\n', CPU_Euler / CPU_ODE);
