clc; clear;

% Initial conditions
y0 = [39000 0 1 5 0 0 0 0];
params = [2.55 1.56 7.65 0.25 0.580 0.001 0.94 0.27 0.5 3.5 1 0.3];

% Time settings
t_end = 90;
dt = 0.01;
tspan = [0 t_end];

%% ---------------- Euler Method ----------------
tic;
[T_euler, y_euler] = Covid_Model_Euler(y0, params, t_end, dt);
cpu_euler = toc;

%% ---------------- ODE45 Method ----------------
tic;
[t_ode45, y_ode45] = ode45(@(t,y) Covid_Model(t,y,params), tspan, y0);
cpu_ode45 = toc;

%% ---------------- Display Results ----------------
fprintf('CPU Time (Euler Method)  : %.6f seconds\n', cpu_euler);
fprintf('CPU Time (ODE45 Method) : %.6f seconds\n', cpu_ode45);

%% ---------------- Bar Plot Comparison ----------------
figure;
bar([cpu_euler cpu_ode45]);
set(gca,'XTickLabel',{'Euler','ODE45'});
ylabel('CPU Time (seconds)');
title('CPU Time Comparison Between Euler and ODE45');
grid on;
