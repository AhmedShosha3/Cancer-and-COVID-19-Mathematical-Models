clear; clc; close all;

Y0 = [39000 ; 0 ;1 ; 5 ;0; 0; 0; 0]; % initial populations
alpha = 0.98; 
% alpha_vals = [0.8 0.9 0.98 1];    % fractional order
params = [2.55^alpha 1.56 7.65^alpha 0.25^alpha 0.580 0.001 0.94^alpha 0.27^alpha 0.5^alpha 3.5^alpha 1^alpha 0.3^alpha];
Tmax = 90;       % simulation time
h = 0.01;         % step size
T = 0:h:Tmax;     % time vector
N = length(T);

rhsFun = @(t,y,params) Covid_RHS_Wrapper(t,y,params);

% FRACTIONAL PECE METHOD
fprintf('Running Fractional PECE Solver...\n');
Y_PECE = zeros(N,8);
Y_PECE(1,:) = Y0';
Y_temp = Y0;


tic;
for k = 1:N-1
    tspan_step = [T(k) T(k+1)];
    [t_step, Y_step] = Frac_Caputo_PECE(alpha, T(k+1)-T(k), h, Y_temp, params, rhsFun);
    Y_PECE(k+1,:) = Y_step(end,:);
    Y_temp = Y_PECE(k+1,:)';
end
cpuPECE = toc;
fprintf('CPU Time (PECE): %.6f seconds\n', cpuPECE);

% FRACTIONAL EULER METHOD

fprintf('Running Fractional Euler Solver...\n');
Y_Euler = zeros(N,8);
Y_Euler(1,:) = Y0';
Y_temp = Y0;

tic;
for k = 1:N-1
    [~, Y_step] = Frac_Caputo_Euler(alpha, h, h, Y_temp, params, rhsFun);
    Y_Euler(k+1,:) = Y_step(end,:);
    Y_temp = Y_Euler(k+1,:)';
end
cpuEuler = toc;
fprintf('CPU Time (Euler): %.6f seconds\n', cpuEuler);

% TIME SERIES COMPARISON
figure('Color','white','Position',[100 100 1200 800]);

% ------- S -------
subplot(3,1,1)
plot(T, Y_PECE(:,1),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,1),'b--','LineWidth',2);
ylabel('S','FontSize',11);
legend('PECE','Euler','Location','best');

title(sprintf('Fractional Methods Comparison (\\alpha = %.2f)', alpha));

% ------- E -------
subplot(3,1,2)
plot(T, Y_PECE(:,2),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,2),'b--','LineWidth',2);
ylabel('E','FontSize',11);
legend('PECE','Euler','Location','best');

% ------- I -------
subplot(3,1,3)
plot(T, Y_PECE(:,3),'r','LineWidth',2); hold on;
plot(T, Y_Euler(:,3),'b--','LineWidth',2);
ylabel('I','FontSize',11);
xlabel('Time','FontSize',11);
legend('PECE','Euler','Location','best');


%figure('Color','white','Position',[100 100 1200 800]);
%% ------- P -------
%subplot(3,1,1)
%plot(T, Y_PECE(:,4),'r','LineWidth',2); hold on;
%plot(T, Y_Euler(:,4),'b--','LineWidth',2);
%ylabel('P','FontSize',11);
%xlabel('Time','FontSize',11);
%legend('PECE','Euler','Location','best');
%%
%% ------- A -------
%subplot(3,1,2)
%plot(T, Y_PECE(:,5),'r','LineWidth',2); hold on;
%plot(T, Y_Euler(:,5),'b--','LineWidth',2);
%ylabel('A','FontSize',11);
%xlabel('Time','FontSize',11);
%legend('PECE','Euler','Location','best');
%%
%% ------- H -------
%subplot(3,1,3)
%plot(T, Y_PECE(:,6),'r','LineWidth',2); hold on;
%plot(T, Y_Euler(:,6),'b--','LineWidth',2);
%ylabel('H','FontSize',11);
%xlabel('Time','FontSize',11);
%legend('PECE','Euler','Location','best');
%%
%figure('Color','white','Position',[100 100 1200 800]);
%% ------- R -------
%subplot(3,1,1)
%plot(T, Y_PECE(:,7),'r','LineWidth',2); hold on;
%plot(T, Y_Euler(:,7),'b--','LineWidth',2);
%ylabel('R','FontSize',11);
%xlabel('Time','FontSize',11);
%legend('PECE','Euler','Location','best');
%%
%% ------- F -------
%subplot(3,1,2)
%plot(T, Y_PECE(:,8),'r','LineWidth',2); hold on;
%plot(T, Y_Euler(:,8),'b--','LineWidth',2);
%ylabel('F','FontSize',11);
%xlabel('Time','FontSize',11);
%legend('PECE','Euler','Location','best');
%%

% Calculate absolute differences
diff_S = abs(Y_PECE(:,1) - Y_Euler(:,1));
diff_E = abs(Y_PECE(:,2) - Y_Euler(:,2));
diff_I = abs(Y_PECE(:,3) - Y_Euler(:,3));
diff_P = abs(Y_PECE(:,4) - Y_Euler(:,4));
diff_A = abs(Y_PECE(:,5) - Y_Euler(:,5));
diff_H = abs(Y_PECE(:,6) - Y_Euler(:,6));
diff_R = abs(Y_PECE(:,7) - Y_Euler(:,7));
diff_F = abs(Y_PECE(:,8) - Y_Euler(:,8));

% CPU TIME COMPARISON SUMMARY
fprintf('\n========================================\n');
fprintf('       CPU TIME COMPARISON\n');
fprintf('========================================\n');
fprintf('Fractional Order (alpha): %.2f\n', alpha);
fprintf('Time Step (h):            %.6f\n', h);
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
fprintf('Max |x_PECE - x_Euler|:   %.6e\n', max(diff_S));
fprintf('Max |y_PECE - y_Euler|:   %.6e\n', max(diff_E));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_I));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_P));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_A));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_H));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_R));
fprintf('Max |z_PECE - z_Euler|:   %.6e\n', max(diff_F));
fprintf('========================================\n');
