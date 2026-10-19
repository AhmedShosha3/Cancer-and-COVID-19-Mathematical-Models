clear; clc; close all;

% Initial conditions
y0 = [10999940 30 15 5 10 0 0 0];

% Real data time points - SHIFTED BY 8 DAYS TO THE RIGHT
time_shift = 17;
treal = [0 2 4 6 8 10 12 14 15 16 17 18 19 20 21 22 23 24 25 26 ...
     27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 ...
     45 46 48 50 52 54 56 58 60 62 64 65] + time_shift;

% Real data values
RealData = [5 8 15 30 60 120 300 700 1700 1500 1800 1600 ...
            2000 2200 2600 2900 3300 3800 3900 3700 3200 ...
            3400 2800 3000 2500 2300 2000 1800 1600 1400 ...
            2000 1800 1900 500 900 600 500 450 420 400 ...
            450 200 150 120 90 70 50 40 30 20];

% Solver settings
alpha_vals = [0.1 0.975];    % Fractional orders to be tested
Tmax = 90;                   % Maximum time (extended to match real data)
h = 0.1;                     % Time step

% Preallocate time vector
T = 0:h:Tmax;
N = length(T);

% Preallocate solution
y = zeros(N,8);
y(1,:) = y0';

% Loop over each alpha value
for i = 1:length(alpha_vals)
    alpha = alpha_vals(i);  % Current alpha value
    params = [ ...
    2.55  ... % beta   (symptomatic transmission)
    1.56  ... % beta_p (super-spreader transmission)
    0.58  ... % l      (hospital contact factor)
    0.25  ... % K      (incubation transition)
    0.58  ... % rho1   (I -> R fraction)
    0.27  ... % rho2   (I -> H fraction)
    1.00  ... % ga     (recovery rate of A)
    0.50  ... % gi     (recovery rate of I)
    0.30  ... % gr     (recovery rate from H)
    0.001 ... % di     (death rate of I)
    0.002 ... % dp     (death rate of P)
    0.005 ... % dh     (death rate of H)
    ];
    
    % rhs wrapper
    rhsFun = @(t,y,params) Covid_RHS_Wrapper(t,y,params);
    
    % Solve the system with the current alpha
    [t, Yfrac] = Frac_Caputo_PECE(alpha, Tmax, h, y(1,:), params, rhsFun);
    
  
    
end

% Plot the model results and real data
t_frac = linspace(0, Tmax, length(t));
y3 = Yfrac(:,3);  % Infected compartment

% Calculate scaling factor to match real data to model scale
% Find the ratio between model peak and real data peak
scale_factor = max(y3) / max(RealData);
RealData_scaled = RealData * scale_factor;

figure
plot(t_frac, y3, 'LineWidth', 2);hold on;
plot(treal, RealData_scaled, 'o', 'MarkerSize', 8, 'LineWidth', 2)
xlabel('t (days)')
ylabel('I(t) - Infected Population')
legend('Fractional model', sprintf('Real Data (scaled by %.2f, shifted +%d days)', scale_factor, time_shift), 'Location', 'best')
grid on
title('COVID-19 Infected Population: Model vs Real Data')

% Display dimensions for verification
fprintf('Length of RealData: %d\n', length(RealData));
fprintf('Length of treal: %d\n', length(treal));
fprintf('Length of model time: %d\n', length(t_frac));
fprintf('Model time range: %.1f to %.1f days\n', t_frac(1), t_frac(end));
fprintf('Time shift applied: %d days\n', time_shift);
fprintf('Scaling factor applied: %.2f\n', scale_factor);
fprintf('Real data max: %d, Model max: %.0f\n', max(RealData), max(y3));