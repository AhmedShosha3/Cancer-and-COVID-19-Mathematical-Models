clear; clc; close all;

%% Initial conditions
y0 = [39000 ; 0 ; 1 ; 5 ; 0 ; 0 ; 0 ; 0];

%% Fractional orders and colors
alpha_vals = [0.8 0.9 0.98 1];
colors = {'m','g','b','r'};
alpha_labels = {'\alpha = 0.8','\alpha = 0.9','\alpha = 0.98','\alpha = 1'};

%% Time settings
Tmax = 90;
h = 0.1;
updateStep = 5;   % update every 5 steps

%% Titles for compartments
titles = {'Susceptible','Exposed','Infected','Super-spreaders', ...
          'Asymptomatic','Hospitalized','Recovered','Fatalities'};

%% Create figures and empty plot handles
hLine = cell(8, length(alpha_vals));

for k = 1:8
    if k == 1
    figure(k);
    hold on;
    xlabel('Time (days)');
    ylabel('Population');
    title(titles{k});
    xlim([0 Tmax]);   % fix x-axis to prevent resizing
    end
end

%% Loop over alpha values
for i = 1:length(alpha_vals)
    
    alpha = alpha_vals(i);
    color = colors{i};

    %% Parameters (fractional powers applied)
    params = [2.55^alpha 1.56 7.65^alpha 0.25^alpha ...
              0.580 0.001 0.94^alpha 0.27^alpha ...
              0.5^alpha 3.5^alpha 1^alpha 0.3^alpha];

    %% RHS function handle
    rhsFun = @(t,y,params) Covid_RHS_Wrapper(t,y,params);

    %% Solve the fractional system
    [t, Yfrac] = Frac_Caputo_Euler(alpha, Tmax, h, y0', params, rhsFun);

    %% Create empty plot handles for this alpha
    for k = 1:8
        if k==1
        figure(k);
        hLine{k,i} = plot(NaN, NaN, color, 'LineWidth', 2);
        end
    end

    %% Real-time update loop (fast using updateStep)
    for n = 1:length(t)
            for k = 1:8
                if k ==1
                idx = 1:updateStep:n;   % decimate points for speed
                set(hLine{k,i},'XData', t(idx),'YData', Yfrac(idx,k));
                end
            end
            drawnow; 
    end
end

for k = 1:8
    if k == 1
    figure(1);
    legend([hLine{1,:}], alpha_labels, 'Location', 'best');
    end
end

