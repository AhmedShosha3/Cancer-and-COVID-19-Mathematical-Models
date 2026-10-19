clc; clear all;

%% Initial conditions and parameters
y0 = [39000 0 1 5 0 0 0 0];
params = [2.55 1.56 7.65 0.25 0.580 0.001 0.94 0.27 0.5 3.5 1 0.3];

%% Time settings
Tmax = 90;
dt = 0.01;
updateStep = 5; % plot every 5 steps

%% Solve the model using Euler method
[T, y] = Covid_Model_Euler(y0, params, Tmax, dt);

%% Compartment titles and colors
compTitles = {'Susceptible','Exposed','Infected','Super-spreaders',...
              'Asymptomatic','Hospitalized','Recovered','Fatalities'};
compColors = {'b','r','y','m','r','g','k','b'};

%% Create figures and initialize plot handles
hLine = cell(8,1);
figHandles = zeros(8,1);

for k = 1:8
    figHandles(k) = figure;
    hold on;
    xlabel('Time (days)');
    ylabel('Population');
    title(compTitles{k});
    xlim([0 Tmax]);
    hLine{k} = plot(NaN, NaN, 'Color', compColors{k}, 'LineWidth', 2);
end

%% Real-time plotting loop
N = length(T);

for n = 1:N
    if mod(n, updateStep) == 0 || n == N
        for k = 1:8
            if k ==1
            figure(figHandles(k));
            idx = 1:updateStep:n;
            set(hLine{k}, 'XData', T(idx), 'YData', y(idx,k));
            end
        end
        drawnow; 
    end
end
