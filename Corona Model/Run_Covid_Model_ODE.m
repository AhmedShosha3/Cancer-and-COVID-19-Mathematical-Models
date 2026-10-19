clc; clear;

%% Initial conditions and parameters
y0 = [39000 0 1 5 0 0 0 0];
params = [2.55 1.56 7.65 0.25 0.580 0.001 0.94 0.27 0.5 3.5 1 0.3];

%% Time settings
tspan = [0 90];
dt = 0.01;
updateStep = 5; % plot every 5 steps

%% Solve using ode45
[t, y] = ode45(@(t,y) Covid_Model(t,y,params), tspan, y0);

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
    xlim([tspan(1) tspan(2)]);
    hLine{k} = plot(NaN, NaN, 'Color', compColors{k}, 'LineWidth', 2);
end

%% Real-time plotting loop
Npoints = length(t);

for n = 1:Npoints
    if mod(n, updateStep) == 0 || n == Npoints
        for k = 1:8
            if k ==1
            figure(figHandles(k));
            idx = 1:updateStep:n;
            set(hLine{k}, 'XData', t(idx), 'YData', y(idx,k));
            end
        end
        drawnow;
    end
end
