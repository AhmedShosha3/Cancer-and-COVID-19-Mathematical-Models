clear; clc; close all;

% Initial conditions
y0 = [43999; 0; 1; 5; 0; 0; 0; 0];

% Solver settings
alpha_vals = [0.8, 0.9, 0.98, 1];    % Fractional orders to be tested
Tmax = 90;                            % Maximum time
h = 0.1;                              % Time step

% State variable names and labels
state_names = {'S', 'E', 'I', 'P', 'A', 'H', 'R', 'D'};
state_full_names = {'Susceptible', 'Exposed', 'Infected', 'Super-spreaders', ...
                    'Asymptomatic', 'Hospitalized', 'Recovered', 'Fatalities'};

% Line colors for each alpha (matching your image)
line_colors = {'b', 'r', 'g', 'm'};  % Blue, Red, Green, Magenta

% Create single figure with 4x2 subplots
figure('Position', [50, 50, 900, 1200]);

% Solve for all alpha values first
solutions = cell(length(alpha_vals), 1);
for i = 1:length(alpha_vals)
    alpha = alpha_vals(i);
    
    % Parameters scaled by alpha
    params = [2.55^alpha, 1.56, 7.65^alpha, 0.25^alpha, 0.580, ...
              0.001, 0.94^alpha, 0.27^alpha, 0.5^alpha, 3.5^alpha, ...
              1^alpha, 0.3^alpha];
    
    % RHS wrapper function
    rhsFun = @(t, y, params) Covid_RHS_Wrapper(t, y, params);
    
    % Solve the fractional system
    [t, Yfrac] = Frac_Caputo_PECE(alpha, Tmax, h, y0, params, rhsFun);
    
    solutions{i}.t = t;
    solutions{i}.Y = Yfrac;
end

% Plot each state variable in its own subplot
for state = 1:8
    subplot(4, 2, state);
    hold on;
    
    % Create legend entries
    legend_entries = cell(1, length(alpha_vals));
    
    % Plot each alpha value
    for i = 1:length(alpha_vals)
        plot(solutions{i}.t, solutions{i}.Y(:, state), ...
             'Color', line_colors{i}, 'LineWidth', 2);
        legend_entries{i} = ['\alpha=' num2str(alpha_vals(i))];
    end
    
    % Formatting
    xlabel('t', 'FontSize', 10);
    ylabel(state_names{state}, 'FontSize', 10);
    title(state_full_names{state}, 'FontSize', 11, 'FontWeight', 'normal');
    
    % Add legend to each subplot
    legend(legend_entries, 'Location', 'Best', 'FontSize', 8);
    grid on;
    box on;
    
    % Set axis limits if needed
    xlim([0 Tmax]);
    ylim_current = ylim;
    if ylim_current(1) < 0
        ylim([0 ylim_current(2)]);
    end
    
    hold off;
end

% Adjust subplot spacing
set(gcf, 'PaperPositionMode', 'auto');