clear; clc; close all;

%% Initial conditions
y0 = [39000 ; 0 ; 1 ; 5 ; 0; 0; 0; 0];

%% Solver settings
alpha_vals = [0.8 0.9 0.98 1];    % Fractional orders
Tmax = 90;                        % Maximum time
h = 0.1;                          % Time step
updateStep = 30;                   % Update every 5 steps

%% Compartment titles
compTitles = {'Susceptible','Exposed','Infected','Super-spreaders', ...
              'Asymptomatic','Hospitalized','Recovered','Fatalities'};

%% Colors and labels
alphaColors = {'m','g','b','r'};
alphaLabels = {'\alpha = 0.8','\alpha = 0.9','\alpha = 0.98','\alpha = 1'};

%% Create figure handles
hLine = cell(8, length(alpha_vals));
figHandles = zeros(8,1);

for k = 1:8
    if k == 1
    figHandles(k) = figure;
    hold on;
    xlabel('Time (days)');
    ylabel('Population');
    title(compTitles{k});
    xlim([0 Tmax]);
    end
end

%% ================= VIDEO SETUP =================
videoName = 'Covid_Fractional_Model.avi';
vidObj = VideoWriter(videoName);
vidObj.FrameRate = 10;
vidObj.Quality = 100;
open(vidObj);
%% ===============================================

%% Loop over alpha values
for i = 1:length(alpha_vals)
    
    alpha = alpha_vals(i);
    color = alphaColors{i};
    
    %% Parameters
    params = [2.55^alpha 1.56 7.65^alpha 0.25^alpha ...
              0.580 0.001 0.94^alpha 0.27^alpha ...
              0.5^alpha 3.5^alpha 1^alpha 0.3^alpha];
    
    %% RHS
    rhsFun = @(t,y,params) Covid_RHS_Wrapper(t,y,params);
    
    %% Solve fractional system
    [t, Yfrac] = Frac_Caputo_PECE(alpha, Tmax, h, y0', params, rhsFun);
    
    %% Initialize plot lines
    for k = 1:8
        if k ==1
        figure(figHandles(k));
        hLine{k,i} = plot(NaN, NaN, 'Color', color, 'LineWidth', 2);
        end
    end
    
    %% Real-time plotting + video recording
    for n = 1:length(t)
        if mod(n, updateStep) == 0 || n == length(t)
            
            idx = 1:updateStep:n;
            
            for k = 1:8
                if k == 1
                set(hLine{k,i}, 'XData', t(idx), 'YData', Yfrac(idx,k));
                end
            end
            
            drawnow;
            
            % Capture frame from first figure (Susceptible)
            frame = getframe(figHandles(1));
            writeVideo(vidObj, frame);
        end
    end
end

%% Add legends
for k = 1:8
    if k == 1
    figure(figHandles(k));
    legend([hLine{k,:}], alphaLabels, 'Location', 'best');
    end
end

%% Close video
close(vidObj);
disp('Video saved successfully.');
