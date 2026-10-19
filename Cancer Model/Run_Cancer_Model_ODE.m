clear; clc; close all;

C = [1 2.5 0.6 1.5 4.5 1 0.2 0.5];
X0 = [0.1; 0.1; 0.1];

tspan = [0 300];
dt = 0.01;

T = tspan(1):dt:tspan(2);
N = length(T);

X = zeros(N,3);
X(1,:) = X0';

%% ===================== 2D VIDEO =====================
v = VideoWriter('myVideo1.mp4','MPEG-4');
v.FrameRate = 60;
v.Quality = 100;
open(v);

fig2D = figure('Color','white');        %%% FIX
h1 = plot(NaN,NaN,'r','LineWidth',2); hold on;
h2 = plot(NaN,NaN,'g','LineWidth',2);
h3 = plot(NaN,NaN,'b','LineWidth',2);
xlabel('Time'); ylabel('Population');
title('Integer-order Cancer Model (ODE)');
axis([0 300 0 1]);

%% ===================== 3D VIDEO =====================
u = VideoWriter('myVideo2.mp4','MPEG-4');
u.FrameRate = 60;
u.Quality = 100;
open(u);

fig3D = figure('Color','white');         %%% FIX
h3d = plot3(NaN,NaN,NaN,'b','LineWidth',2);
xlabel('x'); ylabel('y'); zlabel('z');
title('Integer-order Cancer Model (ODE)');
axis([0 1 0 1 0 1]);
view(-45,20);

%% ===================== SOLVER =====================
disp('Running Integer ODE Solver...');
tic;
X_temp = X0;
updateStep = 250;

for k = 1:N-1

    tspan_step = [T(k) T(k+1)];
    [~, X_step] = ode45(@(t,X) Cancer_Model(t,X,C), tspan_step, X_temp);
    
    X(k+1,:) = X_step(end,:);
    X_temp = X(k+1,:).';

    if mod(k, updateStep) == 0
        
        % ----- 2D update -----
        set(h1,'XData',T(1:k),'YData',X(1:k,1));
        set(h2,'XData',T(1:k),'YData',X(1:k,2));
        set(h3,'XData',T(1:k),'YData',X(1:k,3));
        
        frame1 = getframe(fig2D);     %%% FIX
        writeVideo(v, frame1);

        % ----- 3D update -----
        set(h3d,'XData',X(1:k,1),'YData',X(1:k,2),'ZData',X(1:k,3));
        
        frame2 = getframe(fig3D);     %%% FIX
        writeVideo(u, frame2);

        drawnow;
    end
end

cpuODE = toc;
fprintf('CPU Time (ODE) = %.6f seconds\n', cpuODE);

%% ===================== CLOSE VIDEOS =====================
close(v);     %%% REQUIRED
close(u);     %%% REQUIRED
disp('MP4 videos saved successfully!');
