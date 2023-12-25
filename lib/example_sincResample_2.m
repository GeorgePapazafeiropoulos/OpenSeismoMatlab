%% example sincResample
% Apply sinc resampling for a sinusoidal acceleration time history

%% Sinusoidal motion
% Develop the sinusoidal motion
n = 200;
dt = 1/20;
t = dt*(0:n-1);
T = dt*n;
y = sin(2*pi*15*t/T);

%% Upsampling
% Define new time step for arbitrary upsampling ratio 1.15
dti=dt/1.15;
S1=OpenSeismoMatlab(dt,y,'SINCRESAMPLE',dti);

%% Downsampling
% Define new time step for arbitrary downsampling ratio 1.85
dti=dt*1.85;
S2=OpenSeismoMatlab(dt,y,'SINCRESAMPLE',dti);

%% Plot initial and upsampled acceleration time history
%
figure()
plot(t,y,'-ob');
hold on;
plot(S1.time,S1.acc,'-*r');
plot(S2.time,S2.acc,'-og');
hold off
ylabel('Acceleration (m/s^2)')
xlabel('Time (sec)')
legend('Initial','Upsampled (1.15x)','Downsampled (1.85x)')
drawnow;
pause(0.1)

%% Copyright
%
% Copyright (c) 2018-2023 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%

