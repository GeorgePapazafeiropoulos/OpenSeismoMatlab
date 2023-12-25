%% example sincResample
% Apply sinc resampling for an earthquake acceleration time history

%% Earthquake motion
% Load earthquake data
fid=fopen('elcentro_NS_full.dat','r');
text=textscan(fid,'%f %f');
fclose(fid);
time=text{1,1};
xgtt=text{1,2};
dt=time(2)-time(1);

%% Upsampling
% Define new time step for arbitrary upsampling ratio 1.15
dti=dt/1.15;
S1=OpenSeismoMatlab(dt,xgtt,'SINCRESAMPLE',dti);

%% Plot initial and upsampled acceleration time history
%
figure()
plot(time,xgtt,'-b','LineWidth',1)
hold on
plot(S1.time,S1.acc,'-r','LineWidth',1)
hold off
ylabel('Acceleration (m/s^2)')
xlabel('Time (sec)')
legend('Initial','Upsampled (1.15x)')
drawnow;
pause(0.1)

%% Downsampling
% Define new time step for arbitrary downsampling ratio 1.85
dti=dt*1.85;
S1=OpenSeismoMatlab(dt,xgtt,'SINCRESAMPLE',dti);

%% Plot initial and downsampled acceleration time history
%
figure()
plot(time,xgtt,'-b','LineWidth',1)
hold on
plot(S1.time,S1.acc,'-r','LineWidth',1)
hold off
ylabel('Acceleration (m/s^2)')
xlabel('Time (sec)')
legend('Initial','Downsampled (1.85x)')
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

