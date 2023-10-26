%% example RPSReSp
% Calculate rigid plastic sliding response spectra in OpenSeismoMatlab

%% Generate earthquake motion
% For reproducibility
rng(0)
%%
% Generate earthquake acceleration time history
dt=0.02;
N=10;
a=rand(N,1)-0.5;
b=100*pi*rand(N,1);
c=pi*(rand(N,1)-0.5);
t=(0:dt:(100*dt))';
xgtt=zeros(size(t));
for i=1:N
    xgtt=xgtt+a(i)*sin(b(i)*t+c(i));
end

%% Plot the generated time history
% 
figure()
plot(t,xgtt,'k','LineWidth',1)
ylabel('Acceleration (m/s^2)')
xlabel('Time (sec)')
title('Artificial acceleration time history')
drawnow;
pause(0.1)

%% Setup parameters for RPSReSp function
% Coulomb friction coefficients
CF=linspace(0.05,0.5,1000)';

%%
% Algorithm to be used for the time integration
AlgID='U0-V0-Opt';

%%
% Minimum absolute value of the eigenvalues of the amplification matrix
rinf=1;

%%
% Maximum tolerance of convergence for time integration algorithm
maxtol=0.01;

%%
% Maximum number of iterations per integration time step
jmax=200;

%%
% Infinitesimal acceleration
dak=eps;

%% Calculate spectra
% Apply RPSReSp
[Sd,Sv,Sa]=RPSReSp(dt,xgtt,CF,AlgID,rinf,maxtol,jmax,dak);

%% Plot the spectra
% Displacement spectrum
figure()
plot(CF,Sd,'k','LineWidth',1)
ylabel('Displacement (m)')
xlabel('Coulomb friction coefficient (-)')
drawnow;
pause(0.1)

%%
% Velocity spectrum
figure()
plot(CF,Sv,'k','LineWidth',1)
ylabel('Velocity (m/s)')
xlabel('Coulomb friction coefficient (-)')
drawnow;
pause(0.1)

%%
% Acceleration spectrum
figure()
plot(CF,Sa,'k','LineWidth',1)
ylabel('Acceleration (m/s^2)')
xlabel('Coulomb friction coefficient (-)')
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

