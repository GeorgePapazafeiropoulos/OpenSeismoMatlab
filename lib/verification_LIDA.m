%% verification Linear dynamic response of SDOF oscillator
% Calculate the dynamic response of a linear SDOF oscillator. This example
% verifies Figure 6.6.1 in Chopra for Tn=1sec

%% Reference
% Chopra, A. K. (2020). Dynamics of structures, Theory and Applications to
% Earthquake Engineering, 5th edition. Prenctice Hall.

%% Earthquake motion
% Load earthquake data
dt=0.02;
fid=fopen('elcentro_NS_trunc.dat','r');
text=textscan(fid,'%f %f');
fclose(fid);
xgtt=text{1,2};

%% Setup parameters for LIDA function
% Eigenperiod
Tn=1;

%%
% Critical damping ratio
ksi=0.02;

%%
% Initial displacement
u0=0;

%%
% Initial velocity
ut0=0;

%%
% Algorithm to be used for the time integration
AlgID='U0-V0-Opt';

%%
% Minimum absolute value of the eigenvalues of the amplification matrix
rinf=1;

%% Calculate dynamic response
% Calculate circular eigenfrequency
omega=2*pi/Tn;

%%
% Apply LIDA
[u,ut,utt] = LIDA(dt,xgtt,omega,ksi,u0,ut0,AlgID,rinf);

%% Results
% Maximum displacement in cm
D=max(abs(u))*100

%% Copyright
%
% Copyright (c) 2018-2023 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%

