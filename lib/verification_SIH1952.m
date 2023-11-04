%% verification Spectral Intensity according to Housner (1952)
% Calculate the spectral intensity as defined by Housner (1952) using
% OpenSeismoMatlab

%% Reference
% Housner, G. W. (1952). Intensity of ground motion during strong
% earthquakes. California Institute of Technology, Second technical report,
% Office of naval research, task order 25, project NR-081-095.

%% Description
% Table II on page 15 of the above reference is verified for the case
% of the El Centro, 1940 earthquake.The spectral intensity according to
% this table is 8.35ft for damping 0.0, 2.71ft for damping 0.2 and 1,89ft
% for damping 0.4. These values are verified with the present example

%% Load earthquake data
% Earthquake acceleration time history of the El Centro earthquake will be
% used (El Centro, 1940, El Centro Terminal Substation Building)
fid=fopen('elcentro_NS_trunc.dat','r');
text=textscan(fid,'%f %f');
fclose(fid);
t=text{1,1};
dt=t(2)-t(1);
xgtt=text{1,2};

%% Calculate Housner Spectral Intensity
% Switch
sw='SIH1952';

%%
% First value of damping (Table II in the above reference)
ksi1=0.0;

%%
% Second value of damping (Table II in the above reference)
ksi2=0.2;

%%
% Third value of damping (Table II in the above reference)
ksi3=0.4;

%%
% Apply OpenSeismoMatlab for calculating of the spectral intensity
S1=OpenSeismoMatlab(dt,xgtt,sw,ksi1);
S2=OpenSeismoMatlab(dt,xgtt,sw,ksi2);
S3=OpenSeismoMatlab(dt,xgtt,sw,ksi3);

%% Spectral intensities in ft
% For damping 0.0
S1.SI*3.28084

%%
% For damping 0.2
S2.SI*3.28084

%% 
% For damping 0.4
S3.SI*3.28084

%% Copyright
%
% Copyright (c) 2018-2023 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%

