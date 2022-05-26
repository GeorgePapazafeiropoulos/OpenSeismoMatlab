%% Verification
% This is to verify that OpenSeismoMatlab works properly for all of its
% possible options and types of application

%% Load earthquake data
% Earthquake acceleration time history of the El Centro earthquake will be
% used (El Centro, 1940, El Centro Terminal Substation Building)
fid=fopen('elcentro.dat','r');
text=textscan(fid,'%f %f');
fclose(fid);
t=text{1,1};
dt=t(2)-t(1);
xgtt=9.81*text{1,2};
%% Time histories without baseline correction
% 
sw='timehist';
baselineSw=false;
S1=OpenSeismoMatlab(dt,xgtt,sw,baselineSw);
%%
figure(1)
plot(S1.time,S1.disp,'k','LineWidth',1)
%%
figure(2)
plot(S1.time,S1.vel,'k','LineWidth',1)
%%
figure(3)
plot(S1.time,S1.acc,'k','LineWidth',1)

%% Time histories with baseline correction
% 
sw='timehist';
baselineSw=true;
S2=OpenSeismoMatlab(dt,xgtt,sw,baselineSw);
%% 
figure(4)
plot(S2.time,S2.disp,'k','LineWidth',1)
%% 
figure(5)
plot(S2.time,S2.vel,'k','LineWidth',1)
%% 
figure(6)
plot(S2.time,S2.acc,'k','LineWidth',1)

%% Resample acceleration time history from 0.02 sec to 0.01 sec.
%
sw='resample';
dti=0.01;
S3=OpenSeismoMatlab(dt,xgtt,sw,[],dti);
%% 
figure(7)
plot(S3.time,S3.acc,'k','LineWidth',1)

%% PGA
%
sw='pga';
S4=OpenSeismoMatlab(dt,xgtt,sw);
%%
S4.PGA

%% PGV
%
sw='pgv';
S5=OpenSeismoMatlab(dt,xgtt,sw);
%%
S5.PGV

%% PGD
%
sw='pgd';
S6=OpenSeismoMatlab(dt,xgtt,sw);
%%
S6.PGD

%% Arias intensity and significant duration
%
sw='arias';
S7=OpenSeismoMatlab(dt,xgtt,sw);
%%
S7.Ecum
%% 
figure(8)
plot(S7.time,S7.EcumTH,'k','LineWidth',1)
%%
S7.t_5_95
%%
S7.Td
%%
S7.arias

%% Linear elastic response spectra and pseudospectra
%
sw='es';
ksi=0.05;
T=0.04:0.02:1;
S8=OpenSeismoMatlab(dt,xgtt,sw,[],[],ksi,T);
%% 
figure(9)
plot(S8.Period,S8.PSa,'k','LineWidth',1)
%% 
figure(10)
plot(S8.Period,S8.PSv,'k','LineWidth',1)
%% 
figure(11)
plot(S8.Period,S8.Sd,'k','LineWidth',1)
%% 
figure(12)
plot(S8.Period,S8.Sv,'k','LineWidth',1)
%% 
figure(13)
plot(S8.Period,S8.Sa,'k','LineWidth',1)
%% 
figure(14)
plot(S8.Period,S8.SievABS,'k','LineWidth',1)
%% 
figure(15)
plot(S8.Period,S8.SievREL,'k','LineWidth',1)
%%
S8.PredPSa

%% 
S8.PredPeriod

%% Constant ductility response spectra and pseudospectra
%
sw='cds';
ksi=0.05;
T=0.04:0.02:1;
mu=2;
S9=OpenSeismoMatlab(dt,xgtt,sw,[],[],ksi,T,mu);
%% 
figure(15)
plot(S9.Period,S9.CDPSa,'k','LineWidth',1)
%% 
figure(16)
plot(S9.Period,S9.CDPSv,'k','LineWidth',1)
%% 
figure(17)
plot(S9.Period,S9.CDSd,'k','LineWidth',1)
%% 
figure(18)
plot(S9.Period,S9.CDSv,'k','LineWidth',1)
%% 
figure(19)
plot(S9.Period,S9.CDSa,'k','LineWidth',1)
%% 
figure(20)
plot(S9.Period,S9.fyK,'k','LineWidth',1)
%% 
figure(21)
plot(S9.Period,S9.muK,'k','LineWidth',1)
%% 
figure(22)
plot(S9.Period,S9.iterK,'k','LineWidth',1)

%% Fourier amplitude spectrum and mean period
%
sw='fas';
S10=OpenSeismoMatlab(dt,xgtt,sw);
%% 
figure(23)
plot(S10.freq,S10.FAS,'k','LineWidth',1)
%%
S10.Tm
%%
S10.Fm

%% Copyright
%
% Copyright (c) 2018-2022 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%



