%% Linear elastic response spectra
% Generate the linear elastic response spectra of an earthquake suite using
% OpenSeismoMatlab.

%% Input
% Earthquake motions
eqmotions={'Imperial Valley'; % Imperial valley 1979
    'Kocaeli';
    'Loma Prieta';
    'Northridge';
    'San Fernando';
    'Spitak';
    'Cape Mendocino';
    'ChiChi';
    'El Centro'; % Imperial valley 1940
    'Hollister';
    'Kobe'};
%%
% Set the eigenperiod range for which the response spectra will be
% calculated.
Tspectra=(0.01:0.01:4)';
%%
% Set critical damping ratio of the response spectra to be calculated.
ksi=0.05;
%%
% Set the target ductility (not used here)
mu=2;
%%
% Extract linear elastic response spectra
sw='es';

%% Calculation
% Initialize LERS
LERS=cell(numel(eqmotions),1);
%%
% Calculation of peak values
for i=1:numel(eqmotions)
    % earthquake
    data=load([eqmotions{i},'.dat']);
    t=data(:,1);
    dt=t(2)-t(1);
    xgtt=data(:,2);
    S=OpenSeismoMatlab(dt,xgtt,sw,[],[],ksi,Tspectra,mu);
    LERS{i}=[S.Period,S.Sd,S.PSv,S.PSa];
end

%% Output
% Plot displacement response spectra
Fig1 = figure('units', 'centimeters', 'Position', [0,0, 20/sqrt(2), 20]);
% Scan all subplots
for i=1:numel(eqmotions)
    subplot(4,3,i)
    plot(LERS{i}(:,1),LERS{i}(:,2),'k','LineWidth',1);
    set(gca,'FontName','Times New Roman')
    title(eqmotions{i},'FontName','Times New Roman')
    ylabel('Sd','FontName','Times New Roman')
    xlabel('Period (s)','FontName','Times New Roman')
    axis tight
end
%%
% Plot pseudo-velocity response spectra
Fig2 = figure('units', 'centimeters', 'Position', [0,0, 20/sqrt(2), 20]);
% Scan all subplots
for i=1:numel(eqmotions)
    subplot(4,3,i)
    plot(LERS{i}(:,1),LERS{i}(:,3),'k','LineWidth',1);
    set(gca,'FontName','Times New Roman')
    title(eqmotions{i},'FontName','Times New Roman')
    ylabel('PSv','FontName','Times New Roman')
    xlabel('Period (s)','FontName','Times New Roman')
    axis tight
end
%%
% Plot pseudo-acceleration response spectra
Fig3 = figure('units', 'centimeters', 'Position', [0,0, 20/sqrt(2), 20]);
% Scan all subplots
for i=1:numel(eqmotions)
    subplot(4,3,i)
    plot(LERS{i}(:,1),LERS{i}(:,4),'k','LineWidth',1);
    set(gca,'FontName','Times New Roman')
    title(eqmotions{i},'FontName','Times New Roman')
    ylabel('PSa','FontName','Times New Roman')
    xlabel('Period (s)','FontName','Times New Roman')
    axis tight
end

%% Copyright
%
% Copyright (c) 2018-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%

