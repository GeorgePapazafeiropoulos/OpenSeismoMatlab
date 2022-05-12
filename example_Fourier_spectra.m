%% Fourier spectra
% Generate the Fourier spectra of an earthquake suite using
% OpenSeismoMatlab.

%% Input
% earthquake motions
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
Tspectra=(0.02:0.01:4)';
%%
% Set critical damping ratio of the response spectra to be calculated.
ksi=0.05;
%%
% Set the target ductility (not used here)
mu=2;
%%
% Extract fourier spectra
sw='fas';

%% Calculation
% Initialize Fourier
Fourier=cell(numel(eqmotions),1);
% Calculation of peak values
for i=1:numel(eqmotions)
    % earthquake
    data=load([eqmotions{i},'.dat']);
    t=data(:,1);
    dt=t(2)-t(1);
    xgtt=data(:,2);
    S=OpenSeismoMatlab(dt,xgtt,sw);
    Fourier{i}=[S.freq,S.FAS];
end

%% Output
% Plot Fourier amplitude
Fig1 = figure('units', 'centimeters', 'Position', [0,0, 20/sqrt(2), 20]);
% Scan all subplots
for i=1:numel(eqmotions)
    subplot(4,3,i)
    plot(Fourier{i}(:,1),Fourier{i}(:,2),'k','LineWidth',1);
    set(gca,'FontName','Times New Roman')
    title(eqmotions{i},'FontName','Times New Roman')
    ylabel('FAS','FontName','Times New Roman')
    xlabel('Frequency (Hz)','FontName','Times New Roman')
    axis tight
end

%% Copyright
%
% Copyright (c) 2018-2022 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D.
% * Email: gpapazafeiropoulos@yahoo.gr
%

