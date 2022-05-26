function [f,U] = FASp(dt,xgtt)
%
% Single sided Fourier amplitude spectrum
%
% [F,U] = FASP(DT,XGTT)
%
% Description
%     Fourier amplitude spectrum of an earthquake.
%
% Input parameters
%     DT [double(1 x 1)] is the time step of the input acceleration time
%         history XGTT.
%     XGTT [double(1:numsteps x 1)] is the input acceleration time history.
%         numsteps is the length of the input acceleration time history.
%
% Output parameters
%     F [double(1:2^(nextpow2(length(XGTT))-1) x 1)] is the frequency range
%         in which the Fourier amplitudes are calculated.
%     U [double(1:2^(nextpow2(length(XGTT))-1) x 1)] contains the Fourier
%         amplitudes
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<2
    error('Input arguments less than required')
end
if nargin>2
    error('Input arguments more than required')
end
% required inputs
if ~isscalar(dt)
    error('dt is not scalar')
end
if dt<=0
    error('dt is zero or negative')
end
if ~isvector(xgtt)
    error('xgtt is not vector')
end

%% Calculation
% Nyquist frequency (highest frequency)
Ny = (1/dt)/2; 
% number of points in xgtt
L  = length(xgtt); 
% Next power of 2 from length of xgtt
NFFT = 2^nextpow2(L);
% frequency spacing
df = 1/(NFFT*dt);
% Fourier amplitudes 
U = abs(fft(xgtt,NFFT))*dt; 
% Single sided Fourier amplitude spectrum
U = U(2:Ny/df+1);
% frequency range
f = linspace(df,Ny,Ny/df)'; 
end