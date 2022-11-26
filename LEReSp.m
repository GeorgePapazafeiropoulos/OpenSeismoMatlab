function [PSa,PSv,Sd,Sv,Sa,SievABS,SievREL]=LEReSp(dt,xgtt,T,ksi,dtTol,...
    AlgID,rinf)
%
% Linear Elastic Response Spectra
%
% [PSA,PSV,SD,SV,SA,SIEVABS,SIEVREL]=LERESP(DT,XGTT,T,KSI,DTTOL,...
%     ALGID,RINF)
%
% Description
%     The linear elastic response spectra for a given time-history of
%     constant time step, a given eigenperiod range and a given viscous
%     damping ratio are computed. These spectra include the spectral
%     acceleration, spectral velocity, spectral displacement,
%     pseudoacceleration, pseudovelocity, absolute equivalent input energy
%     velocity and relative equivalent input energy velocity. This function
%     is part of the OpenSeismoMatlab software. It can be used as
%     standalone, however attention is needed for the correctness of the
%     input arguments, since no checks are performed in this function. See
%     the example example_LEReSp.m for more details about how this function
%     can be implemented.
%
% Input parameters
%     DT [double(1 x 1)] is the time step of the input acceleration time
%         history XGTT.
%     XGTT [double(1:numsteps x 1)] is the input acceleration time history.
%         numsteps is the length of the input acceleration time history.
%     T [double(1:numSDOFs x 1)] contains the values of eigenperiods for
%         which the response spectra are requested. numSDOFs is the number
%         of SDOF oscillators being analysed to produce the spectra.
%     KSI [double(1 x 1)] is the fraction of critical viscous damping.
%     DTTOL [double(1 x 1)] is the maximum ratio of the integration time
%         step to the eigenperiod.
%     ALGID [char(1 x :inf)] is the algorithm to be used for the time
%         integration. It can be one of the following strings for superior
%         optimally designed algorithms:
%             'generalized a-method': The generalized a-method (Chung &
%             Hulbert, 1993)
%             'HHT a-method': The Hilber-Hughes-Taylor method (Hilber,
%             Hughes & Taylor, 1977)
%             'WBZ': The Wood–Bossak–Zienkiewicz method (Wood, Bossak &
%             Zienkiewicz, 1980)
%             'U0-V0-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement zero order velocity algorithm
%             'U0-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement zero order
%             velocity algorithm
%             'U0-V1-Opt': Optimal numerical dissipation and dispersion
%             zero order displacement first order velocity algorithm
%             'U0-V1-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) zero order displacement first order
%             velocity algorithm
%             'U0-V1-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) zero order displacement first order
%             velocity algorithm
%             'U1-V0-Opt': Optimal numerical dissipation and dispersion
%             first order displacement zero order velocity algorithm
%             'U1-V0-CA': Continuous acceleration (zero spurious root at
%             the low frequency limit) first order displacement zero order
%             velocity algorithm
%             'U1-V0-DA': Discontinuous acceleration (zero spurious root at
%             the high frequency limit) first order displacement zero order
%             velocity algorithm
%             'Newmark ACA': Newmark Average Constant Acceleration method
%             'Newmark LA': Newmark Linear Acceleration method
%             'Newmark BA': Newmark Backward Acceleration method
%             'Fox-Goodwin': Fox-Goodwin formula
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004).
%
% Output parameters
%     PSA [double(1:numSDOFs x 1)] is the Pseudo Acceleration Spectrum.
%     PSV [double(1:numSDOFs x 1)] is the Pseudo Velocity Spectrum.
%     SD [double(1:numSDOFs x 1)] is the Spectral Displacement.
%     SV [double(1:numSDOFs x 1)] is the Spectral Velocity.
%     SA [double(1:numSDOFs x 1)] is the Spectral Acceleration.
%     SIEVABS [double(1:numSDOFs x 1)] is the equivalent absolute input
%         energy velocity.
%     SIEVREL [double(1:numSDOFs x 1)] is the equivalent relative input
%         energy velocity.
%
% Example
%     dt=0.02;
%     N=10;
%     a=rand(N,1)-0.5;
%     b=100*pi*rand(N,1);
%     c=pi*(rand(N,1)-0.5);
%     t=(0:dt:(100*dt))';
%     xgtt=zeros(size(t));
%     for i=1:N
%         xgtt=xgtt+a(i)*sin(b(i)*t+c(i));
%     end
%     T=logspace(log10(0.02),log10(50),1000)';
%     ksi=0.05;
%     dtTol=0.02;
%     AlgID='U0-V0-Opt';
%     rinf=1;
%     [PSa,PSv,Sd,Sv,Sa,SievABS,SievREL]=LEReSp(dt,xgtt,T,ksi,dtTol,...
%         AlgID,rinf);
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Calculation
% Set integration constants
[w1,w2,w3,W1,W1L1,W2L2,W3L3,W1L4,W2L5,W1L6,l1,l2,l3,l4,l5,rinf] = ...
    TimeIntConstants(AlgID,rinf);
% Initialize
NumSDOF=length(T);
Sd=zeros(NumSDOF,1);
Sv=zeros(NumSDOF,1);
Sa=zeros(NumSDOF,1);
SievABS=zeros(NumSDOF,1);
SievREL=zeros(NumSDOF,1);
% Set the eigenfrequencies of the SDOF population
omega=2*pi./T;
% Flip eigenfrequency vector in order for the half-stepping algorithm
% (HalfStep function) to work from large to small eigenperiods
omega=omega(end:-1:1);
% set initial conditions
u0=0;
ut0=0;
for j=1:length(T)
    omegaj=omega(j);
    % Check if dt/T>dtTol. If yes, then reproduce the time history with the
    % half step
    if dt*omegaj/(2*pi)>dtTol
        xgtt=HalfStep(xgtt);
        dt=dt/2;
    end
    [u,ut,utt] = LIDA(dt,xgtt,omegaj,ksi,u0,ut0,AlgID,rinf);
    % output
    Sd(j)=max(abs(u));
    Sv(j)=max(abs(ut));
    Sa(j)=max(abs(utt));
    SievABS(j)=sqrt(2*dt^2*sum((utt+xgtt).*cumsum(xgtt)));
    SievREL(j)=sqrt(-2*dt^2*sum((xgtt).*cumsum(utt)));
end
% Flip output quantities to be compatible with omega
omega=omega(end:-1:1);
Sd=Sd(end:-1:1);
Sv=Sv(end:-1:1);
Sa=Sa(end:-1:1);
SievABS=SievABS(end:-1:1);
SievREL=SievREL(end:-1:1);
% Calculate pseudovelocity and pseudoacceleration
PSv=Sd.*omega;
PSa=Sd.*omega.^2;

end
