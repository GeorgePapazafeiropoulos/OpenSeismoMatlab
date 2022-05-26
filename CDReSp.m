function [PSa,PSv,Sd,Sv,Sa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,varargin)
%
% Constant Ductility Response Spectra (CDReSp)
%
% [PSA,PSV,SD,SV,SA,FYK,MUK,ITERK]=CDRESP(DT,XGTT,T,KSI,MU,N,TOL,PYSF,...
%     DTTOL,ALGID,RINF)
%
% Description
%     The constant ductility response spectra for a given time-history of
%     constant time step, a given eigenperiod range, a given viscous
%     damping ratio and a given ductility are computed. See section 7.5 in
%     Chopra (2012) and the notes "Inelastic Response Spectra" (CEE 541.
%     Structural Dynamics) by Henri P. Gavin.
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
%     MU [double(1 x 1)] is the target ductility for which the response
%         spectra are calculated.
%     N [double(1 x 1)] is the maximum number of iterations that can be
%         performed until convergence of the calculated ductility to the
%         target ductility is achieved. Default value 30.
%     TOL [double(1 x 1)] is the tolerance for convergence for the target
%         ductility. Default value 0.01.
%     PYSF [double(1 x 1)] is the post-yield stiffness factor, i.e. the
%         ratio of the postyield stiffness to the initial stiffness. PYSF=0
%         is not recommended for simulation of an elastoplastic system. A
%         small positive value is always suggested. PYSF is ignored if
%         MU=1. Default value 0.01.
%     DTTOL [double(1 x 1)] is the tolerance for resampling of the input
%         acceleration time history. For a given eigenperiod T, resampling
%         takes place if DT/T>dtTol. Default value 0.01.
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
%         Default value 'U0-V0-CA'.
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004). Default value 0.
%
% Output parameters
%     PSA [double(1:numSDOFs x 1)] is the Pseudo-Spectral Acceleration.
%     PSV [double(1:numSDOFs x 1)] is the Pseudo-Spectral Velocity.
%     SD [double(1:numSDOFs x 1)] is the Spectral Displacement.
%     SV [double(1:numSDOFs x 1)] is the Spectral Velocity.
%     SA [double(1:numSDOFs x 1)] is the Spectral Acceleration.
%     FYK [double(1:numSDOFs x 1)] is the yield limit that each SDOF must
%         have in order to attain ductility equal to muK.
%     MUK [double(1:numSDOFs x 1)] is the achieved ductility for each
%         period (each SDOF).
%     ITERK [double(1:numSDOFs x 1)] is the number of iterations needed for
%         convergence for each period (each SDOF).
%
% Example
%     %
%     dt=0.02;
%     %
%     N=10;
%     a=rand(N,1)-0.5;
%     b=100*pi*rand(N,1);
%     c=pi*(rand(N,1)-0.5);
%     t=(0:dt:(100*dt))';
%     xgtt=zeros(size(t));
%     for i=1:N
%         xgtt=xgtt+a(i)*sin(b(i)*t+c(i));
%     end
%     %
%     T=(0.04:0.04:4)';
%     %
%     ksi=0.05;
%     %
%     mu=2;
%     %
%     n=50;
%     %
%     tol=0.01;
%     %
%     pysf=0.1;
%     %
%     dtTol=0.02;
%     %
%     AlgID='U0-V0-Opt';
%     %
%     rinf=1;
%     %
%     [CDPSa,CDPSv,CDSd,CDSv,CDSa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,ksi,...
%         mu,n,tol,pysf,dtTol,AlgID,rinf);
%     %
%     figure()
%     plot(T,CDSd)
%     %
%     figure()
%     plot(T,fyK)
%     %
%     figure()
%     plot(T,muK)
%     %
%     figure()
%     plot(T,iterK)
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<3
    error('Input arguments less than required')
end
if nargin>11
    error('Input arguments more than required')
end
% set defaults for optional inputs
optargs = {0.05,2,50,0.01,0.01,0.01,'U0-V0-CA',0};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[ksi,mu,n,tol,pysf,dtTol,AlgID,rinf] = optargs{:};
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
if ~isvector(T)
    error('T is not vector')
end
if any(T<=0)
    error('T must be positive')
end
% optional inputs
if ~isscalar(ksi)
    error('ksi is not scalar')
end
if ksi<=0
    error('ksi is zero or negative')
end
if ~isscalar(mu)
    error('mu is not scalar')
end
if mu<1
    error('mu is lower than unity')
end
if ~isscalar(n)
    error('n is not scalar')
end
if n<10
    warning('n is lower than 10')
end
if ~isscalar(tol)
    error('tol is not scalar')
end
if tol<=0
    error('tol is zero or negative')
end
if ~isscalar(pysf)
    error('redf is not scalar')
end
if pysf<=0
    error('redf is zero or negative')
end
if ~isscalar(dtTol)
    error('dtTol is not scalar')
end
if dtTol<=0
    error('dtTol is zero or negative')
end
if ~isscalar(rinf)
    error('rinf is not scalar')
end
if rinf<0 || rinf>1
    error('rinf is lower than 0 or higher than 1')
end

%% Calculation
% Set integration constants
if all(size(AlgID)==[1,14])
    % define integration constants explicitly
    w1=AlgID(1);
    w2=AlgID(2);
    w3=AlgID(3);
    W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
    % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
    % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
    W1L1=AlgID(4);
    W2L2=AlgID(5);
    W3L3=AlgID(6);
    W1L4=AlgID(7);
    W2L5=AlgID(8);
    W1L6=AlgID(9);
    l1=AlgID(10);
    l2=AlgID(11);
    l3=AlgID(12);
    l4=AlgID(13);
    l5=AlgID(14);
else
    switch AlgID
        case 'U0-V0-Opt'
            % zero-order displacement & velocity overshooting behavior and
            % optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf); % suggested
            w2=15*(3-4*rinf)/(1-4*rinf); % suggested
            w3=-35*(1-rinf)/(1-4*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/2/(1+rinf)^2;
            W1L4=1/(1+rinf);
            W2L5=1/(1+rinf)^2; % suggested
            W1L6=(3-rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-CA'
            % zero-order displacement & velocity overshooting behavior and
            % continuous acceleration
            % rinf must belong to [1/3 1]
            if rinf<1/3
                rinf=1/3;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/3');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-5*rinf)/(3-7*rinf); % suggested
            w2=15*(1-13*rinf)/(3-7*rinf); % suggested
            w3=140*rinf/(3-7*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=(1+3*rinf)/4/(1+rinf);
            W3L3=(1+3*rinf)/4/(1+rinf)^2;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=(1+3*rinf)/2/(1+rinf)^2; % suggested
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-DA'
            % zero-order displacement & velocity overshooting behavior and
            % discontinuous acceleration
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15; % suggested
            w2=45; % suggested
            w3=-35; % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/2/(1+rinf);
            W1L4=1;
            W2L5=1/(1+rinf); % suggested
            W1L6=(3+rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V1-Opt'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'generalized a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-CA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'HHT a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-DA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'WBZ'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U1-V0-Opt'
            % first-order displacement & zero-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-8*rinf+6*rinf^2)/(9-22*rinf+19*rinf^2);
            w2=15*(25-74*rinf+53*rinf^2)/2/(9-22*rinf+19*rinf^2);
            w3=-35*(3-10*rinf+7*rinf^2)/(9-22*rinf+19*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3-rinf)/2/(1+rinf);
            W2L2=1/(1+rinf)^2;
            W3L3=1/(1+rinf)^3;
            W1L4=(3-rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^3;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-CA'
            % first-order displacement & zero-order velocity overshooting
            % behavior and continuous acceleration
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-60*(2-8*rinf+7*rinf^2)/(11-48*rinf+41*rinf^2);
            w2=15*(37-140*rinf+127*rinf^2)/2/(11-48*rinf+41*rinf^2);
            w3=-35*(5-18*rinf+17*rinf^2)/(11-48*rinf+41*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=2*rinf/(1+rinf)^2;
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=4*rinf/(1+rinf)^3;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-DA'
            % first-order displacement & zero-order velocity overshooting behavior
            % and discontinuous acceleration
            % This is the Newmark average acceleration a-form algorithm
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-4*rinf)/(9-11*rinf);
            w2=15*(25-37*rinf)/2/(9-11*rinf);
            w3=-35*(3-5*rinf)/(9-11*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3+rinf)/2/(1+rinf);
            W2L2=1/(1+rinf);
            W3L3=1/(1+rinf)^2;
            W1L4=(3+rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^2;
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'Newmark ACA'
            % Newmark Average Constant Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.25;
            W3L3=0.25;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.25;
            l4=1;
            l5=0.5;
        case 'Newmark LA'
            % Newmark Linear Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/6;
            W3L3=1/6;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/6;
            l4=1;
            l5=0.5;
        case 'Newmark BA'
            % Newmark Backward Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.5;
            W3L3=0.5;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.5;
            l4=1;
            l5=0.5;
        case 'Fox-Goodwin'
            % Fox-Goodwin formula
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/12;
            W3L3=1/12;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/12;
            l4=1;
            l5=0.5;
        otherwise
            error('No appropriate algorithm specified.');
    end
end
% initialize
NumSDOF=length(T);
Sd=nan(NumSDOF,1);
Sv=nan(NumSDOF,1);
Sa=nan(NumSDOF,1);
PSv=nan(NumSDOF,1);
PSa=nan(NumSDOF,1);
fyK=nan(NumSDOF,1);
uyK=nan(NumSDOF,1);
muK=nan(NumSDOF,1);
iterK=zeros(NumSDOF,1);
% natural frequencies of SDOFs to be analysed
omega=2*pi./T;
% Flip eigenfrequency vector in order for the half-stepping algorithm
% (HalfStep function) to work from large to small eigenperiods
omega=omega(end:-1:1);
% set initial conditions
u0=0;
ut0=0;
% set SDOF mass to unity (any value can be specified without change of the
% results)
m=1;
% Maximum tolerance for convergence
tol2=0.01;
% Maximum number of iterations per increment
jmax=200;
% Infinitesimal variation of acceleration
dak=eps;
for j=1:length(T)
    % Step 3, 7.5.3 in Chopra
    omegaj=omega(j);
    % Check if dt/T>dtTol. If yes, then reproduce the time history with the
    % half step
    if dt*omegaj/(2*pi)>dtTol
        xgtt=HalfStep(xgtt);
        dt=dt/2;
    end
    k_hi=m*omegaj^2;
    if mu==1
        k_lo=k_hi;
    else
        k_lo=k_hi*pysf;
    end
    % Step 4, 7.5.3 in Chopra
    [u,ut,utt,p,Es,Ed] = DRHA(k_hi,m,dt,xgtt,ksi,u0,ut0,AlgID,rinf);
    upeak=max(abs(u));
    fpeak=k_hi*upeak;
    
    % Step 5, 7.5.3 in Chopra
    % Find uy in order to achieve ductility equal to mu
    
    % The following steps are setup according to the notes "Inelastic
    % Response Spectra" (CEE 541. Structural Dynamics) by Henri P. Gavin
    
    % Step 4.(b) in Gavin
    % Maximum yield strength coefficient Cymax=2
    maxuy = max(abs(u))*0.999;
    % Calculate minimum ductility demand mumin
    % ductility factor from eq. 7.2.4. in Chopra
    [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
        k_lo,maxuy,ksi,AlgID,u0,ut0,rinf,tol2,jmax,dak);
    umax=max(abs(um));
    fy=k_hi*maxuy;
    fybark=fy/fpeak;
    mumin=(umax/upeak)/fybark;
    % Step 4.(a) in Gavin
    % Minimum yield strength coefficient Cymin=0.01
    minuy = max(abs(u))/(mu*5);
    % Calculate minimum ductility demand mumax
    % ductility factor from eq. 7.2.4. in Chopra
    [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
        k_lo,minuy,ksi,AlgID,u0,ut0,rinf,tol2,jmax,dak);
    umax=max(abs(um));
    fy=k_hi*minuy;
    fybark=fy/fpeak;
    mumax=(umax/upeak)/fybark;
    % Step 4.(c)
    % Hyperbola coefficients [mu=alpha/(uy+beta)]
    alpha=mumin*mumax*(maxuy-minuy)/(mumax-mumin);
    beta=(maxuy*mumin-minuy*mumax)/(mumax-mumin);
    % Step 4.(d)
    % Determine yield strengths corresponding to 1.2*mu and 0.8*mu
    % according to the hyperbola of Step 4.(c)
    uy1=max(alpha/(1.2*mu)-beta,minuy);
    uy2=min(alpha/(0.8*mu)-beta,maxuy);
    % Convergence iterations
    % k = iteration number
    for k=1:n
        % Step 5.(b) Solve for uy1
        [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
            k_lo,uy1,ksi,AlgID,u0,ut0,rinf,tol2,jmax,dak);
        umax=max(abs(um));
        % ductility factor from eq. 7.2.4. in Chopra
        fy=k_hi*uy1;
        fybark=fy/fpeak;
        mu1=(umax/upeak)/fybark;
        % Step 5.(c) Solve for uy2
        [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
            k_lo,uy2,ksi,AlgID,u0,ut0,rinf,tol2,jmax,dak);
        umax=max(abs(um));
        % ductility factor from eq. 7.2.4. in Chopra
        fy=k_hi*uy2;
        fybark=fy/fpeak;
        mu2=(umax/upeak)/fybark;
        % Step 5.(d) Slope of line connecting (uy1,mu1) and (uy2,mu2)
        S=(mu2-mu1)/(uy2-uy1);
        % Step 5.(e) Evaluate Duy
        Duy=min(abs(mu-mu2)/S,0.1*(uy1+uy2));
        % Step 5.(f)
        if (mu-mu2)/S>0
            Duy=-Duy;
        end
        % Step 5.(g)
        if S>0 && mu2<mu
            Duy=-0.5*uy2;
        end
        % Step 5.(h)
        if S>0 && mu2>mu
            Duy=0.1*uy2;
        end
        % Step 5.(i)
        if uy2+Duy<0
            Duy=(minuy-uy2)*(mu-mu2)/(mumax-mu2);
        end
        % Step 5.(j)
        uy1=uy2;
        mu1=mu2;
        % Step 5.(k)
        uy2=uy2+Duy;
        % Step 5.(l)
        [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
            k_lo,uy2,ksi,AlgID,u0,ut0,rinf,tol2,jmax,dak);
        umax=max(abs(um));
        % ductility factor from eq. 7.2.4. in Chopra
        fy=k_hi*uy2;
        fybark=fy/fpeak;
        mu2=(umax/upeak)/fybark;
        % Store the output values of the current iteration
        muK(j)=mu2;
        fyK(j)=fy;
        uyK(j)=uy2;
        iterK(j)=k;
        % Step 5.(m) Check for convergence
        if abs(uy1-uy2)<1e-4*tol || abs(mu1-mu2)<1e-4*tol || abs(mu2-mu)<tol
            break
        end
    end
    % find Sd, Sv, Sa
    Sd(j)=umax;
    Sv(j)=max(abs(umt));
    Sa(j)=max(abs(umtt));
    % find PSv, PSa
    PSv(j)=Sd(j)*omegaj;
    PSa(j)=Sd(j)*omegaj^2;
end
% Flip output quantities to be compatible with omega
PSa=PSa(end:-1:1);
PSv=PSv(end:-1:1);
Sd=Sd(end:-1:1);
Sv=Sv(end:-1:1);
Sa=Sa(end:-1:1);
fyK=fyK(end:-1:1);
muK=muK(end:-1:1);
iterK=iterK(end:-1:1);

end

