function [u,ut,utt,Fs,Ey,Es,Ed,jiter] = NLIDABLKIN(dt,xgtt,m,k_hi,k_lo,uy,...
    varargin)
%
% Non Linear Implicit Dynamic Analysis of a bilinear kinematic hardening
% hysteretic structure with elastic damping
% 
% [U,UT,UTT,FS,EY,ES,ED,JITER] = NLIDABLKIN(DT,XGTT,M,K_HI,K_LO,UY,KSI,...
%     ALGID,U0,UT0,RINF,MAXTOL,JMAX,DAK)
%
% Description
%     General linear implicit direct time integration of second order
%     differential equations of a bilinear elastoplastic hysteretic SDOF
%     dynamic system with elastic damping, with lumped mass.
%     The General Single Step Single Solve (GSSSS) family of algorithms
%     published by X.Zhou & K.K.Tamma (2004) is employed for direct time
%     integration of the general linear or nonlinear structural Single
%     Degree of Freedom (SDOF) dynamic problem. Selection among 9
%     algorithms, all designed according to the above journal article, can
%     be made in this routine. These algorithms encompass the scope of
%     Linear Multi-Step (LMS) methods and are limited by the Dahlquist
%     barrier theorem (Dahlquist,1963).
%
% Input parameters
%     DT [double(1 x 1)] is the time step of the integration
%     XGTT [double(1:NumSteps x 1)] is the acceleration time history which
%         is imposed at the lumped mass of the SDOF structure.
%     M [double(1 x 1)] is the lumped masses of the structure. Define the
%         lumped masses from the top to the bottom, excluding the fixed dof
%         at the base
%     K_HI [double(1 x 1)] is the initial stiffness of the system before
%         its first yield, i.e. the high stiffness. Give the stiffness of
%         each storey from top to bottom.
%     K_LO [double(1 x 1)] is the post-yield stiffness of the system,
%         i.e. the low stiffness. Give the stiffness of each storey from
%         top to bottom.
%     UY [double(1 x 1)] is the yield limit of the stiffness elements of
%         the structure. The element is considered to yield, if the
%         interstorey drift between degrees of freedom i and i+1 exceeds
%         UY(i). Give the yield limit of each storey from top to bottom.
%     KSI [double(1 x 1)] is the ratio of critical viscous damping of the
%         system, assumed to be unique for all damping elements of the
%         structure.
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
%     U0 [double(1 x 1)] is the initial displacement.
%     UT0 [double(1 x 1)] is the initial velocity.
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004).
%     MAXTOL [double(1 x 1)] is the maximum tolerance of convergence of the
%         Full Newton Raphson method for numerical computation of
%         acceleration.
%     JMAX [double(1 x 1)] is the maximum number of iterations per
%         increment. If JMAX=0 then iterations are not performed and the
%         MAXTOL parameter is not taken into account.
%     DAK [double(1 x 1)] is the infinitesimal acceleration for the
%         calculation of the derivetive required for the convergence of the
%         Newton-Raphson iteration.
%
% Output parameters
%     U [double(1 x 1:NumSteps)] is the time-history of displacement
%     UT [double(1 x 1:NumSteps)] is the time-history of velocity
%     UTT [double(1 x 1:NumSteps)] is the time-history of acceleration
%     FS [double(1 x 1:NumSteps)] is the time-history of the internal
%         force of the structure analysed.
%     EY [double(1 x 1:NumSteps)] is the time history of the sum of the
%         energy dissipated by yielding during each time step and the
%         recoverable strain energy of the system (incremental).
%         cumsum(EY)-Es gives the time history of the total energy
%         dissipated by yielding from the start of the dynamic analysis.
%     ES [double(1 x 1:NumSteps)] is the time-history of the recoverable
%         strain energy of the system (total and not incremental).
%     ED [double(1 x 1:NumSteps)] is the time-history of the energy
%         dissipated by viscoelastic damping during each time step
%         (incremental). cumsum(ED) gives the time history of the total
%         energy dissipated from the start of the dynamic analysis.
%     JITER [double(1 x 1:NumSteps)] is the iterations per increment
%
% Notation in the code
%     u=displacement
%     un=displacement after increment n
%     ut=velocity
%     utn=velocity after increment n
%     utt=acceleration
%     uttn=acceleration after increment n
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<6
    error('Input arguments less than required')
end
if nargin>14
    error('Input arguments more than required')
end
% set defaults for optional inputs
optargs = {0.05,'U0-V0-Opt',0,0,1,0.01,200,eps};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak] = optargs{:};
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
if ~isscalar(m)
    error('m is not scalar')
end
if m<=0
    error('m is zero or negative')
end
if ~isscalar(k_hi)
    error('k_hi is not scalar')
end
if k_hi<=0
    error('k_hi is zero or negative')
end
if ~isscalar(k_lo)
    error('k_lo is not scalar')
end
if k_lo<=0
    error('k_lo is zero or negative')
end
if ~isscalar(uy)
    error('uy is not scalar')
end
if uy<=0
    error('uy is zero or negative')
end
% optional inputs
if ~isscalar(ksi)
    error('ksi is not scalar')
end
if ksi<=0
    error('ksi is zero or negative')
end
if ~isscalar(u0)
    error('u0 is not scalar')
end
if ~isscalar(ut0)
    error('ut0 is not scalar')
end
if ~isscalar(rinf)
    error('rinf is not scalar')
end
if rinf<0 || rinf>1
    error('rinf is lower than 0 or higher than 1')
end
if ~isscalar(maxtol)
    error('maxtol is not scalar')
end
if maxtol<=0
    error('maxtol is zero or negative')
end
if ~isscalar(jmax)
    error('jmax is not scalar')
end
if jmax<0
    error('jmax is negative')
end
if floor(jmax)~=jmax
    error('jmax is not integer')
end
if ~isscalar(dak)
    error('dak is not scalar')
end
if dak<=0
    error('dak is zero or negative')
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

%% Calculation
% number of analysis increments
NumSteps=length(xgtt);
% Initialize output
u=zeros(1,NumSteps);
ut=zeros(1,NumSteps);
utt=zeros(1,NumSteps);
Fs=zeros(1,NumSteps);
Ey=zeros(1,NumSteps);
Es=zeros(1,NumSteps);
Ed=zeros(1,NumSteps);
jiter=zeros(1,NumSteps);
% set initial values of displacement, velocity, acceleration (u0,ut0 and
% utt0 respectively) at n=0.
u(1)=u0;
ut(1)=ut0;
% construct lumped mass matrix
M=diag(m,0);
% calculation for first increment
k=k_hi;
d=0;
[FKC0,K0,C0,k,d] = BLKIN(u0,ut0,k_hi,k_lo,uy,M,ksi,k,d);
utt0=-xgtt(1)-M\FKC0;
utt(1)=utt0;
Fs(1)=FKC0;
% initial assignments
FKCn=FKC0;
Kn=K0;
Cn=C0;
un=u0;
utn=ut0;
uttn=utt0;
% integration increments n
for n=1:NumSteps-1
    % effective force
    Feffn1k=-FKCn...
        -Kn*(W1L1*dt*utn+W2L2*dt^2*uttn)...
        -W1L4*dt*Cn*uttn...
        +M*(((1-W1)*xgtt(n)+W1*xgtt(n+1))-uttn);
    % effective mass
    Meffn=W1L6*M+W2L5*dt*Cn+W3L3*dt^2*Kn;
    % initial estimate of da
    dan=Meffn\Feffn1k;
    
    % start iteration number k
    j=1;
    % set initial quotient of variation of da equal to maxtol
    quda=maxtol;
    % full Newton-Raphson iterations k
    while max(abs(quda))>=maxtol && j<=jmax
        % iteration k+1 of increment n+1
        
        % _________________________________________________________________
        %/
        % displacement, velocity, acceleration, internal force, stiffness
        % and damping for uttn+dan
        % calculate the residual Rn1k at uttn+dan
        % update kinematic quantities
        un1k=un+l1*utn*dt+l2*uttn*dt^2+l3*dan*dt^2;
        utn1k=utn+l4*uttn*dt+l5*dan*dt;
        uttn1k=uttn+dan;
        % force due to stiffness and damping
        [FKCn1k,Kn1k,Cn1k,~,~]=BLKIN(un1k,utn1k,k_hi,k_lo,uy,M,ksi,k,d);
        % effective force
        Feffn1k=-FKCn1k...
            -Kn1k*(W1L1*dt*utn1k+W2L2*dt^2*uttn1k)...
            -Cn1k*W1L4*dt*uttn1k...
            +M*((1-W1)*xgtt(n)+W1*xgtt(n+1)-uttn1k);
        % effective mass
        Meffn1k=Kn1k*W3L3*dt^2+Cn1k*W2L5*dt+M*W1L6;
        % residual
        Rn1k=Feffn1k-Meffn1k*dan;
        %\_________________________________________________________________

        % _________________________________________________________________
        %/
        % displacement, velocity, acceleration, internal force, stiffness
        % and damping for uttn+(dan+dak)
        % calculate the derivative at uttn+dan as:
        % dR/da=(dRn1k-Rn1k)/(uttn+(dan+dak)-(uttn+dan))=(dRn1k-Rn1k)/dak
        % update kinematic quantities
        dun1k=un+l1*utn*dt+l2*uttn*dt^2+l3*(dan+dak)*dt^2;
        dutn1k=utn+l4*uttn*dt+l5*(dan+dak)*dt;
        duttn1k=uttn+(dan+dak);

        % force due to stiffness and damping
        [dFKCn1k,dKn1k,dCn1k,~,~]=BLKIN(dun1k,dutn1k,k_hi,k_lo,uy,M,ksi,k,d);
        % effective force
        dFeffn1k=-dFKCn1k...
            -dKn1k*(W1L1*dt*dutn1k+W2L2*dt^2*duttn1k)...
            -dCn1k*W1L4*dt*duttn1k...
            +M*((1-W1)*xgtt(n)+W1*xgtt(n+1)-duttn1k);
        % effective mass
        dMeffn1k=dKn1k*W3L3*dt^2+dCn1k*W2L5*dt+M*W1L6;
        % residual
        dRn1k=dFeffn1k-dMeffn1k*duttn1k;
        %\_________________________________________________________________
        
        % Full Newton-Raphson update:
        % da_new=da-Rn1k/(dR/da)=da*(1-Rn1k/(dRn1k/dak)/da)
        % (to be checked for while loop termination)
        quda=(Rn1k./(dRn1k-Rn1k).*dak)./dan;
        % test if derivative becomes zero
        a=isinf(quda);
        b=isnan(quda);
        if any(a) || any(b)
            break
            %quda=zeros(size(quda));
        end
        % update da
        dan=(1-quda).*dan;
        % update iteration number
        j=j+1;
    end
    
    % _____________________________________________________________________
    %/
    % displacement and its derivatives after iteration k+1 of increment
    % n+1
    un1k=un+l1*utn*dt+l2*uttn*dt^2+l3*dan*dt^2;
    utn1k=utn+l4*uttn*dt+l5*dan*dt;
    uttn1k=uttn+dan;
    % internal force, stiffness and damping after iteration k+1 of
    % increment n+1
    [FKCn1k,Kn1k,Cn1k,k,d] = BLKIN(un1k,utn1k,k_hi,k_lo,uy,M,ksi,k,d);
    %\_____________________________________________________________________
    
    % assignments to output parameters
    u(n+1)=-un1k;
    ut(n+1)=-utn1k;
    utt(n+1)=-uttn1k+xgtt(n+1);
    Fs(n+1)=FKCn1k;
    Ey(n+1)=-(cumsum(FKCn1k-Cn1k*utn1k)+cumsum(FKCn-Cn*utn))/2.*diff([un1k-un;0]);
    Es(n+1)=cumsum(FKCn1k-Cn1k*utn1k).^2./k_hi/2;
    Ed(n+1)=-(cumsum(Cn1k*utn1k)+cumsum(Cn*utn))/2.*diff([un1k-un;0]);
    jiter(n+1)=j-1;
    % assignments for next increment
    FKCn=FKCn1k;
    Kn=Kn1k;
    Cn=Cn1k;
    un=un1k;
    utn=utn1k;
    uttn=uttn1k;
end

end
