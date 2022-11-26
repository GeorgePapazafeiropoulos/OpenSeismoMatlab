function [PSa,PSv,Sd,Sv,Sa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,ksi,mu,n,tol,...
    pysf,dtTol,AlgID,rinf,maxtol,jmax,dak)
%
% Constant Ductility Response Spectra
%
% [PSA,PSV,SD,SV,SA,FYK,MUK,ITERK]=CDRESP(DT,XGTT,T,KSI,MU,N,TOL,...
%     PYSF,DTTOL,ALGID,RINF,MAXTOL,JMAX,DAK)
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
%         target ductility is achieved.
%     TOL [double(1 x 1)] is the tolerance for convergence for the target
%         ductility.
%     PYSF [double(1 x 1)] is the post-yield stiffness factor, i.e. the
%         ratio of the postyield stiffness to the initial stiffness. PYSF=0
%         is not recommended for simulation of an elastoplastic system. A
%         small positive value is always suggested. PYSF is ignored if
%         MU=1.
%     DTTOL [double(1 x 1)] is the tolerance for resampling of the input
%         acceleration time history. For a given eigenperiod T, resampling
%         takes place if DT/T>dtTol.
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
%     T=(0.04:0.04:4)';
%     ksi=0.05;
%     mu=2;
%     n=50;
%     tol=0.01;
%     pysf=0.1;
%     dtTol=0.02;
%     AlgID='U0-V0-Opt';
%     rinf=1;
%     maxtol=0.01;
%     jmax=200;
%     dak=eps;
%     [CDPSa,CDPSv,CDSd,CDSv,CDSa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,ksi,...
%         mu,n,tol,pysf,dtTol,AlgID,rinf,maxtol,jmax,dak);
%     figure()
%     plot(T,CDSd)
%     figure()
%     plot(T,fyK)
%     figure()
%     plot(T,muK)
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

%% Calculation
% Set integration constants
[w1,w2,w3,W1,W1L1,W2L2,W3L3,W1L4,W2L5,W1L6,l1,l2,l3,l4,l5,rinf] = ...
    TimeIntConstants(AlgID,rinf);
% Initialize
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
% Natural frequencies of SDOFs to be analysed
omega=2*pi./T;
% Flip eigenfrequency vector in order for the half-stepping algorithm
% (HalfStep function) to work from large to small eigenperiods
omega=omega(end:-1:1);
% Set initial conditions
u0=0;
ut0=0;
% Set SDOF mass to unity (any value can be specified without change of the
% results)
m=1;

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
        k_lo,maxuy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
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
        k_lo,minuy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
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
            k_lo,uy1,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
        umax=max(abs(um));
        % ductility factor from eq. 7.2.4. in Chopra
        fy=k_hi*uy1;
        fybark=fy/fpeak;
        mu1=(umax/upeak)/fybark;
        % Step 5.(c) Solve for uy2
        [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,xgtt,m,k_hi,...
            k_lo,uy2,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
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
            k_lo,uy2,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
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

