function [DM,IM]=IDA(dt,xgtt,T,lambdaF,IM_DM,m,uy,pysf,ksi,AlgID,u0,ut0,...
    rinf,maxtol,jmax,dak)
%
% Incremental Dynamic Analysis
%
% [DM,IM]=IDA(DT,XGTT,T,LAMBDAF,IM_DM,M,UY,PYSF,KSI,ALGID,U0,UT0,...
%     RINF,MAXTOL,JMAX,DAK)
%
% Description
%     This function performs incremental dynamic analysis of a given
%     acceleration time history and SDOF oscillator.
%
% Input parameters
%     DT [double(1 x 1)] is the time step of the input acceleration time
%         history XGTT.
%     XGTT [double(:inf x 1)] is the input acceleration time history.
%         numsteps is the length of the input acceleration time history.
%     T [double(1 x 1)] contains the eigenperiod of the SDOF system for
%         which the incremental dynamic analysis response curve is
%         requested.
%     LAMBDAF [double(:inf x 1)] contains the values of the scaling factor
%         (lambda factor) for the incremental dynamic analysis. 
%     IM_DM [char(1 x :inf)] is the Intensity Measure (IM) - Damage Measure
%         (DM) pair that is to be calculated from the incremental dynamic
%         analysis. IM_DM can take one of the following values (strings are
%         case insensitive):
%         'SA_MU': Spectral acceleration-ductility
%         'PGD_MU': Peak displacement-ductility
%         'PGV_MU': Peak velocity-ductility
%         'PGA_MU': Peak acceleration-ductility
%         'SA_DISP': Spectral acceleration-displacement
%         'PGD_DISP': Peak displacement-displacement
%         'PGV_DISP': Peak velocity-displacement
%         'PGA_DISP': Peak acceleration-displacement
%         'SA_VEL': Spectral acceleration-velocity
%         'PGD_VEL': Peak displacement-velocity
%         'PGV_VEL': Peak velocity-velocity
%         'PGA_VEL': Peak acceleration-velocity
%         'SA_ACC': Spectral acceleration-acceleration
%         'PGD_ACC': Peak displacement-acceleration
%         'PGV_ACC': Peak velocity-acceleration
%         'PGA_ACC': Peak acceleration-acceleration
%     M [double(1 x 1)] is the mass of the SDOF oscillator.
%     UY [double(1 x 1)] is the yield displacement of the SDOF oscillator.
%     PYSF [double(1 x 1)] is the post-yield stiffness factor, i.e. the
%         ratio of the postyield stiffness to the initial stiffness. PYSF=0
%         is not recommended for simulation of an elastoplastic system. A
%         small positive value is always suggested. PYSF is ignored if
%         MU=1.
%     KSI [double(1 x 1)] is the fraction of critical viscous damping.
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
%     U0 [double(1 x 1)] is the initial displacement of the SDOF
%         oscillator.
%     UT0 [double(1 x 1)] is the initial velocity of the SDOF oscillator.
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004). Default value 0.
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
%     DM [double(:inf x 1)] is the Damage Measure.
%     IM [double(:inf x 1)] is the Intensity Measure.
%
% Example
%     eqmotions={'elcentro'};
%     data=load([eqmotions{1},'.dat']);
%     t=data(:,1);
%     dt=t(2)-t(1);
%     xgtt=data(:,2);
%     sw='ida';
%     T=1;
%     lambdaF=logspace(log10(0.001),log10(10),100);
%     IM_DM='Sa_disp';
%     m=1;
%     uy = 0.082*9.81/(2*pi/T)^2;
%     pysf=0.01;
%     ksi=0.05;
%     S5=OpenSeismoMatlab(dt,xgtt,sw,T,lambdaF,IM_DM,m,uy,pysf,ksi);
%     figure()
%     plot(S5.DM*1000,S5.IM/9.81,'k','LineWidth',1)
%     grid on
%     xlabel('Displacement (mm)')
%     ylabel('Sa(T1,5%)[g]')
%     xlim([0,200])
%     ylim([0,0.7])
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________


%% Calculation
% Initial and postyield stiffnesses
k_hi=m*(2*pi/T)^2;
k_lo=k_hi*pysf;

% Initialize
k=length(lambdaF);
IM_DM=lower(IM_DM);
switch IM_DM
    case 'sa_mu'
        S1=OpenSeismoMatlab(dt,xgtt,'es',T,ksi,AlgID); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.Sa; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um))/uy; %%
        end
    case 'pgd_mu'
        S1=OpenSeismoMatlab(dt,xgtt,'pgd'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGD; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um))/uy; %%
        end
    case 'pgv_mu'
        S1=OpenSeismoMatlab(dt,xgtt,'pgv'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGV; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um))/uy; %%
        end
    case 'pga_mu'
        S1=OpenSeismoMatlab(dt,xgtt,'pga'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGA; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um))/uy; %%
        end
    case 'sa_disp'
        S1=OpenSeismoMatlab(dt,xgtt,'es',T,ksi,AlgID); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.Sa; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um)); %%
        end
    case 'pgd_disp'
        S1=OpenSeismoMatlab(dt,xgtt,'pgd'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGD; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um)); %%
        end
    case 'pgv_disp'
        S1=OpenSeismoMatlab(dt,xgtt,'pgv'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGV; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um)); %%
        end
    case 'pga_disp'
        S1=OpenSeismoMatlab(dt,xgtt,'pga'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGA; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(um)); %%
        end
    case 'sa_vel'
        S1=OpenSeismoMatlab(dt,xgtt,'es',T,ksi,AlgID); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.Sa; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umt)); %%
        end
    case 'pgd_vel'
        S1=OpenSeismoMatlab(dt,xgtt,'pgd'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGD; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umt)); %%
        end
    case 'pgv_vel'
        S1=OpenSeismoMatlab(dt,xgtt,'pgv'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGV; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umt)); %%
        end
    case 'pga_vel'
        S1=OpenSeismoMatlab(dt,xgtt,'pga'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGA; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umt)); %%
        end
    case 'sa_acc'
        S1=OpenSeismoMatlab(dt,xgtt,'es',T,ksi,AlgID); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.Sa; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umtt)); %%
        end
    case 'pgd_acc'
        S1=OpenSeismoMatlab(dt,xgtt,'pgd'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGD; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umtt)); %%
        end
    case 'pgv_acc'
        S1=OpenSeismoMatlab(dt,xgtt,'pgv'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGV; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umtt)); %%
        end
    case 'pga_acc'
        S1=OpenSeismoMatlab(dt,xgtt,'pga'); %
        IM=zeros(length(lambdaF),1);
        DM=zeros(length(lambdaF),1);
        for j=1:k
            IM(j)=lambdaF(j)*S1.PGA; %
            [um,umt,umtt,p,Ey,Es,Ed,iter] = NLIDABLKIN(dt,...
                lambdaF(j)*xgtt,m,k_hi,...
                k_lo,uy,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,dak);
            DM(j)=max(abs(umtt)); %%
        end
    otherwise
        error('No appropriate IM_DM specified.');
end

end
