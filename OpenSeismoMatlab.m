function param=OpenSeismoMatlab(dt,xgtt,sw,varargin)
%
% Seismic parameters and processing of an acceleration time history
%
% Syntax
%     PARAM=OpenSeismoMatlab(DT,XGTT,SW,__)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'PGA')
%     PARAM=OpenSeismoMatlab(DT,XGTT,'PGV')
%     PARAM=OpenSeismoMatlab(DT,XGTT,'PGD')
%     PARAM=OpenSeismoMatlab(DT,XGTT,'ARIAS')
%     PARAM=OpenSeismoMatlab(DT,XGTT,'TIMEHIST',BASELINESW)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'RESAMPLE',DTI)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'ES',T,KSI,ALGID,RINF,DTTOL)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'CDS',T,KSI,MU,PYSF,DTTOL,ALGID,...
%         RINF,MAXTOL,JMAX,DAK)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'FS')
%     PARAM=OpenSeismoMatlab(DT,XGTT,'BUTTERWORTHHIGH',BORDER,FLC)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'BUTTERWORTHLOW',BORDER,FHC)
%     PARAM=OpenSeismoMatlab(DT,XGTT,'IDA',T,LAMBDAF,IM_DM,M,UY,PYSF,...
%         KSI,ALGID,U0,UT0,RINF,MAXTOL,JMAX,DAK)
%     Omit or set as empty ([]) the input arguments for which default
%     values are desired.
%
% Description
%     This function calculates the seismic parameters, develops various
%     spectra and performs various analyses from an acceleration time
%     history. More specifically, it calculates the following:
%     1) Peak ground acceleration
%     2) Peak ground velocity
%     3) Peak ground displacement
%     4) Total cumulative energy and normalized cumulative energy vs time
%     5) Significant duration D_5_95 according to Trifunac & Brady (1975)
%     6) Significant duration D_5_75
%     7) Total Arias intensity (Ia)
%     8) Velocity time history (with baseline correction or not)
%     9) Displacement time history (with baseline correction or not)
%     10) Resampled acceleration time history (i.e. the input acceleration
%        time history with modified time step size)
%     11) Linear elastic pseudo-acceleration response spectrum
%     12) Linear elastic pseudo-velocity response spectrum
%     13) Linear elastic displacement response spectrum
%     14) Linear elastic velocity response spectrum
%     15) Linear elastic acceleration response spectrum
%     16) Constant ductility displacement response spectrum
%     17) Constant ductility velocity response spectrum
%     18) Constant ductility acceleration response spectrum
%     19) Fourier amplitude spectrum
%     20) Mean period (Tm)
%     21) Lowpass Butterworth-filtered acceleration time history
%     22) Highpass Butterworth-filtered acceleration time history
%     23) Incremental Dynamic Analysis (IDA) of SDOF system excited with
%         the input acceleration time history
%     Depending on the value of SW, which determines the type of analysis
%     that OpenSeismoMatlab performs, various additional parameters are
%     needed as input by the user. All possible syntaxes appear above.
%     
% Input parameters
%     DT [double(1 x 1)] is the size of the time step of the input
%         acceleration time history xgtt.
%     XGTT [double(:inf x 1)] is the input acceleration time history.
%     SW [char(1 x :inf)] is a string which determines which parameters,
%         spectras or analyses of the input acceleration time history will
%         be calculated. SW can take one of the following values (strings
%         are case insensitive):
%         'TIMEHIST': the displacement, velocity and acceleration time
%             histories are calculated.
%         'RESAMPLE': the acceleration time history with modified time step
%             size is calculated.
%         'PGA': The peak ground acceleration is calculated.
%         'PGV': The peak ground velocity is calculated.
%         'PGD': The peak ground displacement is calculated.
%         'ARIAS': The total cumulative energy, significant duration
%             D_5_95 according to Trifunac & Brady (1975), significant
%             duration D_5_75 and Arias intensity are calculated.
%         'ES': The linear elastic response spectra and pseudospectra are
%             calculated.
%         'CDS': The constant ductility response spectra are calculated.
%         'FS': The Fourier amplitude spectrum and the mean period are
%             calculated.
%         'BUTTERWORTHHIGH': The high-pass Butterworth filtered
%             acceleration time history is calculated.
%         'BUTTERWORTHLOW': The low-pass Butterworth filtered
%             acceleration time history is calculated.
%         'IDA': Incremental Dynamic Analysis of an elastoplastic SDOF
%             system excited by the input acceleration time history is
%             performed.
% Additional required parameters
%     BASELINESW [logical(1 x 1)] determines if baseline correction will be
%         applied for the calculation of the various time histories.
%     DTI [double(1 x 1)] is the new time step size for resampling of the
%         input acceleration time history.
%     T [double(:inf x 1)] contains the values of eigenperiods for
%         which the response spectra are requested. Its length is the
%         number of SDOF oscillators being analysed to produce the spectra.
%         T must be a vector if SW='ES' or SW='CDS'. T must be scalar if
%         SW='IDA'.
%     KSI [double(1 x 1)] is the fraction of critical viscous damping.
%     ALGID [char(1 x :inf)] is the algorithm to be used for the time
%         integration, if applicable. It can be one of the following
%         strings for superior optimally designed algorithms (strings are
%         case sensitive):
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
%         Default value 'U0-V0-Opt'.
%     MU [double(1 x 1)] is the specified ductility for which the constant
%         ductility response spectra are calculated.
%     PYSF [double(1 x 1)] is the post-yield stiffness factor, i.e. the
%         ratio of the postyield stiffness to the initial stiffness. PYSF=0
%         is not recommended for simulation of an elastoplastic system; a
%         small positive value is always suggested due to numerical
%         reasons. PYSF is ignored if MU=1. Default value 0.01.
%     BORDER [double(1 x 1)] is the order of the Butterworth filter that is
%         applied for filtering of the acceleration time history.
%     FLC [double(1 x 1)] is the low cutoff frequency for the high-pass
%         Butterworth filter.
%     FHC [double(1 x 1)] is the high cutoff frequency for the low-pass
%         Butterworth filter.
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
%     DTTOL [double(1 x 1)] is the tolerance for resampling of the input
%         acceleration time history. For a given eigenperiod T, resampling
%         takes place if DT/T>dtTol. Default value 0.02.
%
% Output parameters
%     PARAM (structure) has the following fields:
%         PARAM.time [double(:inf x 1)] Time
%         PARAM.acc [double(:inf x 1)] Acceleration time history
%         PARAM.vel [double(:inf x 1)] Velocity time history
%         PARAM.disp [double(:inf x 1)] Displacement time history
%         PARAM.PGA [double(1 x 1)] Peak ground acceleration
%         PARAM.PGV [double(1 x 1)] Peak ground velocity
%         PARAM.PGD [double(1 x 1)] Peak ground displacement
%         PARAM.Ecum [double(1 x 1)] Total cumulative energy
%         PARAM.EcumTH [double(:inf x 1)] normalized cumulative
%             energy vs time
%         PARAM.t_5_95 [double(1 x 2)] Time instants at which 5% and 95% of
%             cumulative energy have occurred
%         PARAM.Td_5_95 [double(1 x 1)] Time between when 5% and 95% of
%             cumulative energy has occurred (significant duration
%             according to Trifunac-Brady (1975))
%         PARAM.t_5_75 [double(1 x 2)] Time instants at which 5% and 75% of
%             cumulative energy have occurred
%         PARAM.Td_5_75 [double(1 x 1)] Time between when 5% and 75% of
%             cumulative energy has occurred
%         PARAM.arias [double(1 x 1)] Total Arias intensity (Ia)
%         PARAM.PSa [double(:inf x 1)] Linear elastic pseudo-acceleration
%             response spectrum
%         PARAM.PSv [double(:inf x 1)] Linear elastic pseudo-velocity
%             response spectrum
%         PARAM.Sd [double(:inf x 1)] Linear elastic displacement response
%             spectrum
%         PARAM.Sv [double(:inf x 1)] Linear elastic velocity response
%             spectrum
%         PARAM.Sa [double(:inf x 1)] Linear elastic acceleration response
%             spectrum
%         PARAM.SievABS [double(:inf x 1)] Linear elastic absolute input
%             energy equivalent velocity spectrum
%         PARAM.SievREL [double(:inf x 1)] Linear elastic relative input
%             energy equivalent velocity spectrum
%         PARAM.PredPSa [double(1 x 1)] Predominant acceleration of the PSa
%             spectrum
%         PARAM.PredPeriod [double(1 x 1)] Predominant period of the PSa
%             spectrum
%         PARAM.CDPSa [double(:inf x 1)] Constant ductility
%             pseudo-acceleration response spectrum
%         PARAM.CDPSv [double(:inf x 1)] Constant ductility pseudo-velocity
%             response spectrum
%         PARAM.CDSd [double(:inf x 1)] Constant ductility displacement
%             response spectrum
%         PARAM.CDSv [double(:inf x 1)] Constant ductility velocity
%             response spectrum
%         PARAM.CDSa [double(:inf x 1)] Constant ductility acceleration
%             response spectrum
%         PARAM.fyK [double(:inf x 1)] yield limit that each SDOF must have
%             in order to attain ductility equal to PARAM.muK.
%         PARAM.muK [double(:inf x 1)] achieved ductility for each period
%             (each SDOF).
%         PARAM.iterK [double(:inf x 1)] number of iterations needed for
%             convergence for each period (each SDOF).
%         PARAM.FAS [double(1:2^(nextpow2(length(xgtt))-1) x 1)] Fourier
%             amplitude spectrum
%         PARAM.Tm [double(1 x 1)] Mean period (Tm)
%         PARAM.Fm [double(1 x 1)] Mean frequency (Fm)
%         param.DM [double(:inf x 1)] is the damage measure (DM) of the
%             incremental dynamic analysis
%         param.IM [double(:inf x 1)] is the intensity measure (IM) of the
%             incremental dynamic analysis
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
% Required inputs
if ~isa(dt,'double')
    error('dt is not double');
end
if ~isfinite(dt)
    error('dt is not finite');
end
if ~isscalar(dt)
    error('dt is not scalar')
end
if ~isa(xgtt,'double')
    error('xgtt is not double');
end
if ~isfinite(xgtt)
    error('xgtt is not finite');
end
if ~isvector(xgtt)
    error('xgtt is not vector')
end
if ~ischar(sw)
    error('sw is not a string')
end

%% Calculation
xgtt = xgtt(:);
nxgtt=numel(xgtt);
time=(0:nxgtt-1)'*dt;
sw=lower(sw);
switch sw
    case 'timehist'
        % TIME SERIES
        % Input variables needed:
        % baselineSw
        % Set defaults for optional inputs
        optargs = {true};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [baselineSw] = optargs{:};
        % Check
        if ~islogical(baselineSw)
            error('baselineSw is not true or false')
        end
        if ~isscalar(baselineSw)
            error('baselineSw is not scalar')
        end
        % Engine
        if baselineSw
            [cor_xg,cor_xgt,cor_xgtt] = baselineCorr(time,xgtt);
            param.time=time;
            param.acc=cor_xgtt;
            param.vel=cor_xgt;
            param.disp=cor_xg;
        else
            param.time=time;
            % Acceleration time history
            param.acc = xgtt;
            % Velocity time history
            param.vel = cumtrapz(time,xgtt);
            % Displacement time history
            param.disp = cumtrapz(time,param.vel);
        end
        
    case 'resample'
        % RESAMPLE THE ACCELERATION TIME HISTORY TO INCREASE OR DECREASE
        % THE SIZE OF THE TIME STEP
        % Input variables needed:
        % dti
        % Set defaults for optional inputs
        optargs = {0.01};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [dti] = optargs{:};
        % Check
        if ~isa(dti,'double')
            error('dti is not double');
        end
        if ~isfinite(dti)
            error('dti is not finite');
        end
        if ~isscalar(dti)
            error('dti is not scalar')
        end
        if dti<0
            error('dti is negative')
        end
        % Engine
        [d1,d2] = rat(dt/dti);
        % Resample the acceleration time history
        [xgtt,~] = resample(xgtt,d1,d2);
        NANxgtt=find(isnan(xgtt));
        errxgtt=find(diff(NANxgtt)>1);
        if any(errxgtt)
            error(['Non consecutive NaNs in resampled acceleration ',...
                'time history'])
        end
        if any(NANxgtt)
            xgtt = xgtt(1:NANxgtt(1)-1);
        end
        param.acc = xgtt;
        % Time scale
        param.time=(0:numel(param.acc)-1)'*dti;
        
    case 'pga'
        % PEAK GROUND ACCELERATION
        % Engine
        param.PGA = max(abs(xgtt));
    case 'pgv'
        % PEAK GROUND VELOCITY
        % Engine
        param.PGV = max(abs(cumtrapz(time,xgtt)));
    case 'pgd'
        % PEAK GROUND DISPLACEMENT
        % Engine
        param.PGD = max(abs(cumtrapz(time,cumtrapz(time,xgtt))));
        
    case 'arias'
        param.time=time;
        
        % CUMULATIVE ENERGY
        % Engine
        % time history of cumulative energy
        EcumTH = cumsum(xgtt.^2)*dt;
        % Total cumulative energy at the end of the ground motion
        Ecum = EcumTH(end);
        param.Ecum = Ecum;
        % time history of the normalized cumulative energy
        param.EcumTH = EcumTH/Ecum;
        
        % SIGNIFICANT DURATION
        % Engine
        % elements of the time vector which are within the significant
        % duration 5%-95%
        timed = time(EcumTH>=0.05*Ecum & EcumTH<=0.95*Ecum);
        % starting and ending points of the significant duration
        param.t_5_95 = [timed(1),timed(end)];
        % significant duration
        param.Td_5_95 = timed(end)-timed(1)+dt;
        % elements of the time vector which are within the significant
        % duration 5%-75%
        timed = time(EcumTH>=0.05*Ecum & EcumTH<=0.75*Ecum);
        % starting and ending points of the significant duration
        param.t_5_75 = [timed(1),timed(end)];
        % significant duration
        param.Td_5_75 = timed(end)-timed(1)+dt;
        
        % ARIAS INTENSITY
        % Engine
        % time history of Arias Intensity
        ariasTH = 1/9.81*cumsum(xgtt(time<=param.Td_5_95).^2)*pi*dt/2;
        % Total Arias Intensity at the end of the ground motion
        arias = ariasTH(end);
        param.arias = arias;
        
    case 'es'
        % LINEAR ELASTIC RESPONSE SPECTRA
        % Input variables needed:
        % T
        % ksi
        % AlgID
        % rinf
        % dtTol
        % Set defaults for optional inputs
        optargs = {logspace(log10(0.02),log10(50),1000)',0.05,...
            'U0-V0-Opt',1,0.01};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [T,ksi,AlgID,rinf,dtTol] = optargs{:};
        % Check
        if ~isa(T,'double')
            error('T is not double');
        end
        if ~isfinite(T)
            error('T is not finite');
        end
        if ~isvector(T)
            error('T is not vector')
        end
        
        if ~isa(ksi,'double')
            error('ksi is not double');
        end
        if ~isfinite(ksi)
            error('ksi is not finite');
        end
        if ~isscalar(ksi)
            error('ksi is not scalar')
        end
        if ksi<0
            error('ksi is negative')
        end
        
        if isempty(AlgID)
            error('AlgID is empty')
        end
        if ~ischar(AlgID)
            if ~all(size(AlgID)==[1,14])
                error('AlgID must be a 1x14 vector or string')
            end
        end
        
        if ~isa(rinf,'double')
            error('rinf is not double');
        end
        if ~isfinite(rinf)
            error('rinf is not finite');
        end
        if ~isscalar(rinf)
            error('rinf is not scalar')
        end
        if rinf<0
            error('rinf is negative')
        end
        if rinf>1
            error('rinf is larger than 1')
        end
        
        if ~isa(dtTol,'double')
            error('dtTol is not double');
        end
        if ~isfinite(dtTol)
            error('dtTol is not finite');
        end
        if ~isscalar(dtTol)
            error('dtTol is not scalar')
        end
        if dtTol<0
            error('dtTol is negative')
        end
        
        % Engine
        T=T(:);
        param.Period = T(:);
        [PSa,PSv,Sd,Sv,Sa,SievABS,SievREL]=LEReSp(dt,xgtt,T,ksi,dtTol,...
            AlgID,rinf);
        param.PSa=PSa(:);
        param.PSv=PSv(:);
        param.Sd=Sd(:);
        param.Sv=Sv(:);
        param.Sa=Sa(:);
        param.SievABS=SievABS(:);
        param.SievREL=SievREL(:);
        [a1,a2]=max(PSa(:));
        param.PredPSa=a1;
        param.PredPeriod=T(a2);
        
    case 'cds'
        % CONSTANT DUCTILITY RESPONSE SPECTRA
        % Input variables needed:
        % T
        % ksi
        % mu
        % pysf
        % AlgID
        % rinf
        % maxtol
        % jmax
        % dak
        % dtTol
        % Set defaults for optional inputs
        optargs = {logspace(log10(0.02),log10(50),1000)',0.05,2,...
            0.01,0.01,'U0-V0-Opt',1,0.01,200,eps};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [T,ksi,mu,pysf,dtTol,AlgID,rinf,maxtol,jmax,dak] = optargs{:};
        % Check
        if ~isa(T,'double')
            error('T is not double');
        end
        if ~isfinite(T)
            error('T is not finite');
        end
        if ~isvector(T)
            error('T is not vector')
        end
        
        if ~isa(ksi,'double')
            error('ksi is not double');
        end
        if ~isfinite(ksi)
            error('ksi is not finite');
        end
        if ~isscalar(ksi)
            error('ksi is not scalar')
        end
        if ksi<0
            error('ksi is negative')
        end
        
        if ~isa(mu,'double')
            error('mu is not double');
        end
        if ~isfinite(mu)
            error('mu is not finite');
        end
        if ~isscalar(mu)
            error('mu is not scalar')
        end
        if mu<1
            error('mu is lower than unity')
        end
        
        if ~isa(pysf,'double')
            error('pysf is not double');
        end
        if ~isfinite(pysf)
            error('pysf is not finite');
        end
        if ~isscalar(pysf)
            error('pysf is not scalar')
        end
        if pysf<0
            error('pysf is negative')
        end
        if pysf>1
            error('pysf is larger than 1')
        end
        
        if isempty(AlgID)
            error('AlgID is empty')
        end
        if ~ischar(AlgID)
            if ~all(size(AlgID)==[1,14])
                error('AlgID must be a 1x14 vector or string')
            end
        end
        
        if ~isa(rinf,'double')
            error('rinf is not double');
        end
        if ~isfinite(rinf)
            error('rinf is not finite');
        end
        if ~isscalar(rinf)
            error('rinf is not scalar')
        end
        if rinf<0
            error('rinf is negative')
        end
        if rinf>1
            error('rinf is larger than 1')
        end
        
        if ~isa(maxtol,'double')
            error('maxtol is not double');
        end
        if ~isfinite(maxtol)
            error('maxtol is not finite');
        end
        if ~isscalar(maxtol)
            error('maxtol is not scalar')
        end
        if maxtol<0
            error('maxtol is negative')
        end

        if ~isa(jmax,'double')
            error('jmax is not double');
        end
        if ~isfinite(jmax)
            error('jmax is not finite');
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
        
        if ~isa(dak,'double')
            error('dak is not double');
        end
        if ~isfinite(dak)
            error('dak is not finite');
        end
        if ~isscalar(dak)
            error('dak is not scalar')
        end
        if dak<0
            error('dak is negative')
        end

        if ~isa(dtTol,'double')
            error('dtTol is not double');
        end
        if ~isfinite(dtTol)
            error('dtTol is not finite');
        end
        if ~isscalar(dtTol)
            error('dtTol is not scalar')
        end
        if dtTol<0
            error('dtTol is negative')
        end
        
        % Engine
        T=T(:);
        param.Period = T(:);
        % n [double(1 x 1)] is the maximum number of iterations.
        n=100;
        tol=0.005;
        [CDPSa,CDPSv,CDSd,CDSv,CDSa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,...
            ksi,mu,n,tol,pysf,dtTol,AlgID,rinf,maxtol,jmax,dak);
        param.CDPSa=CDPSa(:);
        param.CDPSv=CDPSv(:);
        param.CDSd=CDSd(:);
        param.CDSv=CDSv(:);
        param.CDSa=CDSa(:);
        param.fyK=fyK(:);
        param.muK=muK(:);
        param.iterK=iterK(:);
        
    case 'ida'
        % INCREMENTAL DYNAMIC ANALYSIS OF SDOF OSCILLATOR
        % Input variables needed:
        % T
        % lambdaF
        % IM_DM
        % m
        % uy
        % pysf
		% ksi
		% AlgID
		% u0
		% ut0
		% rinf
		% maxtol
		% jmax
		% dak
        % Set defaults for optional inputs
        optargs = {1,linspace(0.05,4,20)','Sa_disp',1,0.01,0.01,0.05,...
            'U0-V0-Opt',0,0,1,0.01,200,eps};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [T,lambdaF,IM_DM,m,uy,pysf,ksi,AlgID,u0,ut0,rinf,maxtol,jmax,...
            dak] = optargs{:};
        % Check
        if ~isa(T,'double')
            error('T is not double');
        end
        if ~isfinite(T)
            error('T is not finite');
        end
        if ~isscalar(T)
            error('T is not scalar')
        end
        if T<=0
            error('T must be positive')
        end
        
        if ~isa(lambdaF,'double')
            error('lambdaF is not double');
        end
        if ~isfinite(lambdaF)
            error('lambdaF is not finite');
        end
        if ~isvector(lambdaF)
            error('lambdaF is not vector')
        end
        
        if isempty(IM_DM)
            error('IM_DM is empty')
        end
        if ~ischar(IM_DM)
            error('IM_DM is not of type char')
        end
        if size(IM_DM,1)~=1
            error('IM_DM is not a single row char')
        end
        
        if ~isa(m,'double')
            error('m is not double');
        end
        if ~isfinite(m)
            error('m is not finite');
        end
        if ~isscalar(m)
            error('m is not scalar')
        end
        if m<0
            error('m is negative')
        end
        
        if ~isa(uy,'double')
            error('uy is not double');
        end
        if ~isfinite(uy)
            error('uy is not finite');
        end
        if ~isscalar(uy)
            error('uy is not scalar')
        end
        if uy<0
            error('uy is negative')
        end
        
        if ~isa(pysf,'double')
            error('pysf is not double');
        end
        if ~isfinite(pysf)
            error('pysf is not finite');
        end
        if ~isscalar(pysf)
            error('pysf is not scalar')
        end
        if pysf<0
            error('pysf is negative')
        end
        if pysf>1
            error('pysf is larger than 1')
        end
        
        if ~isa(ksi,'double')
            error('ksi is not double');
        end
        if ~isfinite(ksi)
            error('ksi is not finite');
        end
        if ~isscalar(ksi)
            error('ksi is not scalar')
        end
        if ksi<0
            error('ksi is negative')
        end
        
        if isempty(AlgID)
            error('AlgID is empty')
        end
        if ~ischar(AlgID)
            if ~all(size(AlgID)==[1,14])
                error('AlgID must be a 1x14 vector or string')
            end
        end
        
        if ~isa(u0,'double')
            error('u0 is not double');
        end
        if ~isfinite(u0)
            error('u0 is not finite');
        end
        if ~isscalar(u0)
            error('u0 is not scalar')
        end
        
        if ~isa(ut0,'double')
            error('ut0 is not double');
        end
        if ~isfinite(ut0)
            error('ut0 is not finite');
        end
        if ~isscalar(ut0)
            error('ut0 is not scalar')
        end
        
        if ~isa(rinf,'double')
            error('rinf is not double');
        end
        if ~isfinite(rinf)
            error('rinf is not finite');
        end
        if ~isscalar(rinf)
            error('rinf is not scalar')
        end
        if rinf<0
            error('rinf is negative')
        end
        if rinf>1
            error('rinf is larger than 1')
        end
        
        if ~isa(maxtol,'double')
            error('maxtol is not double');
        end
        if ~isfinite(maxtol)
            error('maxtol is not finite');
        end
        if ~isscalar(maxtol)
            error('maxtol is not scalar')
        end
        if maxtol<0
            error('maxtol is negative')
        end

        if ~isa(jmax,'double')
            error('jmax is not double');
        end
        if ~isfinite(jmax)
            error('jmax is not finite');
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
        
        if ~isa(dak,'double')
            error('dak is not double');
        end
        if ~isfinite(dak)
            error('dak is not finite');
        end
        if ~isscalar(dak)
            error('dak is not scalar')
        end
        if dak<0
            error('dak is negative')
        end

        % Engine
        [DM,IM]=IDA(dt,xgtt,T,lambdaF,IM_DM,m,uy,pysf,ksi,AlgID,u0,ut0,...
            rinf,maxtol,jmax,dak);
        param.DM = DM;
        param.IM = IM;
        
    case 'fs'
        % FOURIER AMPLITUDE SPECTRUM
        % Engine
        [f,U]=FASp(dt,xgtt);
        param.freq = f;
        param.FAS = U;
        
        % MEAN PERIOD AND FREQUENCY
        % Engine
        fi = f(f>0.25 & f<20);
        Ci = U(f>0.25 & f<20);
        Tm = ((Ci(:)'.^2)*(1./fi(:)))/(Ci(:)'*Ci(:));
        param.Tm = Tm;
        Fm = ((Ci(:)'.^2)*(fi(:)))/(Ci(:)'*Ci(:));
        param.Fm = Fm;
        
    case 'butterworthhigh'
        % APPLY HIGH PASS BUTTEWORTH FILTER
        % Input variables needed:
        % bOrder
        % flc
        % Set defaults for optional inputs
        optargs = {4,0.1};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [bOrder,flc] = optargs{:};
        % Check
        if ~isa(bOrder,'double')
            error('bOrder is not double');
        end
        if ~isfinite(bOrder)
            error('bOrder is not finite');
        end
        if bOrder~=floor(bOrder)
            error('bOrder is not integer');
        end
        if ~isa(flc,'double')
            error('flc is not double');
        end
        if ~isfinite(flc)
            error('flc is not finite');
        end
        if ~isscalar(flc)
            error('flc is not scalar')
        end
        if flc<0
            error('flc is negative')
        end
        % Engine
        Nyq=1./(2*dt);
        [b,a] = butter(bOrder,flc/Nyq,'high');
        param.time=time;
        param.acc = filtfilt(b,a,xgtt); % High pass
        
    case 'butterworthlow'
        % APPLY LOW PASS BUTTEWORTH FILTER
        % Input variables needed:
        % bOrder
        % flc
        % Set defaults for optional inputs
        optargs = {4,10};
        nOpts=numel(optargs);
        % Skip any new inputs if they are empty or redundant
        if nargin>3+nOpts
            varargin=varargin(1:nOpts);
        end
        newVals = cellfun(@(x) ~isempty(x), varargin);
        % Overwrite the default values by those specified in varargin
        optargs(newVals) = varargin(newVals);
        % Place optional args in memorable variable names
        [bOrder,fhc] = optargs{:};
        % Check
        if ~isa(bOrder,'double')
            error('bOrder is not double');
        end
        if ~isfinite(bOrder)
            error('bOrder is not finite');
        end
        if bOrder~=floor(bOrder)
            error('bOrder is not integer');
        end
        if ~isa(fhc,'double')
            error('fhc is not double');
        end
        if ~isfinite(fhc)
            error('fhc is not finite');
        end
        if ~isscalar(fhc)
            error('fhc is not scalar')
        end
        if fhc<0
            error('fhc is negative')
        end
        % Engine
        Nyq=1./(2*dt);
        [c,d] = butter(bOrder,fhc/Nyq,'low');
        param.time=time;
        param.acc = filtfilt(c,d,xgtt); % Low pass
        
end
end
