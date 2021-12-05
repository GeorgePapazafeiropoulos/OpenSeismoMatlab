function param=OpenSeismoMatlab(dt,xgtt,varargin)
%
% Seismic parameters of an acceleration time history
%
% PARAM=OPENSEISMOMATLAB(DT,XGTT,SW,BASELINESW,DTI,KSI,T,MU,ALGID)
%
% Description
%     This function calculates the seismic parameters from an acceleration
%     time history. More specifically, it calculates the following:
%     1) Velocity vs time
%     2) Displacement vs time
%     3) Resampled acceleration time history (i.e. the input acceleration
%     time history with modified time step size)
%     4) Peak ground acceleration
%     5) Peak ground velocity
%     6) Peak ground displacement
%     7) Total cumulative energy and normalized cumulative energy vs time
%     8) Significant duration according to Trifunac & Brady (1975)
%     9) Total Arias intensity (Ia)
%     10) Linear elastic pseudo-acceleration response spectrum
%     11) Linear elastic pseudo-velocity response spectrum
%     12) Linear elastic displacement response spectrum
%     13) Linear elastic velocity response spectrum
%     14) Linear elastic acceleration response spectrum
%     15) Constant ductility displacement response spectrum
%     16) Constant ductility velocity response spectrum
%     17) Constant ductility acceleration response spectrum
%     18) Fourier amplitude spectrum
%     19) Mean period (Tm)
%
% Input parameters
%     DT [double(1 x 1)] is the size of the time step of the input
%         acceleration time history xgtt.
%     XGTT [double(1:numsteps x 1)] is the input acceleration time history.
%     SW [char(1 x :inf)] is a string which determines which parameters of
%         the input acceleration time history will be calculated. sw can
%         take one of the following values (strings are case insensitive):
%         'TIMEHIST': the displacement, velocity and acceleration time
%             histories are calculated.
%         'RESAMPLE': the acceleration time history with modified time step
%             size is calculated.
%         'PGA': The peak ground acceleration is calculated.
%         'PGV': The peak ground velocity is calculated.
%         'PGD': The peak ground displacement is calculated.
%         'ARIAS': The total cumulative energy, significant duration
%             according to Trifunac & Brady (1975) and Arias intensity are
%             calculated.
%         'ES': The linear elastic response spectra and pseudospectra are
%             calculated.
%         'CDS': The constant ductility response spectra are calculated.
%         'FAS': The Fourier amplitude spectrum and the mean period are
%             calculated.
%     BASELINESW [logical(1 x 1)] determines if baseline correction will be
%         applied for the calculation of the various output quantities.
%     DTI [double(1 x 1)] is the new time step size for resampling of the
%         input acceleration time history.
%     KSI [double(1 x 1)] is the fraction of critical viscous damping.
%     T [double(1:numSDOFs x 1)] contains the values of eigenperiods for
%         which the response spectra are requested. numSDOFs is the number
%         of SDOF oscillators being analysed to produce the spectra.
%     MU [double(1 x 1)] is the specified ductility for which the response
%         spectra are calculated.
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
%
% Output parameters
%     PARAM (structure) has the following fields:
%         PARAM.vel [double(1:numsteps x 1)] Velocity vs time
%         PARAM.disp [double(1:numsteps x 1)] Displacement vs time
%         PARAM.PGA [double(1 x 1)] Peak ground acceleration
%         PARAM.PGV [double(1 x 1)] Peak ground velocity
%         PARAM.PGD [double(1 x 1)] Peak ground displacement
%         PARAM.Ecum [double(1 x 1)] Total cumulative energy
%         PARAM.EcumTH [double(1:numsteps x 1)] normalized cumulative
%         energy vs time
%         PARAM.t_5_95 [double(1 x 2)] Time instants at which 5% and 95% of
%         cumulative energy have occurred
%         PARAM.Td [double(1 x 1)] Time between when 5% and 95% of
%         cumulative energy has occurred (significant duration according to
%         Trifunac-Brady (1975))
%         PARAM.arias [double(1 x 1)] Total Arias intensity (Ia)
%         PARAM.PSa [double(:inf x 1)] Linear elastic pseudo-acceleration
%         response spectrum
%         PARAM.PSv [double(:inf x 1)] Linear elastic pseudo-velocity
%         response spectrum
%         PARAM.Sd [double(:inf x 1)] Linear elastic displacement response
%         spectrum
%         PARAM.Sv [double(:inf x 1)] Linear elastic velocity response
%         spectrum
%         PARAM.Sa [double(:inf x 1)] Linear elastic acceleration response
%         spectrum
%         PARAM.SievABS [double(:inf x 1)] Linear elastic absolute input
%         energy equivalent velocity spectrum
%         PARAM.SievREL [double(:inf x 1)] Linear elastic relative input
%         energy equivalent velocity spectrum
%         PARAM.PredPSa [double(1 x 1)] Predominant acceleration of the PSa
%         spectrum
%         PARAM.PredPeriod [double(1 x 1)] Predominant period of the PSa
%         spectrum
%         PARAM.CDPSa [double(:inf x 1)] Constant ductility
%         pseudo-acceleration response spectrum
%         PARAM.CDPSv [double(:inf x 1)] Constant ductility pseudo-velocity
%         response spectrum
%         PARAM.CDSd [double(:inf x 1)] Constant ductility displacement
%         response spectrum
%         PARAM.CDSv [double(:inf x 1)] Constant ductility velocity
%         response spectrum
%         PARAM.CDSa [double(:inf x 1)] Constant ductility acceleration
%         response spectrum
%         PARAM.fyK [double(:inf x 1)] yield limit that each SDOF must have
%         in order to attain ductility equal to PARAM.muK.
%         PARAM.muK [double(:inf x 1)] achieved ductility for each period
%         (each SDOF).
%         PARAM.iterK [double(:inf x 1)] number of iterations needed for
%         convergence for each period (each SDOF).
%         PARAM.FAS [double(1:2^(nextpow2(length(xgtt))-1) x 1)] Fourier
%         amplitude spectrum
%         PARAM.Tm [double(1 x 1)] Mean period (Tm)
%         PARAM.Fm [double(1 x 1)] Mean frequency (Fm)
%
%__________________________________________________________________________
% Copyright (c) 2018-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<2
    error('Input arguments less than required')
end
if nargin>9
    error('Input arguments more than required')
end
% set defaults for optional inputs
optargs = {'ES',true,0.01,0.05,logspace(log10(0.02),log10(50),1000)',2,'U0-V0-Opt'};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[sw,baselineSw,dti,ksi,T,mu,AlgID] = optargs{:};
% required inputs
if ~isscalar(dt)
    error('dt is not scalar')
end
if ~isvector(xgtt)
    error('xgtt is not vector')
end
% optional inputs
if ~ischar(sw)
    error('sw is not a string')
end
if ~islogical(baselineSw)
    error('baselineSw is not true or false')
end
if ~isscalar(dti)
    error('dti is not scalar')
end
if dti<0
    error('dti is negative')
end
if ~isscalar(ksi)
    error('ksi is not scalar')
end
if ksi<0
    error('ksi is negative')
end
if ~isvector(T)
    error('time is not vector')
end
if ~isscalar(mu)
    error('mu is not scalar')
end
if mu<1
    error('mu is lower than unity')
end
if ~ischar(AlgID)
    if ~all(size(AlgID)==[1,14])
        error('AlgID must be a 1x14 vector or string')
    end
end

%% Calculation
xgtt = xgtt(:);
nxgtt=numel(xgtt);
time=(0:nxgtt-1)'*dt;
sw=lower(sw);
switch sw
    % TIME SERIES
    case 'timehist'
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
        [d1,d2] = rat(dt/dti);
        % Resample the acceleration time history
        [xgtt,~] = resample(xgtt,d1,d2);
        NANxgtt=find(isnan(xgtt));
        errxgtt=find(diff(NANxgtt)>1);
        if any(errxgtt)
            error('Non consecutive NaNs in resampled acceleration time history')
        end
        if any(NANxgtt)
            xgtt = xgtt(1:NANxgtt(1)-1);
        end
        param.acc = xgtt;
        % Time scale
        param.time=(0:numel(param.acc)-1)'*dti;
        
        % PEAK RESPONSES
    case 'pga'
        % Peak ground acceleration
        param.PGA = max(abs(xgtt));
    case 'pgv'
        % Peak ground velocity
        param.PGV = max(abs(cumtrapz(time,xgtt)));
    case 'pgd'
        % Peak ground displacement
        param.PGD = max(abs(cumtrapz(time,cumtrapz(time,xgtt))));
        
    case 'arias'
        param.time=time;
        
        % CUMULATIVE ENERGY
        % time history of cumulative energy
        EcumTH = cumsum(xgtt.^2)*dt;
        % Total cumulative energy at the end of the ground motion
        Ecum = EcumTH(end);
        param.Ecum = Ecum;
        % time history of the normalized cumulative energy
        param.EcumTH = EcumTH/Ecum;
        
        % SIGNIFICANT DURATION
        % elements of the time vector which are within the significant
        % duration
        timed = time(EcumTH>=0.05*Ecum & EcumTH<=0.95*Ecum);
        % starting and ending points of the significant duration
        param.t_5_95 = [timed(1),timed(end)];
        % significant duration
        param.Td = timed(end)-timed(1)+dt;
        
        % ARIAS INTENSITY
        % time history of Arias Intensity
        ariasTH = 1/9.81*cumsum(xgtt(time<=param.Td).^2)*pi*dt/2;
        % Total Arias Intensity at the end of the ground motion
        arias = ariasTH(end);
        param.arias = arias;
        
    case 'es'
        % LINEAR ELASTIC RESPONSE SPECTRA
        T=T(:);
        param.Period = T(:);
        dtTol=0.02;
        rinf=1;
        [PSa,PSv,Sd,Sv,Sa,SievABS,SievREL]=LEReSp(dt,xgtt,T,ksi,dtTol,AlgID,rinf);
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
        T=T(:);
        param.Period = T(:);
        % n [double(1 x 1)] is the maximum number of iterations.
        n=65;
        tol=0.005;
        dtTol=0.02;
        rinf=1;
        [CDPSa,CDPSv,CDSd,CDSv,CDSa,fyK,muK,iterK]=CDReSp(dt,xgtt,T,ksi,mu,n,tol,...
            dtTol,AlgID,rinf);
        param.CDPSa=CDPSa(:);
        param.CDPSv=CDPSv(:);
        param.CDSd=CDSd(:);
        param.CDSv=CDSv(:);
        param.CDSa=CDSa(:);
        param.fyK=fyK(:);
        param.muK=muK(:);
        param.iterK=iterK(:);
        
    case 'fas'
        % FOURIER AMPLITUDE SPECTRUM
        [f,U]=FASp(dt,xgtt);
        param.freq = f;
        param.FAS = U;
        
        % MEAN PERIOD AND FREQUENCY
        fi = f(f>0.25 & f<20);
        Ci = U(f>0.25 & f<20);
        Tm = ((Ci(:)'.^2)*(1./fi(:)))/(Ci(:)'*Ci(:));
        param.Tm = Tm;
        Fm = ((Ci(:)'.^2)*(fi(:)))/(Ci(:)'*Ci(:));
        param.Fm = Fm;
        
end
end
