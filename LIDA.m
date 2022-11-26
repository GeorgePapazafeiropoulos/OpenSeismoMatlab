function [u,ut,utt] = LIDA(dt,xgtt,omega,ksi,u0,ut0,AlgID,rinf)
%
% Linear Implicit Dynamic Analysis
%
% [U,UT,UTT] = LIDA(DT,XGTT,OMEGA,KSI,U0,UT0,ALGID,RINF)
%
% Description
%     Linear implicit direct time integration of second order differential
%     equation of motion of dynamic response of linear elastic SDOF systems
%     The General Single Step Single Solve (GSSSS) family of algorithms
%     published by X.Zhou & K.K.Tamma (2004) is employed for direct time
%     integration of the general linear or nonlinear structural Single
%     Degree of Freedom (SDOF) dynamic problem. The optimal numerical
%     dissipation and dispersion zero order displacement zero order
%     velocity algorithm designed according to the above journal article,
%     is used in this routine. This algorithm encompasses the scope of
%     Linear Multi-Step (LMS) methods and is limited by the Dahlquist
%     barrier theorem (Dahlquist,1963). The force - displacement - velocity
%     relation of the SDOF structure is linear. This function is part of
%     the OpenSeismoMatlab software. It can be used as standalone, however
%     attention is needed for the correctness of the input arguments, since
%     no checks are performed in this function. See the example
%     example_LIDA.m for more details about how this function can be
%     implemented.
%
% Input parameters
%     DT [double(1 x 1)] is the time step
%     XGTT [double(1:nstep x 1)] is the column vector of the acceleration
%         history of the excitation imposed at the base. nstep is the
%         number of time steps of the dynamic response.
%     OMEGA [double(1 x 1)] is the eigenfrequency of the structure in
%         rad/sec.
%     KSI [double(1 x 1)] is the ratio of critical damping of the SDOF
%         system.
%     U0 [double(1 x 1)] is the initial displacement of the SDOF system.
%     UT0 [double(1 x 1)] is the initial velocity of the SDOF system.
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
%     U [double(1:nstep x 1)] is the time-history of displacement
%     UT [double(1:nstep x 1)] is the time-history of velocity
%     UTT [double(1:nstep x 1)] is the time-history of acceleration
%
% Example (Figure 6.6.1 in Chopra, Tn=1sec)
%     dt=0.02;
%     fid=fopen('elcentro.dat','r');
%     text=textscan(fid,'%f %f');
%     fclose(fid);
%     xgtt=text{1,2};
%     Tn=1;
%     omega=2*pi/Tn;
%     ksi=0.02;
%     u0=0;
%     ut0=0;
%     AlgID='U0-V0-Opt';
%     rinf=1;
%     [u,ut,utt] = LIDA(dt,xgtt,omega,ksi,u0,ut0,AlgID,rinf);
%     D=max(abs(u))/0.0254
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
% Transfer function denominator
Omega=omega*dt;
D=W1L6+2.*W2L5.*ksi.*Omega+W3L3.*Omega.^2;
A31=-Omega.^2./D;
A32=-1./D.*(2.*ksi.*Omega+W1L1.*Omega.^2);
A33=1-1./D.*(1+2.*W1L4.*ksi.*Omega+W2L2.*Omega.^2);
A11=1+l3.*A31;
A12=l1+l3.*A32;
A13=l2-l3.*(1-A33);
A21=l5.*A31;
A22=1+l5.*A32;
A23=l4-l5.*(1-A33);
% Amplification matrix
A=[A11 A12 A13;A21 A22 A23;A31 A32 A33];
% Amplification matrix invariants
A1=A(1,1)+A(2,2)+A(3,3);
A2=A(1,1)*A(2,2)-A(1,2)*A(2,1)+A(1,1)*A(3,3)-A(1,3)*A(3,1)+A(2,2)*A(3,3)-...
    A(2,3)*A(3,2);
A3=A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)+A(1,2)*...
    A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2)-A(1,3)*A(2,2)*A(3,1);
% Transfer function denominator
a=[1 -A1 A2 -A3];
% Transfer function nominator
B1=1./D.*dt^2.*l3.*W1;
B2=1./D.*dt^2.*(l3.*(1-W1)-(A22+A33).*l3.*W1+A12.*l5.*W1+A13.*W1);
B3=1./D.*dt^2.*(-(A22+A33).*l3.*(1-W1)+A12.*l5.*(1-W1)+A13.*(1-W1)+...
    (A22.*A33-A23.*A32).*l3.*W1-(A12.*A33-A13.*A32).*l5.*W1+(A12.*A23-...
    A13.*A22).*W1);
B4=1./D.*dt^2.*((A22.*A33-A23.*A32).*l3.*(1-W1)-(A12.*A33-A13.*A32).*l5.*(1-...
    W1)+(A12.*A23-A13.*A22).*(1-W1));
b=[B1,B2,B3,B4];
% form initial conditions for filter function
% equivalent external force
f=-xgtt;
% stiffness
k=omega.^2;
% damping constants
c=2.*omega.*ksi;
% initial acceleration
utt0=-f(1)-(k*u0+c*ut0);
U_1=A\[u0;dt*ut0;dt^2*utt0];
u_1=U_1(1);
U_2=A\U_1;
u_2=U_2(1);
ypast=[u0,u_1,u_2];
vinit=zeros(1,3);
vinit(3:-1:1) = filter(-a(4:-1:2),1,ypast);
% main dynamic analysis
u=filter(b,a,f,vinit);
% calculate velocity from the following system of equations:
% 1st: the first scalar equation of the matrix equation (60) in X.Zhou &
% K.K.Tamma (2004)
% 2nd: equation of motion (eq.6.12.3 in Chopra:Dynamics of Structures)
C_u=omega^2*A(1,3)*dt^2-A(1,1);
C_f=-A(1,3)*dt^2;
C_ut=A(1,2)*dt-A(1,3)*dt^2*2*ksi*omega;
L=1/D*l3*dt^2*((1-W1)*[0;f(1:end-1)]+W1*f);
ut=(u+C_u*[u0;u(1:end-1)]+C_f*[0;f(1:end-1)]-L)/C_ut;
% calculate acceleration from equation of motion
utt=-omega^2*u-2*ksi*omega*ut;
end

