function [U,V,A,f,Es,Ed] = DRHA(k,m,dt,xgtt,ksi,u0,ut0,AlgID,rinf)
%
% Dynamic Response History Analysis
%
% [U,V,A,F,ES,ED] = DRHA(K,M,DT,XGTT,KSI,U0,UT0,ALGID,RINF)
%
% Description
%     Calculate the dynamic response of a linear MDOF system using modal
%     analysis. This function is part of the OpenSeismoMatlab software. It
%     can be used as standalone, however attention is needed for the
%     correctness of the input arguments, since no checks are performed in
%     this function. See the example example_DRHA.m for more details about
%     how this function can be implemented.
%
% Input parameters
%     K [double(:inf x 1)] is the stiffness of the system.
%     M [double(:inf x 1)] is the lumped masses of the structure.
%     DT [double(1 x 1)] is the time step of the dynamic response history
%         analysis
%     XGTT [double(1:nstep x 1)]: column vector of the acceleration history
%         of the excitation imposed at the base. nstep is the number of
%         time steps of the dynamic response.
%     KSI [double(1 x 1)] is the ratio of critical damping of the SDOF
%         system.
%     U0 [double(:inf x 1)] is the initial displacement of the SDOF system.
%     UT0 [double(:inf x 1)] is the initial velocity of the SDOF system.
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
%     U [double(1 x 1:nstep)]: displacement time history.
%     V [double(1 x 1:nstep)]: velocity time history.
%     A [double(1 x 1:nstep)]: acceleration time history.
%     F [double(1 x 1:nstep)]: equivalent static force time history.
%     ES [double(1 x 1:nstep)]: time-history of the recoverable
%         strain energy of the system (total and not incremental).
%     ED [double(1 x 1:nstep)]: time-history of the energy
%         dissipated by viscoelastic damping during each time step
%         (incremental). cumsum(Ed) gives the time history of the total
%         energy dissipated at dof i from the start of the dynamic
%         analysis.
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________


%% Calculation
% Number of time integration steps
nstep=length(xgtt);
% construct lumped mass matrix
M=diag(m,0);
% Assemble stiffness matrix
ndofs=size(M,1);
K=zeros(ndofs+1);
K=K-diag(k,-1)-diag(k,1);
K(1:end-1,1:end-1)=K(1:end-1,1:end-1)+diag(k,0);
K(2:end,2:end)=K(2:end,2:end)+diag(k,0);
K(end,:)=[];
K(:,end)=[];
% Calculate eigenvalues, eigenvectors and their number
[Eigvec,Eigval]=eig(K,M);
% Assemble damping matrix
C=zeros(ndofs);
for i=1:ndofs
    c=2*ksi*sqrt(Eigval(i,i))*(((M*Eigvec(:,i))*Eigvec(:,i)')*M);
    C=C+c;
end
% Take the eigenvalues in column vector and sort in ascending order of
% eigenfrequency:
D1=diag(Eigval,0);
% Generalized masses Mn for all eigenmodes from eq.(13.1.5) of Chopra
% (2012).
Mn=diag(Eigvec'*M*Eigvec);
% Ln coefficients from eq.(13.1.5) of Chopra (2012).
Ln=Eigvec'*M*ones(ndofs,1);
% Gamman coefficients from eq.(13.1.5) of Chopra (2012).
Gamman=Ln./Mn;
% Eigenperiods of the building
omega=D1.^0.5;
% Initial displacements
u0mod=Eigvec\u0;
% Normalization
u0mod=u0mod./Mn;
% Initial velocities
ut0mod=Eigvec\ut0;
% Normalization
ut0mod=ut0mod./Mn;
% Displacements, velocities and accelerations of the response of the
% eigenmodes of the structure for the given earthquake
neig=size(Eigvec,1);
U=zeros(neig,nstep);
V=zeros(neig,nstep);
A=zeros(neig,nstep);
f=zeros(neig,nstep);
for i=1:neig
    [u,ut,utt] = LIDA(dt,xgtt,omega(i),ksi,u0mod(i),ut0mod(i),AlgID,rinf);
    U=U+Gamman(i)*Eigvec(:,i)*u';
    V=V+Gamman(i)*Eigvec(:,i)*ut';
    A=A+Gamman(i)*Eigvec(:,i)*utt';
    f=f+Gamman(i)*omega(i)^2*(M*Eigvec(:,i))*u';
end
Es=cumsum(K*U).^2./k(:,ones(1,nstep))/2;
Ed=cumsum(C*V).*diff([-V;zeros(1,nstep)])*dt;
end

