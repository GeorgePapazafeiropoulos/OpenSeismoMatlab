function [U,V,A,f,Es,Ed] = DRHA(k,m,dt,xgtt,varargin)
%
% Dynamic Response History Analysis (DRHA) of a SDOF system
%
% [U,V,A,F,ES,ED] = DRHA(K,M,DT,XGTT,KSI,U0,UT0,RINF)
%
% Description
%     Determine the time history of structural response of a SDOF system
%
% Input parameters
%     K [double(1 x 1)] is the stiffness of the system.
%     M [double(1 x 1)] is the lumped masses of the structure.
%     DT [double(1 x 1)] is the time step of the response history analysis
%         from which the response spectrum is calculated
%     XGTT [double(1:nstep x 1)]: column vector of the acceleration history
%         of the excitation imposed at the base. nstep is the number of
%         time steps of the dynamic response.
%     KSI [double(1 x 1)] is the ratio of critical damping of the SDOF
%         system. Default value 0.05.
%     U0 [double(1 x 1)] is the initial displacement of the SDOF system.
%         Default value 0.
%     UT0 [double(1 x 1)] is the initial velocity of the SDOF system.
%         Default value 0.
%     RINF [double(1 x 1)] is the minimum absolute value of the eigenvalues
%         of the amplification matrix. For the amplification matrix see
%         eq.(61) in Zhou & Tamma (2004). Default value 1.
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
% Copyright (c) 2018-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin<4
    error('Input arguments less than required')
end
if nargin>9
    error('Input arguments more than required')
end
% set defaults for optional inputs
optargs = {0.05,0,0,'U0-V0-CA',0};
% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);
% overwrite the default values by those specified in varargin
optargs(newVals) = varargin(newVals);
% place optional args in memorable variable names
[ksi,u0,ut0,AlgID,rinf] = optargs{:};
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
if ~isscalar(k)
    error('k is not scalar')
end
if k<=0
    error('k is zero or negative')
end
if ~isscalar(m)
    error('m is not scalar')
end
if m<=0
    error('m is zero or negative')
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
Ln=Eigvec'*M;
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

