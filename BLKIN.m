function [f,K,C,k_status,d] = BLKIN(u,ut,k_hi,k_lo,uy,M,ksi,k_status,d)
%
% Bilinear elastoplastic hysteretic model with elastic viscous damping
%
% [F,K,C,K_STATUS,D] = BLKIN(U,UT,K_HI,K_LO,UY,M,KSI,K_STATUS,D)
%
% Description
%     Define the internal force vector, tangent stiffness matrix and
%     tangent damping matrix of a bilinear elastoplastic hysteretic
%     structure with elastic damping as a function of displacement and
%     velocity.
%     The MDOF structure modeled with this function consists of lumped
%     masses connected with stiffness and damping elements in series. Each
%     lumped mass has one degree of freedom. The first degree of freedom is
%     at the top of the structure and the last at its fixed base. However,
%     the last degree of freedom is not included in the input arguments of
%     the function, i.e. not contained in ndof, as it is always fixed.
%     The nonlinear stiffness is virtually of the bilinear type, where an
%     initial stiffness and a post-yield stiffness are defined. The
%     unloading or reloading curve of this model are parallel to the
%     initial loading curve, and a hysteresis loop is created by
%     continuously loading and unloading the structure above its yield
%     limit. This behavior can be viewed as hardening of the kinematic
%     type.
%     An appropriate reference for this function definition is Hughes,
%     Pister & Taylor (1979): "Implicit-explicit finite elements in
%     nonlinear transient analysis". This function should be defined in
%     accordance with equations (3.1), (3.2) and (3.3) of this paper. This
%     representation has as special cases nonlinear elasticity and a class
%     of nonlinear “rate-type” viscoelastic materials. Tangent stiffness
%     and tangent damping matrices are the "consistent" linearized
%     operators associated to f in the sense of [Hughes & Pister,
%     "Consistent linearization in mechanics of solids", Computers and
%     Structures, 8 (1978) 391-397].
%
% Input parameters
%     U [double(1 x 1)] is the absolute displacement.
%     UT [double(1 x 1)] is the absolute velocity.
%     K_HI [double(1 x 1)] is the initial stiffness of the system before
%         its first yield, i.e. the high stiffness.
%     K_LO [double(1 x 1)] is the post-yield stiffness of the system, i.e.
%         the low stiffness.
%     UY [double(1 x 1)] is the yield limit of the structure. The structure
%         is considered to yield, if the displacement exceeds uy(i).
%     M [double(1 x 1)] is the lumped mass.
%     KSI [double(1 x 1)] is the ratio of critical viscous damping of the
%         system, assumed to be unique for all damping elements of the
%         structure.
%     K_STATUS [double(1 x 1)] is the is the stiffness vector which takes
%         into account any plastic response of the structure. It is used to
%         record the status of the structure so that it is known before the
%         next application of this function at a next (time) step.
%         Initialize by setting K_STATUS=K_HI.
%     D [double(1 x 1)] is the is the equilibrium displacement vector which
%         takes into account any plastic response of the structure. It is
%         used to record the status of the structure so that it is known
%         before the next application of this function at a next (time)
%         step. Initialize by setting D=zeros(ndof,1).
%
% Output parameters
%     F [double(1 x 1)] is the internal force vector of the structure (sum
%         of forces due to stiffness and damping) at displacement u and
%         velocity ut
%     K [double(1 x 1)] is the tangent stiffness matrix (nonlinear function
%         of displacement u and velocity ut). It is equivalent to the
%         derivative d(f)/d(u)
%     C [double(1 x 1)] is the tangent damping matrix (nonlinear function
%         of displacement u and velocity ut). It is equivalent to the
%         derivative d(f)/d(u)
%     K_STATUS [double(1 x 1)] is the is the stiffness vector which takes
%         into account any plastic response of the structure. It is used to
%         record the status of the structure so that it is known before the
%         next application of this function at a next (time) step.
%     D [double(1 x 1)] is the is the equilibrium displacement vector which
%         takes into account any plastic response of the structure. It is
%         used to record the status of the structure so that it is known
%         before the next application of this function at a next (time)
%         step.
%
% Verification:
%     u=0:0.2:4;
%     ut=0.001*ones(1,numel(u));
%     u=[u,u(end:-1:1)];
%     ut=[ut,-ut];
%     u=[u,-u];
%     ut=[ut,ut(end:-1:1)];
%     u=[u u];
%     ut=[ut ut];
%     k_hi=1000;
%     k_lo=1;
%     uy=2;
%     M=1;
%     ksi=0.05;
%     k=k_hi;
%     d=0;
%     f=zeros(1,numel(u));
%     for i=1:numel(u)
%         [f(i),K,C,k,d] = BLKIN(u(i),ut(i),k_hi,k_lo,uy,M,ksi,k,d);
%     end
%     figure()
%     plot(u,f)
%
%__________________________________________________________________________
% Copyright (c) 2018-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%%
% Elastic tangent stiffness matrix
K=k_hi;
% Elastic tangent damping matrix
C=2*ksi*sqrt(K*M);
% force from stiffness (not damping) of the current storey
fK=k_status*(u(1)-d);
% eq.(46) in ...
fy=k_lo*(u(1))+(k_hi-k_lo)*(uy.*sign(ut(1)));
% check for yielding or load reversal
if k_status==k_hi && ut(1)>0 && fK>fy
    % check for yielding
    % the system has just exceeded its positive yield force level
    k_status=k_lo;
    d=(1-k_hi/k_lo)*uy;
elseif k_status==k_hi && ut(1)<0 && fK<fy
    % check for yielding
    % the system has just exceeded its negative yield force level
    k_status=k_lo;
    d=(k_hi/k_lo-1)*uy;
elseif k_status==k_lo && fK*(ut(1))<0
    % check for load reversal
    % the system reloads from negative ultimate displacement or unloads
    % from positive ultimate displacement
    k_status=k_hi;
    d=(u(1))-k_lo/k_hi*(u(1)-d);
end
fK_bak=k_status*(u(1)-d);
% Update the elastic tangent stiffness matrix
K=k_status;
% internal force due to stiffness and damping
f = fK_bak + C*ut;

end