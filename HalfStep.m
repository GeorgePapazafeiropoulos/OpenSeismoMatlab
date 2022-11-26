function uNew = HalfStep(u)
%
% Reproduce signal with half time step
%
% UNEW = HALFSTEP(U)
%
% Input parameters
%     U [double(1:n x 1)] is the input signal with time step dt.
%
% Output parameters
%     UNEW [double(1:n x 1)] is the output signal with time step dt/2.
%
% Verification:
%     u=0.2:0.2:4;
%     uNew=HalfStep(u);
%     figure()
%     plot((1:numel(u)),u)
%     hold on
%     plot((1:0.5:numel(u)),uNew)
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

%% Initial checks
if nargin>1
    error('Input arguments more than required')
end
% required inputs
if ~isvector(u)
    error('u is not vector')
end
t=false;
if size(u,2)~=1
    u=u(:);
    t=true;
end
%% Calculation
a=[([0;u(1:end-1)]+u)/2,u]';
uNew=a(:);
uNew(1)=[];
if t
    uNew=uNew';
end

end