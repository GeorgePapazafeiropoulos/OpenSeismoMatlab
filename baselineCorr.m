function [cor_xg,cor_xgt,cor_xgtt] = baselineCorr(t,xgtt)
%
% Baseline correction of acceleration time history
%
% [COR_XG,COR_XGT,COR_XGTT] = BASELINECORR(T,XGTT)
%
% Description
%     Linear baseline correction is performed for an uncorrected
%     acceleration time history. Initially, first order fitting (straight
%     line) is performed on the acceleration time history and the fitting
%     line is subrtacted from the acceleration time history, giving thus
%     the first correction. Afterwards, the first correction of the
%     acceleration is integrated to obtain the velocity, and then first
%     order fitting (straight line) is performed on this velocity time
%     history. The gradient of the straight fitting line is then subtracted
%     from the first correction of the acceleration time history, giving
%     thus the second correction of the acceleration time history. The
%     second correction of the acceleration time history is then integrated
%     to give the corrected velocity and displacement time histories.
%
% Input parameters
%     T [double(1:numsteps x 1)] is the time vector of the input
%         acceleration time history XGTT. numsteps is the length of the
%         input acceleration time history.
%     XGTT [double(1:nstep x 1)]: column vector of the acceleration history
%         of the excitation imposed at the base. nstep is the number of
%         time steps of the dynamic response.
%
% Output parameters
%     COR_XG [double(1:nstep x 1)]: time-history of displacement
%     COR_XGT [double(1:nstep x 1)]: time-history of velocity
%     COR_XGTT [double(1:nstep x 1)]: time-history of acceleration
%
% Example
%     fid=fopen('elcentro.dat','r');
%     text=textscan(fid,'%f %f');
%     fclose(fid);
%     time=text{1,1};
%     xgtt1=text{1,2};
%     dt=time(2)-time(1);
%     xgt1 = cumtrapz(time,xgtt1);
%     xg1 = cumtrapz(time,xgt1);
%     [xg2, xgt2, xgtt2] = baselineCorr(time,xgtt1)
%     figure()
%     plot(time,xgtt1)
%     hold on
%     plot(time,xgtt2)
%     figure()
%     plot(time,xgt1)
%     hold on
%     plot(time,xgt2)
%     figure()
%     plot(time,xg1)
%     hold on
%     plot(time,xg2)
%
%__________________________________________________________________________
% Copyright (c) 2018-2022
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D.
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

% degree of fitting
n=1;
% fit a polynomial trend to the acceleration time history
c1 = polyfit(t,xgtt,n);
xgttFit = polyval(c1,t);
% Correct the time histories
cor_xgtt=xgtt-xgttFit;
cor_xgt = cumtrapz(t,cor_xgtt);
cor_xg = cumtrapz(t,cor_xgt);



end