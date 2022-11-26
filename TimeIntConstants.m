function [w1,w2,w3,W1,W1L1,W2L2,W3L3,W1L4,W2L5,W1L6,l1,l2,l3,l4,l5,rinf] = ...
    TimeIntConstants(AlgID,rinf)
% Set integration constants
if all(size(AlgID)==[1,14])
    % define integration constants explicitly
    w1=AlgID(1);
    w2=AlgID(2);
    w3=AlgID(3);
    W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
    % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
    % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
    W1L1=AlgID(4);
    W2L2=AlgID(5);
    W3L3=AlgID(6);
    W1L4=AlgID(7);
    W2L5=AlgID(8);
    W1L6=AlgID(9);
    l1=AlgID(10);
    l2=AlgID(11);
    l3=AlgID(12);
    l4=AlgID(13);
    l5=AlgID(14);
else
    switch AlgID
        case 'U0-V0-Opt'
            % zero-order displacement & velocity overshooting behavior and
            % optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf); % suggested
            w2=15*(3-4*rinf)/(1-4*rinf); % suggested
            w3=-35*(1-rinf)/(1-4*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/2/(1+rinf)^2;
            W1L4=1/(1+rinf);
            W2L5=1/(1+rinf)^2; % suggested
            W1L6=(3-rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-CA'
            % zero-order displacement & velocity overshooting behavior and
            % continuous acceleration
            % rinf must belong to [1/3 1]
            if rinf<1/3
                rinf=1/3;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/3');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-5*rinf)/(3-7*rinf); % suggested
            w2=15*(1-13*rinf)/(3-7*rinf); % suggested
            w3=140*rinf/(3-7*rinf); % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=(1+3*rinf)/4/(1+rinf);
            W3L3=(1+3*rinf)/4/(1+rinf)^2;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=(1+3*rinf)/2/(1+rinf)^2; % suggested
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V0-DA'
            % zero-order displacement & velocity overshooting behavior and
            % discontinuous acceleration
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15; % suggested
            w2=45; % suggested
            w3=-35; % suggested
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/2/(1+rinf);
            W1L4=1;
            W2L5=1/(1+rinf); % suggested
            W1L6=(3+rinf)/2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U0-V1-Opt'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'generalized a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % This is the generalized a-method (Chung & Hulbert, 1993)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(1-4*rinf);
            w2=15*(3-4*rinf)/(1-4*rinf);
            w3=-35*(1-rinf)/(1-4*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1/(1+rinf);
            W2L2=1/2/(1+rinf);
            W3L3=1/(1+rinf)^3;
            W1L4=1/(1+rinf);
            W2L5=(3-rinf)/2/(1+rinf)^2;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-CA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'HHT a-method'
            % zero-order displacement & first-order velocity overshooting
            % behavior and continuous acceleration
            % This is the Hilber-Hughes-Taylor method (Hilber, Hughes &
            % Taylor, 1977)
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15*(1-2*rinf)/(2-3*rinf);
            w2=15*(2-5*rinf)/(2-3*rinf);
            w3=-35*(1-3*rinf)/2/(2-3*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=2*rinf/(1+rinf);
            W2L2=rinf/(1+rinf);
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=2*rinf/(1+rinf);
            W2L5=rinf*(3-rinf)/(1+rinf)^2;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U0-V1-DA'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'WBZ'
            % zero-order displacement & first-order velocity overshooting
            % behavior and discontinuous acceleration
            % This is the Wood–Bossak–Zienkiewicz method (Wood, Bossak &
            % Zienkiewicz, 1980)
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/2;
            W3L3=1/(1+rinf)^2;
            W1L4=1;
            W2L5=(3-rinf)/2/(1+rinf);
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'U1-V0-Opt'
            % first-order displacement & zero-order velocity overshooting
            % behavior and optimal numerical dissipation and dispersion
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1; % mid-point rule a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-8*rinf+6*rinf^2)/(9-22*rinf+19*rinf^2);
            w2=15*(25-74*rinf+53*rinf^2)/2/(9-22*rinf+19*rinf^2);
            w3=-35*(3-10*rinf+7*rinf^2)/(9-22*rinf+19*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3-rinf)/2/(1+rinf);
            W2L2=1/(1+rinf)^2;
            W3L3=1/(1+rinf)^3;
            W1L4=(3-rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^3;
            W1L6=(2-rinf)/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-CA'
            % first-order displacement & zero-order velocity overshooting
            % behavior and continuous acceleration
            % rinf must belong to [1/2 1]
            if rinf<1/2
                rinf=1/2;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1/2');
            end
            if rinf>1
                rinf=1; % Newmark average acceleration a-form algorithm
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-60*(2-8*rinf+7*rinf^2)/(11-48*rinf+41*rinf^2);
            w2=15*(37-140*rinf+127*rinf^2)/2/(11-48*rinf+41*rinf^2);
            w3=-35*(5-18*rinf+17*rinf^2)/(11-48*rinf+41*rinf^2);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(1+3*rinf)/2/(1+rinf);
            W2L2=2*rinf/(1+rinf)^2;
            W3L3=2*rinf/(1+rinf)^3;
            W1L4=(1+3*rinf)/2/(1+rinf);
            W2L5=4*rinf/(1+rinf)^3;
            W1L6=1;
            l1=1;
            l2=1/2;
            l3=1/2/(1+rinf);
            l4=1;
            l5=1/(1+rinf);
        case 'U1-V0-DA'
            % first-order displacement & zero-order velocity overshooting behavior
            % and discontinuous acceleration
            % This is the Newmark average acceleration a-form algorithm
            % rinf must belong to [0 1]
            if rinf<0
                rinf=0;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 0');
            end
            if rinf>1
                rinf=1;
                warning('Minimum absolute eigenvalue of amplification matrix is set to 1');
            end
            w1=-30*(3-4*rinf)/(9-11*rinf);
            w2=15*(25-37*rinf)/2/(9-11*rinf);
            w3=-35*(3-5*rinf)/(9-11*rinf);
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=(3+rinf)/2/(1+rinf);
            W2L2=1/(1+rinf);
            W3L3=1/(1+rinf)^2;
            W1L4=(3+rinf)/2/(1+rinf);
            W2L5=2/(1+rinf)^2;
            W1L6=2/(1+rinf);
            l1=1;
            l2=1/2;
            l3=1/(1+rinf)^2;
            l4=1;
            l5=(3-rinf)/2/(1+rinf);
        case 'Newmark ACA'
            % Newmark Average Constant Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.25;
            W3L3=0.25;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.25;
            l4=1;
            l5=0.5;
        case 'Newmark LA'
            % Newmark Linear Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/6;
            W3L3=1/6;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/6;
            l4=1;
            l5=0.5;
        case 'Newmark BA'
            % Newmark Backward Acceleration method
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=0.5;
            W3L3=0.5;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=0.5;
            l4=1;
            l5=0.5;
        case 'Fox-Goodwin'
            % Fox-Goodwin formula
            w1=-15;
            w2=45;
            w3=-35;
            W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); % definition
            % W2=(1/3+w1/4+w2/5+w3/6)/(1+w1/2+w2/3+w3/4); % definition
            % W3=(1/4+w1/5+w2/6+w3/7)/(1+w1/2+w2/3+w3/4); % definition
            W1L1=1;
            W2L2=1/12;
            W3L3=1/12;
            W1L4=0.5;
            W2L5=0.5;
            W1L6=1;
            l1=1;
            l2=0.5;
            l3=1/12;
            l4=1;
            l5=0.5;
        otherwise
            error('No appropriate AlgID specified.');
    end
end
end

