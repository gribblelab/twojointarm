function aparams = get_aparams()

% from Winter book:
% M=75; % kg
% H=1.8; % m
% l1=.188*H; l2=.253*H;
% m1=.028*M; m2=.022*M;
% rog1=.322*l1; rog2=.468*l2;
% i1=m1*rog1^2; i2=m2*rog2^2;

aparams.l = [0.3384, 0.4554]; % link lengths s,e (metres)
aparams.r = [0.1692, 0.2277]; % radius of gyration (metres)
aparams.m = [2.10, 1.65];     % mass (kg)
aparams.i = [0.025, 0.075];   % moment of inertia (kg*m*m)

end
