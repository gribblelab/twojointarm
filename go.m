%% GET ARM ANTHROPOMETRIC PARAMETERS
arm_params = get_aparams(); % get limb anthropometric parameters

%% KINEMATICS
% Add noise to joint angles and examine distribution of hand positions

A1 = [55,90] * pi/180; % 55,90 degrees shoulder, elbow
H1 = joints_to_hand(A1, arm_params);

n_perts = 500;
A = zeros(n_perts,2);
H = zeros(n_perts,2);
r = 1.0 * pi/180; % random perts: gaussian with 1 degree sd
for i=1:n_perts
   A(i,:) = A1 + randn(1,2)*r;
   H(i,:) = joints_to_hand(A(i,:), arm_params);
end

figure
subplot(1,2,1)
plot(A(:,1)*180/pi, A(:,2)*180/pi, '.b');
hold on
plot(A1(1)*180/pi, A1(2)*180/pi, 'rs');
xlabel('SHOULDER ANGLE (deg)')
ylabel('ELBOW ANGLE (deg)')
axis equal
grid on
subplot(1,2,2)
plot(H(:,1)*1000, H(:,2)*1000, '.b');
hold on
plot(H1(1)*1000, H1(2)*1000, 'rs');
xlabel('HAND X (mm)')
ylabel('HAND Y (mm)')
axis equal
grid on

%% DYNAMICS
% Add noise to joint torques and examine distribution of trajectory endpoint

A1 = [55,90]*pi/180; % 55,90 degrees shoulder, elbow
H1 = joints_to_hand(A1, arm_params);
H2 = H1 + [0, 0.15]; % 15 cm movement straight ahead
mt = 0.500;          % movement time (sec)
sr = 100;            % sample rate (Hz)
npts = mt*sr+1;      % number of time points
dt = 1/sr;           % delta time (sec) between each sample

% get a minimum-jerk hand trajectory
[t,H,Hd,Hdd] = minjerk(H1,H2,mt,npts);

% get corresponding desired joint angles velocities and accelerations
[A,Ad,Add] = hand_to_joints(H,Hd,Hdd,arm_params);

% compute required joint torques
Q = inverse_dynamics(A,Ad,Add,arm_params);

% run a forward simulation using those joint torques Q
A0 = A(1,:);        % initial joint angles
Ad0 = Ad(1,:);      % initial joint velocities
state0 = [A0, Ad0]; % initial state parameters

% add stuff to parameters for simulation
sim_params.t = t;
sim_params.Q = Q;

% set up the state derivative function
fdfun = @(t,state) forward_dynamics(state, t, arm_params, sim_params);

% run the forward simulation
[t,state] = ode45(fdfun, t, state0);

% unpack the state into angles and velocities
A_sim = state(:,1:2);
Ad_sim = state(:,3:4);

figure
subplot(2,2,1)
plot(t,A*180/pi);
hold on
plot(t,A_sim*180/pi);
xlabel('TIME (sec)')
ylabel('JOINT ANGLES (deg)')
subplot(2,2,2)
plot(t,Ad*180/pi);
xlabel('TIME (sec)')
ylabel('JOINT VELOCITIES (deg/s)')
subplot(2,2,3)
plot(t,Add*180/pi)
xlabel('TIME (sec)')
ylabel('JOINT ACCELERATIONS (deg/s/s)')
subplot(2,2,4)
plot(t,Q)
xlabel('TIME (sec)')
ylabel('JOINT MUSCLE TORQUES (Nm)')

%% Forward simulations with random amplification or diminishment of shoulder and elbow torques

H = joints_to_hand(A,arm_params);

% re-do forward simulations and randomly (gaussian) amplify or diminish
% joint torques within a given range, and replot. In particular example
% hand endpoint distribution

figure;
subplot(1,2,1)
plot(H(1,1)*1000,H(1,2)*1000,'r.')
hold on
plot(H(end,1)*1000,H(end,2)*1000,'rs')
xlabel('X (mm)')
xlabel('Y (mm)')
subplot(1,2,2)
xlabel('X (mm)')
xlabel('Y (mm)')
hold on

nperts = 500;
for i=1:nperts
    QQ(:,1) = Q(:,1) * (0.05 * randn + 1.0); % scale by mean 0.0, sd 0.05
    QQ(:,2) = Q(:,2) * (0.05 * randn + 1.0);
    A0 = A(1,:);
    Ad0 = Ad(1,:);
    sim_params.t = t;
    sim_params.Q = QQ;
    fdfun = @(t,state) forward_dynamics(state, t, arm_params, sim_params);
    [t,state] = ode45(fdfun, t, state0);
    A_sim = state(:,1:2);
    Ad_sim = state(:,3:4);
    H_sim = joints_to_hand(A_sim, arm_params);
    subplot(1,2,1)
    plot(H_sim(:,1)*1000,H_sim(:,2)*1000,'b-')
    subplot(1,2,2)
    plot(H_sim(end,1)*1000,H_sim(end,2)*1000,'b.')    
end
plot(H(end,1)*1000,H(end,2)*1000,'rs')
subplot(1,2,1); axis equal
subplot(1,2,2); axis equal


