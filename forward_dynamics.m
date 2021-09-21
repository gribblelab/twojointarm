function stated = forward_dynamics(state,t,aparams,sparams)

a1  = state(1);
a2  = state(2);
a1d = state(3);
a2d = state(4);
Q = interp1(sparams.t, sparams.Q, t)';
[M,C] = compute_dynamics_terms(state(1:2), state(3:4),aparams);
A = [a1;a2];
Ad = [a1d;a2d];
J = arm_jacobian(A,aparams);
Hd = J*Ad;
ACC = inv(M) * (Q - C);
stated = [a1d a2d ACC(1) ACC(2)]';

end
