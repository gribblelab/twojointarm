function [M,C] = compute_dynamics_terms(A,Ad,aparams)

a1 = A(1);
a2 = A(2);
a1d = Ad(1);
a2d = Ad(2);
l1 = aparams.l(1);
l2 = aparams.l(2);
m1 = aparams.m(1);
m2 = aparams.m(2);
i1 = aparams.i(1);
i2 = aparams.i(2);
r1 = aparams.r(1);
r2 = aparams.r(2);
M = [0 0;0 0];
M(1,1) = i1 + i2 + (m1*r1*r1) + (m2*((l1*l1) + (r2*r2) + (2*l1*r2*cos(a2))));
M(1,2) = i2 + (m2*((r2*r2) + (l1*r2*cos(a2))));
M(2,1) = M(1,2);
M(2,2) = i2 + (m2*r2*r2);
C = [0; 0];
C(1) = -(m2*l1*a2d*a2d*r2*sin(a2)) - (2*m2*l1*a1d*a2d*r2*sin(a2));
C(2) = m2*l1*a1d*a1d*r2*sin(a2);

end