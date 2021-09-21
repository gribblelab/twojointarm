function Jd = arm_jacobian_d(A,Ad,aparams)

l1 = aparams.l(1);
l2 = aparams.l(2);
Jd = zeros(2,2);
Jd(1,1) = -l1*cos(A(1))*Ad(1) - l2*(Ad(1) + Ad(2))*cos(A(1) + A(2));
Jd(1,2) = -l2*(Ad(1) + Ad(2))*cos(A(1) + A(2));
Jd(2,1) = -l1*sin(A(1))*Ad(1) - l2*(Ad(1) + Ad(2))*sin(A(1) + A(2));
Jd(2,2) = -l2*(Ad(1) + Ad(2))*sin(A(1) + A(2));

end