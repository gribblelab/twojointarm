function J = arm_jacobian(A,aparams)

l1 = aparams.l(1);
l2 = aparams.l(2);
J = zeros(2,2);
J(1,1) = -l1*sin(A(1)) - l2*sin(A(1)+A(2));
J(1,2) = -l2*sin(A(1)+A(2));
J(2,1) = l1*cos(A(1)) + l2*cos(A(1)+A(2));
J(2,2) = l2*cos(A(1)+A(2));

end