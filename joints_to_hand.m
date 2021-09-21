function [H,E] = joints_to_hand(A, arm_params)

l1 = arm_params.l(1);
l2 = arm_params.l(2);
E = [l1*cos(A(:,1))             , l1*sin(A(:,1))       ];
H = E + [l2 * cos(A(:,1)+A(:,2)), l2*sin(A(:,1)+A(:,2))];

end