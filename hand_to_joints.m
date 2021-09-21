function [A,Ad,Add] = hand_to_joints(H,Hd,Hdd,arm_params)

l1 = arm_params.l(1);
l2 = arm_params.l(2);
n = size(H,1);
A = zeros(n,2);
Ad = zeros(n,2);
Add = zeros(n,2);

A(:,2) = acos((H(:,1).^2 + H(:,2).^2 - l1^2 - l2^2) / (2*l1*l2));
A(:,1) = atan2(H(:,2),H(:,1)) - atan2(l2*sin(A(:,2)),l1+(l2*cos(A(:,2))));
for i=1:n
    J = arm_jacobian(A(i,:),arm_params);
    Ad(i,:) = (inv(J)*Hd(i,:)')';
    Jd = arm_jacobian_d(A(i,:),Ad(i,:),arm_params);
    b = Hdd(i,:) - (Jd*Ad(i,:)')';
    Add(i,:) = (inv(J)*b')';
end

end