function Q = inverse_dynamics(A,Ad,Add,aparams)

n = size(A,1);
Q = zeros(n,2);
for i=1:n
   [M,C] = compute_dynamics_terms(A(i,:),Ad(i,:),aparams);
   ACC = Add(i,:)';
   J = arm_jacobian(A(i,:),aparams);
   Hd = J*Ad(i,:)';
   Q(i,:) = (M*ACC + C)';
end

end