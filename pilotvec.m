function [final_value]=pilotvec(legth_tau)

A_mat=randn(2*legth_tau);
while (rank(A_mat)~=2*legth_tau)
A_mat=randn(2*legth_tau);
end
Orth_A=orth(A_mat);

final_value=(1/sqrt(2))*(Orth_A(:,1:legth_tau) +1i*Orth_A(:,1+legth_tau:2*legth_tau));

end