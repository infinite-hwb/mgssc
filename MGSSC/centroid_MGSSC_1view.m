function [C1, C2, Lambda1, Lambda2, A] = centroid_MGSSC_1view(n, K, C1, C2, C_centroid, Lambda1, Lambda2, lambda, mu)

% updating A

inv_A = inv(2*K+(mu(1)+mu(2))*eye(n));
A = inv_A * (2*K+mu(1)*(C1-Lambda1/mu(1))+mu(2)*(C2-Lambda2/mu(2)));


% updating C2 and C3
C1_new = soft_thresh(A+Lambda1/mu(1),lambda(1)/mu(1));
C1_new = C1_new - diag(diag(C1_new));
C2_new = 1/(2*lambda(2)+mu(2))*(2*lambda(2)*C_centroid+mu(2)*A+Lambda2);

%update variables
C1 = C1_new;
C2 = C2_new;

% updating Lagrange multipliers

Lambda1 = Lambda1 + mu(1) * (A - C1 + diag(diag(C1)));
Lambda2 = Lambda2 + mu(2) * (A - C2);
