function [rho, Ux,Uy,Usqr, E,p,a,H] = F_uExpand_V1(u, gamma)

rho = u(1,:);
Ux  = u(2,:)./rho;
Uy  = u(3,:)./rho;
E   = u(4,:)./rho;
Usqr = Ux.^2 + Uy.^2;
p = (E - 0.5 * Usqr) .* rho * (gamma - 1);
a = sqrt(gamma * p ./ rho );
H = E + p./rho;