function [rho,velo,p] = f_cons2prim(U,gamma,d)
rho = U(1,:);
velo = U(2:1+d,:)./rho;
Ek = 0.5 * rho .* sum(velo.^2,1);
E = U(2+d,:);
p = (E - Ek) * (gamma - 1);

end