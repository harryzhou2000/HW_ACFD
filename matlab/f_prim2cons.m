function U = f_prim2cons(rho,velo,p,gamma,d)

U = nan(d+2,size(rho,2));
U(1,:) = rho;
U(2:1+d,:) = velo .* rho;

Ek = 0.5 * rho .* sum(velo.^2,1);
E = p/(gamma - 1) + Ek;
U(end,:) = E;

end