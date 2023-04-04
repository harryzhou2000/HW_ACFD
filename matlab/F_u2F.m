function F = F_u2F(u,gamma)

rho = u(1,:);
Ux = u(2,:)./rho;
Uy = u(3,:)./rho;
E   = u(4,:)./rho;
Usqr = Ux.^2 + Uy.^2;
p = (E - 0.5 * Usqr) .* rho * (gamma - 1);

F = [u(2,:);...
    u(2,:).*Ux + p;...
    u(3,:).*Ux;...
    (u(4,:) + p).*Ux];


