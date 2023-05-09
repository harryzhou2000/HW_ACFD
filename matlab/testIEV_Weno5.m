
see = 50;
sup = 4;
Nx = 10 * sup;
Ny = 10 * sup;
C = 0.1;
Tmax = 2;
G.p = 1e200;
G.eps = 1e-6; 
G.mapping = 1;
% eps 10* 2 4 8 16 32 
% 1e-6    3.1230530e-01 5.2543401e-02 4.4590657e-03 2.5012196e-04 1.1233154e-05 3.7072570e-07
% 1e-3    2.9349690e-01 3.2229027e-02 2.4662220e-03 8.1062336e-05 1.9013022e-06 4.4152115e-08
% 1e200   9.5435399e-02 1.2498203e-02 6.0787613e-04 2.0860055e-05 6.6494326e-07 2.0892949e-08
% 1e-6m   2.2693522e-01 2.9438651e-02 9.2472577e-04 2.2989118e-05 6.6682832e-07 2.0893558e-08
% 1e-6Roe 

nv = 4;

Lx = 10;
Ly = 10;

xs = linspace(0,Lx-Lx/Nx,Nx);
ys = linspace(0,Ly-Ly/Ny,Ny);

[xm,ym] = meshgrid(xs,ys);
xm = permute(xm,[2,1]);
ym = permute(ym,[2,1]);

ax = 1;
ay = 1;



gamma = 1.4;

chi = 5;
rm = sqrt((xm - 5).^2 + (ym - 5).^2);
dux = chi/2/pi * exp((1-rm.^2)/2) .* -(ym - 5);
duy = chi/2/pi * exp((1-rm.^2)/2) .*  (xm - 5);
dT = -(gamma - 1)/(8*gamma*pi^2)*chi^2 * exp(1-rm.^2);
T = dT + 1;
ux = dux + 1;
uy = duy + 1;
S = 1;
rho = (T./S).^(1/(gamma-1));
p = T.*rho;

p(1,:) = 1;
p(:,1) = 1;
rho(1,:) = 1;
rho(:,1) = 1;
ux(:,1) = 1;
ux(1,:) = 1;
uy(:,1) = 1;
uy(1,:) = 1;
E = 0.5 * (ux.^2 + uy.^2).*rho + p / (gamma -1);

rmend = sqrt((mod(xm - 7 + 5,10)-5).^2 + (mod(ym - 7 + 5,10)-5).^2);
dTend = -(gamma - 1)/(8*gamma*pi^2)*chi^2 * exp(1-rmend.^2);
rhoEnd = (dTend + 1).^(1/(gamma-1));




us = zeros(Nx,Ny,nv);

us(:,:,1) = rho;
us(:,:,2) = ux.*rho;
us(:,:,3) = uy.*rho;
us(:,:,4) = E;
% us(:,:,1) = cos(xm * 2 *pi) .* cos(ym * 2*pi);

G.xm = xm;
G.ym = ym;
G.hx = Lx/Nx;
G.hy = Ly/Ny;
G.ithis = 1:Nx;
G.iri = circshift(G.ithis, -1);
G.ile = circshift(G.ithis,  1);
G.jthis = 1:Ny;
G.jup = circshift(G.jthis, -1);
G.jlo = circshift(G.jthis,  1);
M.ax = ax;
M.ay = ay;
M.gamma = gamma;

dt = G.hx / ax * C;


Niter = round(Tmax/dt);

u = us;

%%
u = gpuArray(u);

for iter = 1:Niter
    R0 = frhs(u ,G,M);
    u1 = u + dt/2 * R0;
    R1 = frhs(u1 ,G,M);
    u2 = u + dt/2 * R1;
    R2 = frhs(u2, G,M);
    u3 = u + dt   * R2;
    R3 = frhs(u3, G,M);
    unew = u + dt/6 * (R0 + 2*R1 + 2*R2 + R3);

%     R0 = frhs(u ,G,M);
    
    
    
    u = unew;
    fprintf('iter %d\n', iter);
    if (mod(iter,see) == 0 || iter == Niter)
        surf(G.xm,G.ym,u(:,:,1),'EdgeColor','none');
        view(0,90)
        drawnow;
    end
end
% err0 = norm(u(:) - us(:), 1) * G.hx * G.hy;

p = (gamma-1) * (u(:,:,4) - 0.5 * (u(:,:,2).^2 + u(:,:,3).^2)./u(:,:,1));
T = p./u(:,:,1);
dT = T-1;

rhoc = u(:,:,1);

% err0 = norm(rho(:) - rhoc(:), 1);
rhoSee = rhoc(xm>5 & xm<9 & ym > 5 & ym <9);
rhoEndSee = rhoEnd(xm>5 & xm<9 & ym > 5 & ym <9);
err0 = sum(abs(rhoSee(:) - rhoEndSee(:))) * G.hx * G.hy;
fprintf("err0 = %.7e\n", err0);

%%
surf(G.xm,G.ym,rhoc - rhoEnd,'EdgeColor','none');
        view(0,90)
        drawnow;
colorbar;
title("\rho error, N = " + string(Nx) + sprintf(", \\epsilon = %.1e", G.eps));
xlabel('x'); ylabel('y');



function dudt = frhs(u, G, M)
p = G.p;
eps = G.eps;
    
    [uLe, uRi] = F_interpi_weno5(u, eps, p, G.mapping);
    [uLo, uUp] = F_interpi_weno5(permute(u,[2,1,3]), eps, p, G.mapping);
    uLo = permute(uLo, [2,1,3]);
    uUp = permute(uUp, [2,1,3]);
    
%     uLe = u;
%     uRi = u;
%     uUp = u;
%     uLo = u;
    
    
    fLe_uRi = uLe;
    fLe_uLe = uRi(G.ile,:,:);
    fRi_uLe = uRi;
    fRi_uRi = uLe(G.iri,:,:);

    fLo_uUp = uLo;
    fLo_uLo = uUp(:,G.jlo,:);
    fUp_uLo = uUp;
    fUp_uUp = uLo(:,G.jup,:);
    
%     fLe_f = RSx(fLe_uLe, fLe_uRi, M.gamma);
%     fRi_f = RSx(fRi_uLe, fRi_uRi, M.gamma);
%     fLo_f = RSy(fLo_uLo, fLo_uUp, M.gamma);
%     fUp_f = RSy(fUp_uLo, fUp_uUp, M.gamma);
%     dudt = (fLe_f - fRi_f) / G.hx + (fLo_f - fUp_f) / G.hy;
    
    
    [fx, lam1, lam2, lam3] = invfluxX(u, M.gamma);
    lamaxx = max(abs(lam1),abs(lam3));
    [fy, lam1, lam2, lam3] = invfluxY(u, M.gamma);
    lamaxy = max(abs(lam1),abs(lam3));
    lammax = max(lamaxx, lamaxy);
    lammaxLe = max(lammax, lammax(G.ile,:));
    lammaxRi = max(lammax, lammax(G.iri,:));
    lammaxLo = max(lammax, lammax(:,G.jlo));
    lammaxUp = max(lammax, lammax(:,G.jup));
    
    [fxLe,fxRi] = F_interpi_weno5(fx, eps, p, G.mapping);
    [fyLo,fyUp] = F_interpi_weno5(permute(fy,[2,1,3]), eps, p, G.mapping);
    fyLo = permute(fyLo, [2,1,3]);
    fyUp = permute(fyUp, [2,1,3]);
    
    fLe_fxRi = fxLe;
    fLe_fxLe = fxRi(G.ile,:,:);
    fRi_fxLe = fxRi;
    fRi_fxRi = fxLe(G.iri,:,:);
    
    fLo_fyUp = fyLo;
    fLo_fyLo = fyUp(:,G.jlo,:);
    fUp_fyLo = fyUp;
    fUp_fyUp = fyLo(:,G.jup,:);
    
    fLo = 0.5 * (fLo_fyLo + fLo_fyUp - 1 * lammaxLo.*(fLo_uUp - fLo_uLo));
    fUp = 0.5 * (fUp_fyLo + fUp_fyUp - 1 * lammaxUp.*(fUp_uUp - fUp_uLo));
    fLe = 0.5 * (fLe_fxLe + fLe_fxRi - 1 * lammaxLe.*(fLe_uRi - fLe_uLe));
    fRi = 0.5 * (fRi_fxLe + fRi_fxRi - 1 * lammaxRi.*(fRi_uRi - fRi_uLe));
    
    dudt = (fLe - fRi) / G.hx + (fLo - fUp) / G.hy;
    
    
 


end

function [f, lam1, lam2, lam3] = invfluxX(u, gamma)
    rho = u(:,:,1);
    ux = u(:,:,2)./rho;
    uy = u(:,:,3)./rho;
    E = u(:,:,4);
    p = (gamma-1) * (E - 0.5 * rho .* (ux.^2 + uy.^2));
    
    f = u.*ux;
    f(:,:,2) = f(:,:,2) + p;
    f(:,:,4) = f(:,:,4) + p.*ux;
    a = sqrt(gamma * p ./rho);
    lam1 = ux - a;
    lam2 = ux;
    lam3 = ux + a;
end

function [f, lam1, lam2, lam3] = invfluxY(u, gamma)
    rho = u(:,:,1);
    ux = u(:,:,2)./rho;
    uy = u(:,:,3)./rho;
    E = u(:,:,4);
    p = (gamma-1) * (E - 0.5 * rho .* (ux.^2 + uy.^2));
    
    f = u.*uy;
    f(:,:,3) = f(:,:,3) + p;
    f(:,:,4) = f(:,:,4) + p.*uy;
    a = sqrt(gamma * p ./rho);
    lam1 = uy - a;
    lam2 = uy;
    lam3 = uy + a;
end



% function f = RSx(uL,uR,gamma)
%     [fL, l1,l2,l3] = invfluxX(uL,gamma);
%     lambdaL = max(abs(l1),abs(l3));
%     [fR, l1,l2,l3] = invfluxX(uR,gamma);
%     lambdaR = max(abs(l1),abs(l3));
%     lammax = max(lambdaL, lambdaR);
%     f = 0.5 * (fL + fR - lammax.*(uR-uL));
% end
% 
% function f = RSy(uL,uR,gamma)
%     [fL, l1,l2,l3] = invfluxY(uL,gamma);
%     lambdaL = max(abs(l1),abs(l3));
%     [fR, l1,l2,l3] = invfluxY(uR,gamma);
%     lambdaR = max(abs(l1),abs(l3));
%     lammax = max(lambdaL, lambdaR);
%     f = 0.5 * (fL + fR - lammax.*(uR-uL));
% end

function f = RSx(uL,uR,gamma)
    uL = permute(uL, [3,2,1]);
    uR = permute(uR, [3,2,1]);
    sz = size(uL);
    uL = reshape(uL,4,[]);
    uR = reshape(uR,4,[]);
    f = RoeSolver2(uL,uR,gamma);
    f = reshape(f, sz);
    f = permute(f, [3,2,1]);
end

function f = RSy(uL,uR,gamma)
    uL = permute(uL, [3,2,1]);
    uR = permute(uR, [3,2,1]);
    sz = size(uL);
    uL = reshape(uL,4,[]);
    uR = reshape(uR,4,[]);
    uL = uL([1,3,2,4],:);
    uR = uR([1,3,2,4],:);
    f = RoeSolver2(uL,uR,gamma);
    f = f([1,3,2,4],:);
    f = reshape(f, sz);
    f = permute(f, [3,2,1]);
end




