
see = 10;
Nx = 100;
Ny = 1;

eps = 1e-6;
p = 2e200;
nv = 2;

xs = linspace(-1/2+0.5/Nx,1/2-0.5/Nx,Nx);
ys = linspace(0.5,0.5,Ny);

[xm,ym] = meshgrid(xs,ys);
xm = permute(xm,[2,1]);
ym = permute(ym,[2,1]);

ax = 1;
ay = 0;

us = zeros(Nx,Ny,nv);

us(:,:,1) = (xm - 0).^2 < 0.25^2;
% us(:,:,1) = cos(xm * 2 *pi) .* cos(ym * 2*pi);

G.xm = xm;
G.ym = ym;
G.hx = 1/Nx;
G.hy = 1/Ny;
G.ithis = 1:Nx;
G.iri = circshift(G.ithis, -1);
G.ile = circshift(G.ithis,  1);
G.jthis = 1:Ny;
G.jup = circshift(G.jthis, -1);
G.jlo = circshift(G.jthis,  1);
G.eps = eps;
G.p = p;
M.ax = ax;
M.ay = ay;

dt = G.hx / ax * 0.5;
Tmax = 0.1;

Niter = round(Tmax/dt);

u0 = us;
u = us;
u = gpuArray(u);

figure(1);clf;cla;

for iter = 1:Niter
    %     R0 = frhs(u ,G,M);
    %     u1 = u + dt/2 * R0;
    %     R1 = frhs(u1 ,G,M);
    %     u2 = u + dt/2 * R1;
    %     R2 = frhs(u2, G,M);
    %     u3 = u + dt   * R2;
    %     R3 = frhs(u3, G,M);
    %     unew = u + dt/6 * (R0 + 2*R1 + 2*R2 + R3);
    R0 = frhs(u ,G,M);
    u1 = u + dt * R0;
    R1 = frhs(u1 ,G,M);
    u2 = 3/4 * u + 1/4*u1 + dt /4 * R1;
    R2 = frhs(u2, G,M);
    unew = 1/3*u + 2/3*u2 + 2*dt/3 * R2;
    
    u = unew;
    
    if (mod(iter,see) == 0 || iter == Niter)
        plot(G.xm,u(:,:,1))
        view(0,90)
        drawnow;
    end
end
err0 = norm(u(:) - us(:), 1) * G.hx * G.hy;
fprintf("err0 = %.7e\n", err0);



hold on;
plot(G.xm, double(abs(mod(G.xm - ax * Tmax + 0.5, 1) -0.5) < 0.25));
xlabel('x');
ylabel('u');
if( p <1e100)
    title(sprintf("WENO5 \\epsilon = %.e",eps));
else
    title(sprintf("ENO3 t = %g", Tmax));
end
legend('solution','exact');
grid on; grid minor;
ylim([-0.2,1.2]);


function dudt = frhs(u, G, M)

[uLe, uRi] = F_interpi_weno5(u, G.eps, G.p, 0);
[uLo, uUp] = F_interpi_weno5(permute(u,[2,1,3]), G.eps, G.p, 0);
uLo = permute(uLo, [2,1,3]);
uUp = permute(uUp, [2,1,3]);

fLe_uRi = uLe;
fLe_uLe = uRi(G.ile,:,:);
fRi_uLe = uRi;
fRi_uRi = uLe(G.iri,:,:);

fLo_uUp = uLo;
fLo_uLo = uUp(:,G.jlo,:);
fUp_uLo = uUp;
fUp_uUp = uLo(:,G.jup,:);

fLe_f = RS(fLe_uLe, fLe_uRi, M.ax);
fRi_f = RS(fRi_uLe, fRi_uRi, M.ax);
fLo_f = RS(fLo_uLo, fLo_uUp, M.ay);
fUp_f = RS(fUp_uLo, fUp_uUp, M.ay);

dudt = (fLe_f - fRi_f) / G.hx + (fLo_f - fUp_f) / G.hy;

end


function f = RS(uL,uR,a)
f = 0.5 * (a * uL + a * uR - abs(a) * (uR - uL));

end




