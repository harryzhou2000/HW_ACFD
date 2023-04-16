N = 128;


pnum = 5;
[rhoL,uL,pL,rhoR,uR,pR] = getRiemannProblem(pnum);


gamma = 1.4;

[tMax, vMax, xEnd, rhoEnd, uEnd, pEnd, D] = Riemann(pnum, gamma);
UL = f_prim2cons(rhoL,uL,pL,gamma,1);
UR = f_prim2cons(rhoR,uR,pR,gamma,1);

xs = linspace(-1,1,N+1);
xc = (xs(2:end) + xs(1:end-1))/2;
xc = reshape(xc, 1, []);
h = xs(2) - xs(1);

ithis = 1:N;
ileft = circshift(ithis,1);
iright = circshift(ithis,-1);

u0 = (xc < 0) .* UL + (xc > 0) .* UR;

dt = 0.5 * h/vMax;

u = u0;
t = 0;

figure(3); clf; hold off;

for iter = 1:10000000
    dtc = dt;
    tnext = t + dtc;
    ifend = false;
    if tnext >= tMax
        dtc = tMax - t;
        ifend = true;
    end
    
    F = EulerRiemannFluxExact_1D(u(:,ileft),u,gamma,1e-4);
    dudt = (F - F(:,iright)) / h;
    dudt(:,1) = 0;
    dudt(:,end) = 0;
    u = u + dtc * dudt;
    
    t = t + dt;
    if ifend
        break;
    end
    if mod(iter, 10 == 0)
       plot(xc,u(1,:),'-*');
       drawnow;
    end
    
end

[rho,velo,p] = f_cons2prim(u,gamma,1);

%%

figure(4); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1200,400]);
t = tiledlayout(1,3,'TileSpacing','Compact' );
title(t, sprintf("Problem %d, t=%.4g", pnum, tMax));


nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, rhoEnd,'DisplayName', 'Exact')
plot(xc, rho,'x','DisplayName', 'Godunov')
xlabel('x');
ylabel('\rho');
grid on;

nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, uEnd,'DisplayName', 'Exact')
plot(xc, velo,'x','DisplayName', 'Godunov')
xlabel('x');
ylabel('u');
grid on;

nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, pEnd,'DisplayName', 'Exact')
plot(xc, p,'x','DisplayName', 'Godunov')
xlabel('x');
ylabel('p');
grid on;
l = legend;
l.Location = 'best';





print(gcf,sprintf("p%d_gn.png", pnum),'-dpng','-r300')




