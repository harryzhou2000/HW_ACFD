function [tMax, vMax, xEnd, rhoEnd, uEnd, pEnd, D] = Riemann(pnum, gamma, ifdraw)


N = 1024;



[rhoL,uL,pL,rhoR,uR,pR] = getRiemannProblem(pnum);


aL = sqrt(gamma * pL / rhoL);
aR = sqrt(gamma * pR / rhoR);





UL = f_prim2cons(rhoL,uL,pL,gamma,1);
UR = f_prim2cons(rhoR,uR,pR,gamma,1);
[F,D] = EulerRiemannFluxExact_1D(UL,UR,gamma);

vMax = max(abs(uL), abs(uR)) + max(aL, aR);
vMax = max(abs([D.sLL, D.sLR, D.sRR, D.sRL]));

tMax = 0.8/vMax;
if(pnum == 7)
   tMax = 0.24; 
end


xs = linspace(-1,1,N+1);
if pnum == 7
xs = linspace(-0.5,0.5,N+1);
end
ts = linspace(0,tMax,N+1);

xEnd = xs;



[xm,tm] = meshgrid(xs,ts);
rhom = nan(size(xm));
um = rhom;
pm = rhom;


for i = 1:numel(xm)
    xxx = xm(i);
    ttt = tm(i);
    xLL = D.sLL * ttt;
    xLR = D.sLR * ttt;
    xRL = D.sRL * ttt;
    xRR = D.sRR * ttt;
    xC = D.Vm * ttt;
    
    if(xxx < xLL)
        rhom(i) = D.rhoL;
        um(i) = D.VL;
        pm(i) = D.pL;
    elseif(xxx < xLR)
        um(i) = 2/(gamma+1) * xxx/(ttt + 1e-15) + (gamma-1)/(gamma+1) * D.VL + 2/(gamma+1) * D.aL;
        a = -(gamma-1)/(gamma+1) * xxx/(ttt + 1e-15) + (gamma-1)/(gamma+1) * D.VL + 2/(gamma+1) * D.aL;
        rhom(i) = (a/D.aL).^(2/(gamma-1)) * D.rhoL;
        pm(i) = (a/D.aL).^(2*gamma/(gamma-1)) * D.pL;
        if( imag(rhom(i)))
            error('imag');
        end
        
    elseif(xxx < xC)
        um(i) = D.Vm;
        pm(i) = D.pm;
        rhom(i) = D.rhomL;
    elseif(xxx == xC)
        um(i) = D.Vm;
        pm(i) = D.pm;
        rhom(i) = (D.rhomL + D.rhomR)/2;
    elseif(xxx <= xRL)
        um(i) = D.Vm;
        pm(i) = D.pm;
        rhom(i) = D.rhomR;
    elseif(xxx <= xRR)
        um(i) = 2/(gamma+1) * xxx/(ttt + 1e-15) + (gamma-1)/(gamma+1) * D.VR - 2/(gamma+1) * D.aR;
        a = +(gamma-1)/(gamma+1) * xxx/(ttt + 1e-15) - (gamma-1)/(gamma+1) * D.VR + 2/(gamma+1) * D.aR;
        rhom(i) = (a/D.aR).^(2/(gamma-1)) * D.rhoR;
        pm(i) = (a/D.aR).^(2*gamma/(gamma-1)) * D.pR;
        if( imag(rhom(i)))
            error('imag');
        end
    else
        rhom(i) = D.rhoR;
        um(i) = D.VR;
        pm(i) = D.pR;
    end
    
    if( mod(i, 1000) == 0)
        fprintf("i = %d\n", i);
        
    end
    
end

D

rhoEnd = rhom(end,:);
uEnd = um(end,:);
pEnd = pm(end,:);



if ifdraw
    %%
    close all;
    figure(1); clf; set(gca, 'FontName', 'Times New Roman');
    set(gcf,'Position',[100,100,1200,400]);
    t = tiledlayout(1,3,'TileSpacing','Compact' );
    title(t, sprintf("Problem %d", pnum));
    
    
    nexttile; hold on;set(gca, 'FontName', 'Times New Roman');
    p = pcolor(xm,tm,rhom);
    p.LineStyle = 'none';
    colorbar;
    xlabel('x');
    ylabel('t');
    title('\rho');
    if abs(D.rhomL - D.rhomR) > 1e-10
        plot(ts * D.Vm, ts, '--k', 'LineWidth', 2);
    end
    if D.pm > D.pL + 1e-10
        plot(ts * D.sLL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pL - 1e-10
        plot(ts * D.sLL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sLR, ts, ':k', 'LineWidth', 2);
    end
    if D.pm > D.pR + 1e-10
        plot(ts * D.sRL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pR - 1e-10
        plot(ts * D.sRL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sRR, ts, ':k', 'LineWidth', 2);
    end
    axis tight;
    
    nexttile; hold on;set(gca, 'FontName', 'Times New Roman');
    p = pcolor(xm,tm,um);
    p.LineStyle = 'none';
    colorbar;
    xlabel('x');
    ylabel('t'); set(gca,'YTick',[], 'YLabel',[]);
    title('u');
    if abs(D.rhomL - D.rhomR) > 1e-10
        plot(ts * D.Vm, ts, '--k', 'LineWidth', 2);
    end
    if D.pm > D.pL + 1e-10
        plot(ts * D.sLL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pL - 1e-10
        plot(ts * D.sLL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sLR, ts, ':k', 'LineWidth', 2);
    end
    if D.pm > D.pR + 1e-10
        plot(ts * D.sRL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pR - 1e-10
        plot(ts * D.sRL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sRR, ts, ':k', 'LineWidth', 2);
    end
    axis tight;
    
    nexttile;  hold on;set(gca, 'FontName', 'Times New Roman');
    p = pcolor(xm,tm,pm);
    p.LineStyle = 'none';
    colorbar;
    xlabel('x');
    ylabel('t'); set(gca,'YTick',[], 'YLabel',[]);
    title('p');
    if abs(D.rhomL - D.rhomR) > 1e-10
        plot(ts * D.Vm, ts, '--k', 'LineWidth', 2);
    end
    if D.pm > D.pL + 1e-10
        plot(ts * D.sLL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pL - 1e-10
        plot(ts * D.sLL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sLR, ts, ':k', 'LineWidth', 2);
    end
    if D.pm > D.pR + 1e-10
        plot(ts * D.sRL, ts, '-k', 'LineWidth', 2);
    elseif D.pm < D.pR - 1e-10
        plot(ts * D.sRL, ts, ':k', 'LineWidth', 2);
        plot(ts * D.sRR, ts, ':k', 'LineWidth', 2);
    end
    axis tight;
    
    
    print(gcf,sprintf("p%d.png", pnum),'-dpng','-r300')
    
    %%
    figure(2); clf; set(gca, 'FontName', 'Times New Roman');
    set(gcf,'Position',[100,100,1200,400]);
    t = tiledlayout(1,3,'TileSpacing','Compact' );
    title(t, sprintf("Problem %d, t=%.4g", pnum, tMax));
    
    
    nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
    plot(xs, rhom(end,:))
    xlabel('x');
    ylabel('\rho');
    grid on;
    
    nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
    plot(xs, um(end,:))
    xlabel('x');
    ylabel('u');
    grid on;
    
    nexttile; hold on; set(gca, 'FontName', 'Times New Roman');
    plot(xs, pm(end,:))
    xlabel('x');
    ylabel('p');
    grid on;
    
    
    
    
    
    print(gcf,sprintf("p%d_l.png", pnum),'-dpng','-r300')
    
end


