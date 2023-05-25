N = 200;
pnum = 7;

ARSType = 1;
ARSFix = 0;
if ARSType == 1
    if ARSFix == 0
        Name = 'Roe';
    elseif ARSFix == 1
        Name = 'Roe Entropy Fixed';
    elseif ARSFix == 2
        Name = 'Roe Entropy Fixed V1';
    else
        error('input');
    end
elseif ARSType == 2
    Name = 'HLL';
elseif ARSType == 3
    Name = 'HLLC';
else
    error('input');
end

CFL = 0.5;



if pnum <= 7
    [rhoL,uL,pL,rhoR,uR,pR] = getRiemannProblem(pnum);
end

gamma = 1.4;

[tMax, vMax, xEnd, rhoEnd, uEnd, pEnd, D] = Riemann(pnum, gamma, false);



tMax = tMax * 0.5;
if pnum == 7
    tMax = 0.24;
end
UL = f_prim2cons(rhoL,uL,pL,gamma,1);
UR = f_prim2cons(rhoR,uR,pR,gamma,1);

xs = linspace(-0.5,0.5,N+1);

xc = (xs(2:end) + xs(1:end-1))/2;
xc = reshape(xc, 1, []);
h = xs(2) - xs(1);

ithis = 1:N;
ileft = circshift(ithis,1);
iright = circshift(ithis,-1);

u0 = (xc < 0) .* UL + (xc > 0) .* UR;

dt = CFL * h/vMax;

ARSFix0 = ARSFix;

for iT = 0:2
    dtM = 1;
    O2damp = 2;
    if iT > 0
        O2damp = iT - 1 + 1;
        
        
        if iT == 1
            cName = 'TVD';
        elseif iT == 2
            cName = 'WENO5';
        end
    end
    
    u = u0;
    t = 0;
    
    figure(3); clf; hold off;
    
    for iter = 1:10000000
        dtc = dt * dtM;
        tnext = t + dtc;
        ifend = false;
        if tnext >= tMax
            dtc = tMax - t;
            ifend = true;
        end
        
        dudt0 = frhs(u, gamma, ARSType, ARSFix, h, h, h, ithis, ileft,iright, O2damp);
        u1 = u + dtc * dudt0;
        dudt1 = frhs(u1, gamma, ARSType, ARSFix, h, h, h, ithis, ileft,iright, O2damp);
        u2 = 3/4*u + 1/4*u1 + 1/4*dtc*dudt1;
        dudt2 = frhs(u2, gamma, ARSType, ARSFix, h, h, h, ithis, ileft,iright, O2damp);
        u3 = 1/3*u + 2/3*u2 + 2/3*dtc*dudt2;
        u = u3;
        
        t = t + dtc;
        if ifend
            break;
        end
        if mod(iter, 10) == 0
            plot(xc,u(1,:),'-*');
            drawnow;
        end
        
    end
    if iT ==0
        [rho,velo,p] = f_cons2prim(u,gamma,1);
    else
        [rhos{iT},velos{iT},ps{iT}] = f_cons2prim(u,gamma,1);
        names{iT} = cName;
    end
    
    
end



%%

figure(4); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1200,400]);
% t = tiledlayout(1,3,'TileSpacing','Compact' );
if pnum < 6
    % title(t, sprintf("Problem %d, %s Flux, t=%.4g, N=%d", pnum, Name, tMax, N));
    sgtitle(t, sprintf("Problem %d, %s Flux, t=%.4g, N=%d", pnum, Name, tMax, N));
end

% nexttile;
subplot(1,3,1);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, rhoEnd,'DisplayName', 'Exact')
plot(xc, rho,'x','DisplayName', 'WENO5')
xlabel('x');
ylabel('\rho');
grid on;

% nexttile;
subplot(1,3,2);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, uEnd,'DisplayName', 'Exact')
plot(xc, velo,'x','DisplayName', 'WENO5')
xlabel('x');
ylabel('u');
grid on;

% nexttile;
subplot(1,3,3);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, pEnd,'DisplayName', 'Exact')
plot(xc, p,'x','DisplayName', 'WENO5')
xlabel('x');
ylabel('p');
grid on;
l = legend;
l.Location = 'best';





print(gcf,sprintf("p%d_%s_N%d.png", pnum,Name, N),'-dpng','-r300')

%%
lineStyles = {'.-','.-','r'};


figure(5); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1200,400]);
% t = tiledlayout(1,3,'TileSpacing','Compact' );
if pnum < 6
    % title(t, sprintf("Problem %d, t=%.4g, N=%d", pnum, tMax, N));
    sgtitle(t, sprintf("Problem %d, t=%.4g, N=%d", pnum, tMax, N));
end

% nexttile;
subplot(1,3,1);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, rhoEnd,'DisplayName', 'Exact')
for iT = 1:2
    plot(xc, rhos{iT},lineStyles{iT},'DisplayName', names{iT})
end
xlabel('x');
ylabel('\rho');
grid on;

% nexttile;
subplot(1,3,2);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, uEnd,'DisplayName', 'Exact')
for iT = 1:2
    plot(xc, velos{iT},lineStyles{iT},'DisplayName', names{iT})
end
xlabel('x');
ylabel('u');
grid on;

% nexttile;
subplot(1,3,3);
hold on; set(gca, 'FontName', 'Times New Roman');
plot(xEnd, pEnd,'DisplayName', 'Exact')
for iT = 1:2
    plot(xc, ps{iT},lineStyles{iT},'DisplayName', names{iT})
end
xlabel('x');
ylabel('p');
grid on;
l = legend;
l.Location = 'best';





print(gcf,sprintf("p%d_SUM_N%d.png", pnum, N),'-dpng','-r300')
%%
function dudt = frhs(u, gamma, ARSType, ARSFix, hs, hl, hr, ithis, ileft,iright, O2damp)


if O2damp <= 1
    [rhoRoe,veloRoe,pRoe] = f_cons2prim(u,gamma,1);
    assert(all(pRoe > 0));
    assert(all(rhoRoe > 0));
    aRoe = sqrt(gamma * pRoe ./ rhoRoe);
    asqrRoe = aRoe.^2;
    vsqrRoe =  sum(veloRoe.^2,1);
    HRoe = asqrRoe/(gamma-1) + 0.5 *vsqrRoe;
    
    UXL = (u - u(:,ileft))./hl;
    UXR = (u(:,iright) - u)./hr;
    WXL = EulerCons2Char(veloRoe,vsqrRoe,aRoe,asqrRoe,HRoe,UXL,1,gamma);
    WXR = EulerCons2Char(veloRoe,vsqrRoe,aRoe,asqrRoe,HRoe,UXR,1,gamma);
    
    WX = F_TVD_Slope(WXL,WXR) * O2damp;
    UX = EulerChar2Cons(WX, 1, 1, 1,veloRoe,vsqrRoe,aRoe,asqrRoe,HRoe,1,gamma);
    
    
    
    UX(:,1) = 0;
    UX(:,end) = 0;
    uL = u - UX.*hs * 0.5;
    uR = u + UX.*hs * 0.5;
else
    if(ARSType ~= 1)
        shape0 = size(u);
        [uL,uR] = F_interpi_weno5(permute(reshape(u,shape0(1),1,shape0(2)),[3,2,1]), 1e-7,2,1);
        uL = reshape(permute(uL,[3,2,1]),shape0);
        uR = reshape(permute(uR,[3,2,1]),shape0);
    end
end




if ARSType == 1 && O2damp > 1
    FL = EulerFluxFDi(u(:,ileft), u, gamma, 1, ARSFix);
else
    FL = EulerFluxARS(uR(:,ileft),uL,gamma,1,ARSType, ARSFix);
end

dudt = (FL - FL(:,iright)) ./ hs;
dudt(:,1:2) = 0;
dudt(:,end-1:end) = 0;

end


