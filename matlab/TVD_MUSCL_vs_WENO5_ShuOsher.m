pnum = 0;

NRef = 10000;
%  [rhoR,veloR,pR,name,xc] = getShuOsher(NRef,2);
%  save('ShuOsherRefA.mat','rhoR','veloR','pR');
%%

RefSH = load('ShuOsherRefA.mat');
Refx = linspace(0 + 0.5/numel(RefSH.pR), 10 - 0.5/numel(RefSH.pR), numel(RefSH.pR));
%%
Ns = [200, 400, 800, 5000];

rhos = {};
velos = {};
ps = {};
name = {};
xcs = {};

for N = Ns
   [rho,velo,p,name,xc] = getShuOsher(N,1);
   rhos{end+1} = rho;
   velos{end+1} = velos;
   ps{end+1} = p;
   names{end+1} = name;
   xcs{end+1} = xc;
end

for N = Ns
   [rho,velo,p,name,xc] = getShuOsher(N,2);
   rhos{end+1} = rho;
   velos{end+1} = velos;
   ps{end+1} = p;
   names{end+1} = name;
   xcs{end+1} = xc;
end


%%


lineStyles = {'.-','.-','r'};


figure(5); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1000,500]);

% t = tiledlayout(2,2,'TileSpacing','Compact' );


% title(t, sprintf("Problem %d, t=%.4g, N=%d", pnum, tMax, N));

for i = 1:numel(Ns)

% nexttile;
subplot(2,2,i);
hold on; set(gca, 'FontName', 'Times New Roman');

title(sprintf('N = %d', Ns(i)));
plot(Refx, RefSH.rhoR,'DisplayName', 'Reference');
plot(xcs{i},  rhos{i}, '.-', 'DisplayName', 'TVD',...
    'MarkerIndices',1:10:numel(xcs{i}));
plot(xcs{i + numel(Ns)},  rhos{i + numel(Ns)}, 'x-', 'DisplayName', 'WENO5',...
    'MarkerIndices',1:10:numel(xcs{i}));
L = legend();
L.Location = 'best';

xlabel('x');
ylabel('\rho');
grid on;

end
print(gcf,sprintf("SH_SUM0.png"),'-dpng','-r300')

%%
figure(6); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1000,500]);

% t = tiledlayout(2,2,'TileSpacing','Compact' );

% title(t, sprintf("Problem %d, t=%.4g, N=%d", pnum, tMax, N));

for i = 1:numel(Ns)

% nexttile;
subplot(2,2,i);
hold on; set(gca, 'FontName', 'Times New Roman');
title(sprintf('N = %d', Ns(i)));
plot(Refx, RefSH.rhoR,'DisplayName', 'Reference');
plot(xcs{i},  rhos{i}, '.-', 'DisplayName', 'TVD',...
    'MarkerIndices',1:10:numel(xcs{i}));
plot(xcs{i + numel(Ns)},  rhos{i + numel(Ns)}, 'x-', 'DisplayName', 'WENO5',...
    'MarkerIndices',1:10:numel(xcs{i}));
L = legend();
L.Location = 'best';
xlim([5,7.5]);
ylim([3,5]);

xlabel('x');
ylabel('\rho');
grid on;

end
print(gcf,sprintf("SH_SUM1.png"),'-dpng','-r300')

%%
figure(7); clf; set(gca, 'FontName', 'Times New Roman');
set(gcf,'Position',[100,100,1000,500]);

% t = tiledlayout(2,2,'TileSpacing','Compact' );

% title(t, sprintf("Problem %d, t=%.4g, N=%d", pnum, tMax, N));

for i = 1:numel(Ns)

% nexttile;
subplot(2,2,i);
hold on; set(gca, 'FontName', 'Times New Roman');
title(sprintf('N = %d', Ns(i)));
plot(Refx, RefSH.rhoR,'DisplayName', 'Reference');
plot(xcs{i},  rhos{i}, '.-', 'DisplayName', 'TVD');
plot(xcs{i + numel(Ns)},  rhos{i + numel(Ns)}, 'x-', 'DisplayName', 'WENO5');
L = legend();
L.Location = 'best';
xlim([7.3,7.5]);
ylim([0.5,4]);

xlabel('x');
ylabel('\rho');
grid on;

end
print(gcf,sprintf("SH_SUM2.png"),'-dpng','-r300')


%%%


function [rho,velo,p,name,xc] = getShuOsher(N, O2damp)


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
see = 100;



gamma = 1.4;

% [tMax, vMax, xEnd, rhoEnd, uEnd, pEnd, D] = Riemann(pnum, gamma, false);
tMax = 1.8;
vMax = 4.5657;





xs = linspace(0,10,N+1);

xc = (xs(2:end) + xs(1:end-1))/2;
xc = reshape(xc, 1, []);
h = xs(2) - xs(1);

rho0 = 1 + 0.2 * sin(5 * (xc-5));
ux0 = rho0 * 0;
p0 = ones(size(rho0));
rho0(xc<1) = 3.857143;
ux0(xc<1) = 2.629369;
p0(xc<1) = 10.333333;


ithis = 1:N;
ileft = circshift(ithis,1);
iright = circshift(ithis,-1);

u0 = f_prim2cons(rho0,ux0,p0,1.4,1);

dt = CFL * h/vMax;

ARSFix0 = ARSFix;



u = u0;

t = 0;

figure(3); clf; hold off;


% u = gpuArray(u);
% h = gpuArray(h);
% gamma = gpuArray(gamma);
% ithis = gpuArray(ithis);
% ileft = gpuArray(ileft);
% iright = gpuArray(iright);

for iter = 1:10000000
    dtc = dt;
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
    if mod(iter, see) == 0
        plot(xc,u(1,:),'-*');
        drawnow;
    end
    
end
[rho,velo,p] = f_cons2prim(u,gamma,1);

if(O2damp == 1)
    scheme_name = "TVD";
elseif(O2damp == 2)
    scheme_name = "WENO5";
else
error("O2damp bad");    
end

name = scheme_name;

end



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