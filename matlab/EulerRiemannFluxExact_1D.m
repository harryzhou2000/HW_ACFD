function [F, D] = EulerRiemannFluxExact_1D(UL, UR, gamma, tol)
if(nargin < 4)
    tol = 1e-10;
end
d = 1;
N = size(UL,2);


[rhoL,veloL,pL] = f_cons2prim(UL,gamma,d);
[rhoR,veloR,pR] = f_cons2prim(UR,gamma,d);
aL = sqrt(gamma * pL ./ rhoL);
aR = sqrt(gamma * pR ./ rhoR);
rveloL = veloL.*rhoL;
rveloR = veloR.*rhoR;


[AL, BL] = eulerRiemannConstants(pL, rhoL, gamma);
[AR, BR] = eulerRiemannConstants(pR, rhoR, gamma);



VL = veloL(1,:);
VR = veloR(1,:);
dV = VR - VL;
pM = (pL + pR) / 2;


pm = nan(1,N);
options = optimoptions('fsolve','Display','off','MaxIterations', 30,...
    'OptimalityTolerance', tol, 'FunctionTolerance', tol);
for i = 1:N
    options.TypicalX = pM(:, i);
    pm(i) = fsolve(@(pc) ...
        f_(pc, AL(:,i), BL(:,i), aL(:,i), pL(:,i), gamma) + ...
        f_(pc, AR(:,i), BR(:,i), aR(:,i), pR(:,i), gamma) + ...
        dV(:,i)...
        , pM(:, i), options);
    if(abs(imag(pm(i))) > tol * pM(i) )
        error('bad sol');
    end
end
pm = real(pm);

Vm = (VL + VR) / 2 + (f_(pm, AR, BR, aR, pR, gamma) - f_(pm, AL, BL, aL, pL, gamma))/2;

rhomL = frho(pm, pL, rhoL, gamma);
rhomR = frho(pm, pR, rhoR, gamma);
amL = sqrt(gamma * pm./rhomL);
amR = sqrt(gamma * pm./rhomR);
rvelomL = frho(pm, pL, rveloL, gamma);
rvelomR = frho(pm, pR, rveloR, gamma);

srLL = VL - aL;
srLR = Vm - amL;
srLL = min(srLL,srLR);
srRL = Vm + amR;
srRR = VR + aR;
srRR = max(srRL,srRR);

sSL = -fshockSpeedR(pm,-Vm,amL,pL, gamma);
sSR =  fshockSpeedR(pm, Vm,amR,pR, gamma);

sLL = (pm > pL) .* sSL + (pm <= pL) .* srLL;
sLR = (pm > pL) .* sSL + (pm <= pL) .* srLR;
sRL = (pm > pR) .* sSR + (pm <= pR) .* srRL;
sRR = (pm > pR) .* sSR + (pm <= pR) .* srRR;

VrL = (gamma-1)/(gamma+1) * VL + 2/(gamma + 1) * aL;
VrR = (gamma-1)/(gamma+1) * VR - 2/(gamma + 1) * aR;
rhorL = abs(VrL ./ aL).^(2/(gamma-1)) .* rhoL; % using abs(Vr) == a in rF
rhorR = abs(VrR ./ aR).^(2/(gamma-1)) .* rhoR;
prL = (rhorL ./ rhoL).^gamma .* pL;
prR = (rhorR ./ rhoR).^gamma .* pR;

rvelorL = (VrL ./ aL).^(2/(gamma-1)) .* rveloL;
rvelorR = (VrR ./ aR).^(2/(gamma-1)) .* rveloR;

good = sLL <= sLR & sLR <= Vm & Vm <= sRL & sRL <= sRR;
if(~all(good))
   error('bad velocities'); 
end

r_left = double(sLL > 0);
r_r_left = double(sLL <= 0 & sLR > 0);
r_m_left = double(sLR <= 0 & Vm > 0);
r_m_m = double(Vm == 0);
r_m_right = double(Vm < 0 & sRL >= 0);
r_r_right = double(sRL < 0 & sRR >= 0);
r_right = double(sRR < 0);

rhoF = r_left.*rhoL + r_r_left.*rhorL + ...
    r_m_left.*rhomL + r_m_m.*(rhomL + rhomR)/2 + r_m_right.*rhomR + ...
    r_r_right.*rhorR + r_right.*rhoR;

VF = r_left.*VL + r_r_left.*VrL + ...
    r_m_left.*Vm + r_m_m.*(Vm + Vm)/2 + r_m_right.*Vm + ...
    r_r_right.*VrR + r_right.*VR;

pF = r_left.*pL + r_r_left.*prL + ...
    r_m_left.*pm + r_m_m.*(pm + pm)/2 + r_m_right.*pm + ...
    r_r_right.*prR + r_right.*pR;

rveloF = r_left.*rveloL + r_r_left.*rvelorL + ...
    r_m_left.*rvelomL + r_m_m.*(rvelomL + rvelomR)/2 + r_m_right.*rvelomR + ...
    r_r_right.*rvelorR + r_right.*rveloR;

veloF = rveloF ./ rhoF;
veloF(1,:) = VF;

EkF = sum(veloF.^2,1) * 0.5 .* rhoF;
EF = EkF + pF / (gamma -1);

F = nan(d+2,N);
F(1,:) = rhoF .* VF;
F(2:1+d,:) = veloF .* rhoF .*VF;
F(2,:) = F(2,:) + pF;
F(end,:) = EF .* VF + pF .* VF;



D.rhoL = rhoL;
D.rhoR = rhoR;
D.VL = VL;
D.VR = VR;
D.rhomL = rhomL;
D.rhomR = rhomR;
D.pm = pm;
D.Vm = Vm;
D.sLL = sLL;
D.sLR = sLR;
D.sRL = sRL;
D.sRR = sRR;
D.pL = pL;
D.pR = pR;
D.aL = aL;
D.aR = aR;









    

    function [A,B] = eulerRiemannConstants(p, rho, gamma)
        A = 2./((gamma+1) * rho);
        B = (gamma-1)/(gamma+1) * p;
    end

    function f = f_(pm, A, B, a, p, gamma)
        fshock = (pm - p) .* (A ./ (pm + B)).^0.5;
        powerR = (gamma-1)/(2*gamma);
        fraref = 2*a/(gamma-1) .* ((pm./p).^powerR - 1);
        fshock(isnan(fshock)) = 0;
        fraref(isnan(fraref)) = 0;
        
        f = (pm > p) .* fshock + (pm <= p) .*fraref;
        
    end

    function rhom = frho(pm, p, rho, gamma)
        gammaR = (gamma+1)/(gamma-1);
        rhoShock = (p./pm + gammaR) ./ (p./pm * gammaR + 1) .* rho;
        rhoRareF = (pm./p).^(1/gamma) .* rho;
        
        rhom = (pm > p) .* rhoShock + (pm <= p) .* rhoRareF;
        
    end

    function phi = fshockSpeedR(pm,Vm,amR,pR,gamma)
       phi = Vm + amR .* sqrt((gamma+1)/2/gamma * pR./pm + (gamma-1)/2/gamma);
        
    end


end