function F = EulerFluxARS(UL, UR, gamma, d, type,fix)


% 1: Roe



[rhoL,veloL,pL] = f_cons2prim(UL, gamma,d);
[rhoR,veloR,pR] = f_cons2prim(UR, gamma,d);
FL = invF(UL,veloL,pL);
FR = invF(UR,veloR,pR);
HL = (UL(d+2,:) + pL)./rhoL;
HR = (UR(d+2,:) + pR)./rhoR;

sqrtRhoL = sqrt(rhoL);
sqrtRhoR = sqrt(rhoR);
rhoRoe = sqrtRhoL .* sqrtRhoR;
veloRoe = (sqrtRhoL .* veloL + sqrtRhoR .* veloR) ./ (sqrtRhoL + sqrtRhoR);
HRoe = (sqrtRhoL .* HL + sqrtRhoR .* HR) ./ (sqrtRhoL + sqrtRhoR);
vsqrRoe = sum(veloRoe.^2,1);
vmRoe = sqrt(vsqrRoe);
asqrRoe = (gamma-1) * (HRoe - 0.5 * vsqrRoe);
if any(asqrRoe < 0)
    warning('imag a');
end
aRoe = sqrt(asqrRoe);
lam0 = abs(veloRoe(1,:) - aRoe);
lam123 = abs(veloRoe(1,:));
lam4 = abs(veloRoe(1,:) + aRoe);

asqrL = gamma * pL ./ rhoL;
asqrR = gamma * pR ./ rhoR;

if fix == 1
    % HY fix
    hyfixT = (vmRoe + aRoe) * 0.05;
    lam0 = fix_hy(hyfixT, lam0);
    lam4 = fix_hy(hyfixT, lam4);
    %
elseif fix == 0
elseif fix == 2
    lamMax = max(abs(veloL(1,:)) + sqrt(asqrL),abs(veloR(1,:)) + sqrt(asqrR));
    fixT = lamMax  * 0.4;
    fixed = lam123 < fixT;
    lam123(fixed) = lamMax(fixed);
    lam0(fixed) = lamMax(fixed);
    lam4(fixed) = lamMax(fixed);
    
    hyfixT = (vmRoe + aRoe) * 0.05;
    lam0 = fix_hy(hyfixT, lam0);
    lam4 = fix_hy(hyfixT, lam4);
else
    error('input');
end


incU = UR - UL;


if type == 1 %Roe
    
    alpha23 = incU(3:1+d,:) - veloRoe(2:d,:) .* incU(1,:);
    incU4b = incU(d+2,:) - dot(alpha23, veloRoe(2:d,:), 1);
    alpha1 = (gamma-1)./asqrRoe .* (...
        incU(1,:).*(HRoe - veloRoe(1,:).^2) + veloRoe(1,:) .* incU(2,:) -incU4b);
    alpha0 = (incU(1,:) .* (veloRoe(1,:) + aRoe) - incU(2,:) - aRoe .* alpha1) ./ (2 * aRoe);
    alpha4 = incU(1,:) - (alpha0 + alpha1);
    
    dF0 = alpha0.*lam0 + alpha1.*lam123 + alpha4.*lam4;
    dF1 = (veloRoe(1,:) - aRoe).*alpha0.*lam0 + veloRoe(1,:).*alpha1.*lam123 + (veloRoe(1,:) + aRoe).*alpha4.*lam4;
    dF23 = veloRoe(2:d,:).*alpha0.*lam0 + veloRoe(2:d,:).*alpha1.*lam123 + veloRoe(2:d,:).*alpha4.*lam4 + ...
        alpha23.*lam123;
    dF4 = (HRoe - veloRoe(1,:).*aRoe).*alpha0.*lam0 + (HRoe + veloRoe(1,:).*aRoe).*alpha4.*lam4 + ...
        0.5 * vsqrRoe .* alpha1.*lam123 + dot(veloRoe(2:d,:),alpha23,1).*lam123;
    
    incF = [dF0;dF1;dF23;dF4];
    
    F = (FL+FR - incF)/2;
    return;
    
elseif type == 2 || type == 3 %HLL
    
%     eta2 = 0.5 * (sqrtRhoL .* sqrtRhoR) ./ (sqrtRhoL + sqrtRhoR).^2;
%     dsqr = (asqrL.*sqrtRhoL + asqrR.*sqrtRhoR)  ./ (sqrtRhoL + sqrtRhoR) + eta2.*(veloL(1,:) - veloR(1,:)).^2;
%     if any(dsqr < 0)
%         error('imag d');
%     end
%     aHLL = sqrt(dsqr);
%     SL = veloRoe(1,:) - aHLL;
%     SR = veloRoe(1,:) + aHLL;
%     SS = veloRoe(1,:);
    pS = 0.5 * (pL + pR) - 0.5 * (veloR(1,:) - veloL(1,:)) .* rhoRoe .* aRoe;
    pS = max(0,pS);
    qL = HLLC_q(pL,pS);
    qR = HLLC_q(pR,pS);
    SL = veloL(1,:) - sqrt(asqrL).*qL;
    SR = veloR(1,:) + sqrt(asqrR).*qR;
    SS = pR - pL + rhoL .*veloL(1,:).*(SL - veloL(1,:)) - rhoR .*veloR(1,:).*(SR - veloR(1,:));
    SS = SS./(rhoL.*(SL-veloL(1,:)) - rhoR.*(SR-veloR(1,:)));
    
    
    
    
    if type == 2 %HLL
        r_L = double(0<=SL);
        r_R = double(SR<=0);
        r_M = double(~(r_L | r_R));
        F = FL.*r_L + FR.*r_R + (SR.*FL-SL.*FR+SL.*SR.*(UR-UL))./(SR-SL).*r_M;
        return;
    elseif type == 3 %HLLC
        %         pSL = HLLC_PS(veloL, pL, rhoL, SL, SS);
        %         pSR = HLLC_PS(veloR, pR, rhoR, SR, SS);
        FSL = HLLC_FS(UL,veloL, rhoL, pL, FL, SL, SS);
        FSR = HLLC_FS(UR,veloR, rhoR, pR, FR, SR, SS);
        FSM = (FSL + FSR)/2;
        r_L = double(0<=SL);
        r_R = double(SR<=0);
        r_SL = double(SL<0 & 0<SS);
        r_SM = double(SS==0);
        r_SR = double(SS<0 & 0<SR);
        
        F = FL.*r_L + FR.*r_R + FSL.*r_SL + FSR.*r_SR + FSM.*r_SM;
        return;
    else
         error('input');
    end
else
     error('input');
end





    function v = fix_hy(hyT, v)
        ifFix = double(v < hyT);
        v = ifFix .* (v.^2 + hyT.^2)./(2 * hyT) + (1-ifFix) .* v;
        
    end

    function F = invF(U, velo, p)
        F = U.*velo(1,:);
        F(2,:) = F(2,:) + p;
        F(d+2,:) = F(d+2,:) + p.*velo(1,:);
    end

    function pS = HLLC_pS(velo, p, rho, S, SS)
        pS = p + rho.*(S - velo(1,:)).*(SS - velo(1,:));
    end

    function FS = HLLC_FS(U,velo, rho, p, F, S, SS)
        FS = SS.*(S.*U - F);
        iA = S.*(p + rho.*(S - velo(1,:)).*(SS - velo(1,:)));
        FS(2,:) = FS(2,:) + iA;
        FS(2+d,:) = FS(2+d,:) + iA.*SS;
        FS = FS./(S - SS);
    end

    function q = HLLC_q(p, pS)
       noshock = pS <= p;
       qshock = sqrt(1 + (gamma+1)/2/gamma *(pS./p - 1) );
       qshock(noshock) = 1;
       q = qshock;
    end


end


