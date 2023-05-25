function F = EulerFluxFDi(UL, UR, gamma, d, fix)


% "UR" is actually UC

% 1: Roe
sz = size(UL);
shape0 = sz;
% ni = sz(2);
% i0 = 1:ni;
% im1 = circshift(i0, 1);
% im2 = circshift(i0, 2);
% ip1 = circshift(i0,-1);
% ip2 = circshift(i0,-2);

% UL = reshape(UL, sz(1),[]);
% UR = reshape(UR, sz(1),[]);

[rhoL,veloL,pL] = f_cons2prim(UL, gamma,d);
[rhoR,veloR,pR] = f_cons2prim(UR, gamma,d);
% FL = invF(UL,veloL,pL);
% FR = invF(UR,veloR,pR);
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


% incU = UR - UL;




% alpha23 = incU(3:1+d,:) - veloRoe(2:d,:) .* incU(1,:);
% incU4b = incU(d+2,:) - dot(alpha23, veloRoe(2:d,:), 1);
% alpha1 = (gamma-1)./asqrRoe .* (...
%     incU(1,:).*(HRoe - veloRoe(1,:).^2) + veloRoe(1,:) .* incU(2,:) -incU4b);
% alpha0 = (incU(1,:) .* (veloRoe(1,:) + aRoe) - incU(2,:) - aRoe .* alpha1) ./ (2 * aRoe);
% alpha4 = incU(1,:) - (alpha0 + alpha1);


[aR, aL] = F_interpi_weno5_char(...
    permute(reshape(UR,shape0(1),1,shape0(2)),[3,2,1]),...
    permute(reshape(UL,shape0(1),1,shape0(2)),[3,2,1]),...
    @(U) permute(reshape(getAlphas(reshape(permute(U,[3,2,1]),shape0)),shape0(1),1,shape0(2)),...
    [3,2,1]), 1e-6, 2, 1);
aR = reshape(permute(aR,[3,2,1]),shape0);
aL = reshape(permute(aL,[3,2,1]),shape0);

% fprintf("%e\n",norm(getFs(alphasL,1,1,1)- UL))

ULRec = getFs(aL, 1,1,1);
URRec = getFs(aR, 1,1,1);
[rhoL,veloL,pL] = f_cons2prim(ULRec, gamma,d);
[rhoR,veloR,pR] = f_cons2prim(URRec, gamma,d);
FL = invF(ULRec,veloL,pL);
FR = invF(URRec,veloR,pR);
incF = getFs(aR - aL, lam0, lam123, lam4);

F = (FL+FR - incF)/2;

%
% [aR, aL] = F_interpi_weno5_char1(...
%     permute(reshape(UR,shape0(1),1,shape0(2)),[3,2,1]),...
%     permute(reshape(UL,shape0(1),1,shape0(2)),[3,2,1]), 0.001, 2, 0);
% URRec = reshape(permute(aR,[3,2,1]),shape0);
% ULRec = reshape(permute(aL,[3,2,1]),shape0);
%
% F = EulerFluxARS(ULRec,URRec,gamma,d,1,fix);
return;






    function alphas = getAlphas(incU)
        calpha23 = incU(3:1+d,:) - veloRoe(2:d,:) .* incU(1,:);
        cincU4b = incU(d+2,:) - dot(calpha23, veloRoe(2:d,:), 1);
        calpha1 = (gamma-1)./asqrRoe .* (...
            incU(1,:).*(HRoe - veloRoe(1,:).^2) + veloRoe(1,:) .* incU(2,:) - cincU4b);
        calpha0 = (incU(1,:) .* (veloRoe(1,:) + aRoe) - incU(2,:) - aRoe .* calpha1) ./ (2 * aRoe);
        calpha4 = incU(1,:) - (calpha0 + calpha1);
        alphas = [calpha0;calpha1;calpha23;calpha4];
    end

    function incF = getFs(alphas,lam0, lam123, lam4)
        alpha0 = alphas(1,:);
        alpha1 = alphas(2,:);
        alpha23 = alphas(3:1+d,:);
        alpha4 = alphas(2+d,:);
        dF0 = alpha0.*lam0 + alpha1.*lam123 + alpha4.*lam4;
        dF1 = (veloRoe(1,:) - aRoe).*alpha0.*lam0 + veloRoe(1,:).*alpha1.*lam123 + (veloRoe(1,:) + aRoe).*alpha4.*lam4;
        dF23 = veloRoe(2:d,:).*alpha0.*lam0 + veloRoe(2:d,:).*alpha1.*lam123 + veloRoe(2:d,:).*alpha4.*lam4 + ...
            alpha23.*lam123;
        dF4 = (HRoe - veloRoe(1,:).*aRoe).*alpha0.*lam0 + (HRoe + veloRoe(1,:).*aRoe).*alpha4.*lam4 + ...
            0.5 * vsqrRoe .* alpha1.*lam123 + dot(veloRoe(2:d,:),alpha23,1).*lam123;
        
        incF = [dF0;dF1;dF23;dF4];
        
        
        
        
        
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


