function F = RoeSolver2(uL, uR, gamma)

[rhoL, UxL,UyL,UsqrL, EL,pL,aL,HL] = F_uExpand_V1(uL, gamma);
[rhoR, UxR,UyR,UsqrR, ER,pR,aR,HR] = F_uExpand_V1(uR, gamma);

rhoLsqrt = sqrt(rhoL);
rhoRsqrt = sqrt(rhoR);
rhoLRsqrtSum = rhoLsqrt + rhoRsqrt;

UxRoe = (rhoLsqrt .* UxL + rhoRsqrt.* UxR)./rhoLRsqrtSum;
UyRoe = (rhoLsqrt .* UyL + rhoRsqrt.* UyR)./rhoLRsqrtSum;
HRoe  = (rhoLsqrt .* HL + rhoRsqrt.* HR)./rhoLRsqrtSum;
UsqrRoe = UxRoe.^2 + UyRoe.^2;
asqrRoe = (gamma-1) * (HRoe - 0.5 * UsqrRoe);
if(sum(asqrRoe<0)>0)
   error('imag a');
end
aRoe = sqrt(asqrRoe);

lam1Roe = UxRoe - aRoe;
lam2Roe = UxRoe;
lam3Roe = UxRoe;
lam5Roe = UxRoe + aRoe;

%HY
deltaEF = 0.05 * (sqrt(UsqrRoe) + aRoe);
lam1Fix = abs(lam1Roe) < deltaEF;
lam5Fix = abs(lam5Roe) < deltaEF;
lam1Roe(lam1Fix) = sign(lam1Roe(lam1Fix)).*(lam1Roe(lam1Fix).^2+deltaEF(lam1Fix).^2)./deltaEF(lam1Fix)/2;
lam5Roe(lam5Fix) = sign(lam5Roe(lam5Fix)).*(lam5Roe(lam5Fix).^2+deltaEF(lam5Fix).^2)./deltaEF(lam5Fix)/2;
%HY

rev1 = [ones(size(lam1Roe)) ; lam1Roe             ; UyRoe               ; HRoe - UxRoe.*aRoe];
rev2 = [ones(size(lam1Roe)) ; UxRoe               ; UyRoe               ; 0.5 * UsqrRoe     ];
rev3 = [zeros(size(lam1Roe)); zeros(size(lam1Roe)); ones(size(lam1Roe)) ; UyRoe             ];
rev5 = [ones(size(lam1Roe)) ; lam5Roe             ; UyRoe               ; HRoe + UxRoe.*aRoe];

uinc = uR - uL;

alpha3 = uinc(3,:) - UyRoe .* uinc(1,:);
incU5c = uinc(4,:) - ( uinc(3,:) - UyRoe .* uinc(1,:) ) .* UyRoe;
alpha2 = (gamma - 1)./(asqrRoe) .* (uinc(1,:) .* (HRoe - UxRoe.^2) + UxRoe .* uinc(2,:) - incU5c);
alpha1 = 0.5./aRoe .* ( uinc(1,:).*lam5Roe - uinc(2,:) - aRoe.*alpha2 );
alpha5 = uinc(1,:) - (alpha1 + alpha2);



% Pike
% incP = pR - pL;
% rhoRoe = sqrt(rhoR.*rhoL);
% incUx = UxR - UxL;
% incUy = UyR - UyL;
% 
% alpha1 = 0.5./aRoe .* (incP - rhoRoe.*aRoe.*incUx);
% alpha2 = uinc(1,:) - incP./aRoe.^2;
% alpha3 = rhoRoe .* incUy;
% alpha5 = 0.5./aRoe .* (incP + rhoRoe.*aRoe.*incUx);
% Pike


FL = F_u2F(uL,gamma);
FR = F_u2F(uR,gamma);


F = 0.5 * (FL + FR)  - 0.5 * (...
    alpha1 .* abs(lam1Roe) .* rev1 + ...
    alpha2 .* abs(lam2Roe) .* rev2 + ...
    alpha3 .* abs(lam3Roe) .* rev3 + ...
    alpha5 .* abs(lam5Roe) .* rev5);








