function incU = EulerChar2Cons(incW, lam0, lam123, lam4, veloRoe,vsqrRoe,aRoe,asqrRoe,HRoe,d, gamma)



alpha0 = incW(1,:);
alpha1 = incW(2,:);
alpha23 = incW(3:1+d,:);
alpha4 = incW(2+d,:);

dF0 = alpha0.*lam0 + alpha1.*lam123 + alpha4.*lam4;
dF1 = (veloRoe(1,:) - aRoe).*alpha0.*lam0 + veloRoe(1,:).*alpha1.*lam123 + (veloRoe(1,:) + aRoe).*alpha4.*lam4;
dF23 = veloRoe(2:d,:).*alpha0.*lam0 + veloRoe(2:d,:).*alpha1.*lam123 + veloRoe(2:d,:).*alpha4.*lam4 + ...
    alpha23.*lam123;
dF4 = (HRoe - veloRoe(1,:).*aRoe).*alpha0.*lam0 + (HRoe + veloRoe(1,:).*aRoe).*alpha4.*lam4 + ...
    0.5 * vsqrRoe .* alpha1.*lam123 + dot(veloRoe(2:d,:),alpha23,1).*lam123;

incU = [dF0;dF1;dF23;dF4];