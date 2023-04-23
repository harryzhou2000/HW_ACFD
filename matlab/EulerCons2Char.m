function incW = EulerCons2Char(veloRoe,vsqrRoe,aRoe,asqrRoe,HRoe, incU, d, gamma)




alpha23 = incU(3:1+d,:) - veloRoe(2:d,:) .* incU(1,:);
incU4b = incU(d+2,:) - dot(alpha23, veloRoe(2:d,:), 1);
alpha1 = (gamma-1)./asqrRoe .* (...
    incU(1,:).*(HRoe - veloRoe(1,:).^2) + veloRoe(1,:) .* incU(2,:) -incU4b);
alpha0 = (incU(1,:) .* (veloRoe(1,:) + aRoe) - incU(2,:) - aRoe .* alpha1) ./ (2 * aRoe);
alpha4 = incU(1,:) - (alpha0 + alpha1);

incW = [alpha0;alpha1;alpha23;alpha4];