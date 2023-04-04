function [fL, fR] = F_interpi_weno5(f, eps, p)



ni = size(f,1);
nj = size(f,2);
nv = size(f,3);

i0 = 1:ni;
im1 = circshift(i0, 1);
im2 = circshift(i0, 2);
ip1 = circshift(i0,-1);
ip2 = circshift(i0,-2);

fR0 = (+ 2 * f(im2 ,:,:) - 7 * f(im1 ,:,:) + 11 * f(i0   ,:,:))/6;
fR1 = (- 1 * f(im1 ,:,:) + 5 * f(i0  ,:,:) + 2  * f(ip1  ,:,:))/6;
fR2 = (+ 2 * f(i0  ,:,:) + 5 * f(ip1 ,:,:) - 1  * f(ip2  ,:,:))/6;

fL0 = (+ 2 * f(ip2 ,:,:) - 7 * f(ip1 ,:,:) + 11 * f(i0   ,:,:))/6;
fL1 = (- 1 * f(ip1 ,:,:) + 5 * f(i0  ,:,:) + 2  * f(im1  ,:,:))/6;
fL2 = (+ 2 * f(i0  ,:,:) + 5 * f(im1 ,:,:) - 1  * f(im2  ,:,:))/6;

% wL = repmat(reshape([1 6 3]'/10, 1,1,3),ni,nj,1);
% wR = wL;
if(eps < 1e100)
    
    betaR0 = 13/12*(f(im2 ,:,:) - 2 * f(im1 ,:,:) + f(i0  ,:,:)).^2 + 1/4*(f(im2 ,:,:) - 4 * f(im1 ,:,:) + 3 * f(i0  ,:,:)).^2;
    betaR1 = 13/12*(f(im1 ,:,:) - 2 * f(i0  ,:,:) + f(ip1 ,:,:)).^2 + 1/4*(f(ip1 ,:,:) - f(im1 ,:,:)).^2;
    betaR2 = 13/12*(f(i0  ,:,:) - 2 * f(ip1 ,:,:) + f(ip2 ,:,:)).^2 + 1/4*(f(ip2 ,:,:) - 4 * f(ip1 ,:,:) + 3 * f(i0  ,:,:)).^2;
    
    betaL0 = 13/12*(f(ip2 ,:,:) - 2 * f(ip1 ,:,:) + f(i0  ,:,:)).^2 + 1/4*(f(ip2 ,:,:) - 4 * f(ip1 ,:,:) + 3 * f(i0  ,:,:)).^2;
    betaL1 = 13/12*(f(ip1 ,:,:) - 2 * f(i0  ,:,:) + f(im1 ,:,:)).^2 + 1/4*(f(im1 ,:,:) - f(ip1 ,:,:)).^2;
    betaL2 = 13/12*(f(i0  ,:,:) - 2 * f(im1 ,:,:) + f(im2 ,:,:)).^2 + 1/4*(f(im2 ,:,:) - 4 * f(im1 ,:,:) + 3 * f(i0  ,:,:)).^2;
    
    alphaR0 = 1/10 ./ (eps + betaR0).^p;
    alphaR1 = 6/10 ./ (eps + betaR1).^p;
    alphaR2 = 3/10 ./ (eps + betaR2).^p;
    alphaL0 = 1/10 ./ (eps + betaL0).^p;
    alphaL1 = 6/10 ./ (eps + betaL1).^p;
    alphaL2 = 3/10 ./ (eps + betaL2).^p;
else
    alphaR0 = 1/10;
    alphaR1 = 6/10;
    alphaR2 = 3/10;
    alphaL0 = 1/10;
    alphaL1 = 6/10;
    alphaL2 = 3/10;
end


alphaSR = alphaR0 + alphaR1 + alphaR2;
alphaSL = alphaL0 + alphaL1 + alphaL2;

omegaR0 = alphaR0 ./ alphaSR;
omegaR1 = alphaR1 ./ alphaSR;
omegaR2 = alphaR2 ./ alphaSR;
omegaL0 = alphaL0 ./ alphaSL;
omegaL1 = alphaL1 ./ alphaSL;
omegaL2 = alphaL2 ./ alphaSL;

fR = omegaR0 .* fR0 + omegaR1 .* fR1 + omegaR2 .*fR2;
fL = omegaL0 .* fL0 + omegaL1 .* fL1 + omegaL2 .*fL2;


