function [fL, fR] = F_interpi_weno5_char1(fLi,fRi, eps, p, mapping)
if nargin < 5
    mapping = 0;
end



ni = size(fLi,1);
nj = size(fLi,2);
nv = size(fLi,3);

i0 = 1:ni;
im1 = circshift(i0, 1);
im2 = circshift(i0, 2);
ip1 = circshift(i0,-1);
ip2 = circshift(i0,-2);

fm2R = fRi(im2 ,:,:);
fm1R = fRi(im1 ,:,:);
fccR = fRi(i0  ,:,:);
fp1R = fRi(ip1 ,:,:);
fp2R = fRi(ip2 ,:,:);

fm2L = fLi(im2 ,:,:);
fm1L = fLi(im1 ,:,:);
fccL = fLi(i0  ,:,:);
fp1L = fLi(ip1 ,:,:);
fp2L = fLi(ip2 ,:,:);





fR0 = (+ 2 * fm2R - 7 * fm1R + 11 * fccR)/6;
fR1 = (- 1 * fm1R + 5 * fccR + 2  * fp1R)/6;
fR2 = (+ 2 * fccR + 5 * fp1R - 1  * fp2R)/6;

fL0 = (+ 2 * fp2L - 7 * fp1L + 11 * fccL)/6;
fL1 = (- 1 * fp1L + 5 * fccL + 2  * fm1L)/6;
fL2 = (+ 2 * fccL + 5 * fm1L - 1  * fm2L)/6;

% wL = repmat(reshape([1 6 3]'/10, 1,1,3),ni,nj,1);
% wR = wL;
if(eps < 1e100)
    
    betaR0 = 13/12*(fm2R - 2 * fm1R + fccR).^2 + 1/4*(fm2R - 4 * fm1R + 3 * fccR).^2;
    betaR1 = 13/12*(fm1R - 2 * fccR + fp1R).^2 + 1/4*(fp1R - fm1R).^2;
    betaR2 = 13/12*(fccR - 2 * fp1R + fp2R).^2 + 1/4*(fp2R - 4 * fp1R + 3 * fccR).^2;
    
    betaL0 = 13/12*(fp2L - 2 * fp1L + fccL).^2 + 1/4*(fp2L - 4 * fp1L + 3 * fccL).^2;
    betaL1 = 13/12*(fp1L - 2 * fccL + fm1L).^2 + 1/4*(fm1L - fp1L).^2;
    betaL2 = 13/12*(fccL - 2 * fm1L + fm2L).^2 + 1/4*(fm2L - 4 * fm1L + 3 * fccL).^2;
    % betaL0 = betaR2;
    % betaL1 = betaR1;
    % betaL2 = betaR0;
    
    if(p < 1e100)
        alphaR0 = 1/10 ./ (eps + betaR0).^p;
        alphaR1 = 6/10 ./ (eps + betaR1).^p;
        alphaR2 = 3/10 ./ (eps + betaR2).^p;
        alphaL0 = 1/10 ./ (eps + betaL0).^p;
        alphaL1 = 6/10 ./ (eps + betaL1).^p;
        alphaL2 = 3/10 ./ (eps + betaL2).^p;
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         alphaR0 = double(betaR0 <= betaR1 & betaR0 <= betaR2);
        %         alphaR1 = double((betaR1 <= betaR0 & betaR1 <= betaR2) & (~alphaR0));
        %         alphaR2 = double((betaR2 <= betaR0 & betaR2 <= betaR1) & (~alphaR0) & (~alphaR1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% %ENO3 vanilla
        
        diffL = abs(f(im1 ,:,:) - f(i0  ,:,:));
        diffR = abs(f(ip1 ,:,:) - f(i0  ,:,:));
        useL = diffL < diffR;
        diffLL = abs(f(im2 ,:,:) + f(i0 ,:,:) - 2 * f(im1 ,:,:));
        diffRR = abs(f(ip2 ,:,:) + f(i0 ,:,:) - 2 * f(ip1 ,:,:));
        diffMM = abs(f(ip1 ,:,:) + f(im1 ,:,:) - 2 * f(i0  ,:,:));
        
        useLL =  useL  & (diffLL <  diffMM);
        useMM = (useL  & (diffLL >= diffMM)) | ((~useL) & (diffRR >= diffMM));
        useRR =(~useL) & (diffRR <  diffMM);
        
        alphaR0 = double(useLL);
        alphaR1 = double(useMM);
        alphaR2 = double(useRR);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        alphaL0 = alphaR2;
        alphaL1 = alphaR1;
        alphaL2 = alphaR0;
        
        if any(alphaL0 + alphaL1 + alphaL2 - 1)
            error('bad eno');
        end
        
    end
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

if(eps<100 && mapping)
    alphaR0 = Henrick_mapping(omegaR0, 1/10);
    alphaR1 = Henrick_mapping(omegaR1, 6/10);
    alphaR2 = Henrick_mapping(omegaR2, 3/10);
    alphaL0 = Henrick_mapping(omegaL0, 1/10);
    alphaL1 = Henrick_mapping(omegaL1, 6/10);
    alphaL2 = Henrick_mapping(omegaL2, 3/10);
    
    alphaSR = alphaR0 + alphaR1 + alphaR2;
    alphaSL = alphaL0 + alphaL1 + alphaL2;
    
    omegaR0 = alphaR0 ./ alphaSR;
    omegaR1 = alphaR1 ./ alphaSR;
    omegaR2 = alphaR2 ./ alphaSR;
    omegaL0 = alphaL0 ./ alphaSL;
    omegaL1 = alphaL1 ./ alphaSL;
    omegaL2 = alphaL2 ./ alphaSL;
    
end

fR = omegaR0 .* fR0 + omegaR1 .* fR1 + omegaR2 .*fR2;
fL = omegaL0 .* fL0 + omegaL1 .* fL1 + omegaL2 .*fL2;

    function am = Henrick_mapping(w,wi)
        am = w.*(wi + wi^2 - 3*wi*w+w.^2)./(wi^2 +w * (1-2*wi));
    end

end


