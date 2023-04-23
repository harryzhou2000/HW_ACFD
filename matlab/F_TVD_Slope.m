function d = F_TVD_Slope(dL,dR)

d = (sign(dL) + sign(dR)) .* abs(dL.*dR)./ (abs(dL + dR) + 1e-300); %VL
% d = 0.5 * (sign(dL) + sign(dR)) .* min(abs(dL),abs(dR)); %MM