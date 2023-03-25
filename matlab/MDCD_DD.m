syms k alpha beta
assume([k,alpha,beta],'real');

km = 1/1i * (...
    (-alpha - beta)/2 * exp(-3i * k) + (+alpha - beta)/2 * exp( 3i * k) + ...
    (24*alpha + 36*beta + 1)/12 * exp(-2i * k) + (-24*alpha + 36*beta - 1)/12 * exp(2i * k) + ...
    (-15*alpha - 45*beta -4)/6 * exp(-1i * k) + ( 15*alpha - 45*beta +4)/6 * exp( 1i * k) + ...
    10*beta);

RKM = simplify(expand(real(km)));
IKM = simplify(expand(imag(km)));

latex(RKM)


fRKM = matlabFunction(symfun(RKM,[k,alpha,beta]));
fIKM = matlabFunction(symfun(IKM,[k,alpha,beta]));
ks = linspace(0,pi, 129);
%% RE
alphas = linspace(0,0.1,6);
clf; cla; hold on;
for alphac = alphas
    kms = fRKM(ks,alphac, 0);
    plot(ks,kms,'.-','DisplayName', sprintf('\\alpha=%g',alphac));
end
plot(ks,ks,'DisplayName' ,'Exact');
xlabel('k');
ylabel('Re(k'')');
legend;
grid on;




%% IM
betas = linspace(0,0.1,6);
clf; cla; hold on;
for betac = betas
    kms = fIKM(ks,0, betac);
    plot(ks,kms,'.-','DisplayName', sprintf('\\beta=%g',betac));
end

plot(ks,ks*0,'DisplayName' ,'Exact');
xlabel('k');
ylabel('Im(k'')');
legend;
grid on;