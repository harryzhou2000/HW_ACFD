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

%% EPS
ks = linspace(0,pi, 10001);
epss = 0.01;
alphas = linspace(0,0.2,10001);
kss = nan * alphas;

for i = 1:numel(alphas)
    alphac = alphas(i);
    errr  = ks - fRKM(ks,alphac, 0);
    errgood = abs(errr) <= epss;
    errChange = [diff(double(errgood))==-1, false];
    kmm = min(ks(errChange));
    kss(i) = kmm;
end

figure(3);
clf; cla; hold on;
plot(alphas,kss)

alphac = 1/30 + 0.01;
kmax =  pi - acos((12*alpha - 1)/(18*alpha));
kmmax = subs(RKM,k,kmax);
s1 = symfun(kmmax - kmax - epss,alpha);
solalpha = double(solve(s1));
plot([solalpha,solalpha], [0.5,1.5],'--');
legend('x^*','\alpha_{opt}');

xlabel('\alpha'); ylabel('k^*'); grid on;

%% RE opt
ks = linspace(0,pi, 129);
alphas = [solalpha,linspace(0,0.1,6)];
clf; cla; hold on;
for alphac = alphas
    kms = fRKM(ks,alphac, 0);
    if(alphac == solalpha)
       plot(ks,kms-ks,'x-','DisplayName', sprintf('\\alpha=\\alpha_{opt}'));
    else
        plot(ks,kms-ks,'.-','DisplayName', sprintf('\\alpha=%.5g',alphac));
    end
    
end
% plot(ks,ks*0,'DisplayName' ,'Exact');
xlabel('k');
ylabel('Re(k'')-k');
legend;
grid on;
plot(ks,ks*0+epss, '-k' ,'DisplayName', sprintf('k+%g', epss));
plot(ks,ks*0-epss, '-k' ,'DisplayName', sprintf('k+%g', epss));
ylim([-epss,epss]*1.5);

