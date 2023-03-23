clear;
syms kappa
syms c sigma
K = 1 - c*(3/2-2*exp(-i*kappa)+1/2*exp(-2i*kappa)) + sigma * (exp(i*kappa)+exp(-i*kappa)-2);
assume(kappa,'real')
assume([c,sigma],'positive')
RK = real(K)
IK = imag(K)
K2 = simplify(expand(RK^2 + IK^2))
fK2 = matlabFunction(symfun(K2,[kappa,c,sigma]))

kaps = linspace(0,pi,1024);
cs = linspace(0,1,256);
ss = linspace(0,1,256);

[cm, km,sm] = meshgrid(cs,kaps,ss);

K2m = fK2(km,cm,sm);
K2mMax = max(K2m,[],1);
K2mMax = reshape(K2mMax, 256, 256);
cm2 = reshape(max(cm, [],1), 256, 256);
sm2 = reshape(max(sm, [],1), 256, 256);

%%
cla;
pc = pcolor(cm2,sm2,1-double(K2mMax<=2.0+10*eps));
pc.LineStyle = 'none';
pc.DisplayName = '|K|<1';
view(2)
xlabel('c');
ylabel('\sigma');
title('|K|<1');


hold on;
plot(cs, cs.^2/2,'r-', 'DisplayName', 'c^2/2', 'LineWidth',1);
plot(cs, 0.5-cs,'g-', 'DisplayName', '0.5-c', 'LineWidth',1);
legend;
colormap('gray');
box on;
set(gca,'LineWidth',2);
grid on;