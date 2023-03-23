
%%
N= 100;
Tend = 0.2;
cla; hold on;
for sigma = [0.1,0.5,1]
[us,xs,t, dt] = getHT1D_FTCS(sigma, N, Tend);assert((t-Tend)<1e-10);
plot(xs,us,'DisplayName',"\sigma="+string(sigma))
end
ylim([-1,1] * 4e-4);
legend
xlabel('x');ylabel('u');title("T="+string(Tend));

%%
sigma = 0.1;
Tend = 0.1;
N = 10;
[us,xs,t, dt] = getHT1D_FTCS(sigma, N, Tend);assert((t-Tend)<1e-10);
xsold = xs;usold = us;

finingErrs = [];
dxs = [];

for i = 1:4
    N = N * 2;
[us,xs,t, dt] = getHT1D_FTCS(sigma, N, Tend);assert((t-Tend)<1e-10);
err = norm(usold - us(1:2:end),2)/norm(usold,2);
finingErrs(i) = err;
dxs(i) = xs(2) - xs(1);


xsold = xs;usold = us;

end
cla; hold on;
plot(dxs,finingErrs, 'o','DisplayName','ErrR')
plot(dxs, finingErrs(1).*(dxs./dxs(1) * 1.2).^2, '--','DisplayName', 'O2');
set(gca,'XScale','log');
set(gca,'YScale','log');
legend
grid on;
xlabel('hx');ylabel('Error');
xticks(sort(dxs))
xtickformat('%.2e')

%%
function [us,xs, t,dt] = getHT1D_FTCS(sigma, N, tend)

gamma =1;
xs = linspace(0,1,N+1);
dx = xs(2) -xs(1);
dt = sigma * dx^2 / gamma;

ithis = 2:N;
ile = mod(ithis-1 - 1, N + 1) + 1;
iri = mod(ithis+1 - 1, N + 1) + 1;

u0 = sin(2*pi*xs) + 0.1*sin(20*pi*xs);
fL = @(t) 0;
fR = @(t) 0;

niter = round(tend/dt);
t = 0;
us = u0;
for i = 1:niter
    us(1) = fL(t);
    us(end) = fR(t);
    unew = us;
    unew(ithis) = us(ithis) + dt * (us(ile) + us(iri) - 2 * us(ithis))/dx^2;
    us = unew;
    t = t+ dt;
end

end



