
clear;

N = 100;
Tmax = .1;
a = 1;
xs = linspace(-0.5,0.5-1/N,N);

this = 1:N;
right = circshift(this,-1);
left = circshift(this,1);
left2 = circshift(this,2);

cmax = (-1 + sqrt(1 + 1))/1;

c = 0.5;
sigma = c^2/2;
% us = sin(xs *2* pi) + 0.4 * sin(xs * 10* pi);
us0  =double(xs>=-0.25 & xs<=0.25);

dt = c * (xs(2)-xs(1)) / a;
cla;
clf; hold on;

%% UW1
us = us0;
for iter = 1:round(Tmax/dt)
    duC = - c * (us-us(left));
    us = us + duC;
end
plot(xs,us, '.-','DisplayName', 'UW1');
drawnow;
%% LW
us = us0;
for iter = 1:round(Tmax/dt)
    duC = - c * (0.5*us(right)-0.5 * us(left));
    duD = sigma * (us(left) + us(right) - 2*us(this));
    us = us + duC + duD;
end
plot(xs,us, 'o-', 'DisplayName', 'LW');
drawnow;
%% WB
us = us0;
for iter = 1:round(Tmax/dt)
    duC = - c * (1.5*us - 2 * us(left) + 0.5 * us(left2));
    duD = sigma * (us(left2) + us(this) - 2*us(left));
    us = us + duC + duD;
end
plot(xs,us, 'd-', 'DisplayName', 'WB');
drawnow;
%% Ana
xsnew = xs - a * Tmax;
xsnew = mod(xsnew + 0.5, 1) - 0.5;
ua = (xsnew) >= -0.25 & (xsnew) <= 0.25;
ua = double(ua);
plot(xs,ua, '-','DisplayName', 'analytical');
drawnow;
%%
L = legend;
L.Location = 'northeastoutside';
grid on;
xlabel('x');
ylabel('u');
title(sprintf('t = %g', Tmax));
