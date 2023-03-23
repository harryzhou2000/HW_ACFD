
clear;

N = 128;
xs = linspace(0,1-1/N,N);

this = 1:N;
right = circshift(this,-1);
left = circshift(this,1);
left2 = circshift(this,2);

cmax = (-1 + sqrt(1 + 1))/1;

c = cmax;
sigma = c^2/2;
us = sin(xs *2* pi) + 0.4 * sin(xs * 10* pi);


for iter = 1:1000
    duC = - c * (1.5*us - 2 * us(left) + 0.5 * us(left2));
    duD = sigma * (us(left) + us(right) - 2*us);
    us = us + duC + duD;
    
    if(mod(iter, 10) ==0)
        cla;
        clf;
        plot(xs,us);
        drawnow;
    end
end

