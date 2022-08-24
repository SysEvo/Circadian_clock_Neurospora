clear; close all; clc;

x0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; 
tM = 0:0.05:200;   
[t,x] = ode45(@modelo_complejo,tM,x0);
x1 = 24;

% Plot ODE results
[pks1,locs1] = findpeaks(x(:,8));
[pks2,locs2] = findpeaks(x(:,12));

plot(t,x(:,[8,12]),'-', t(locs1),pks1, "bo", t(locs2),pks2, "bo",...
    120, linspace(2,6,20), 'r.', ...
    144, linspace(2,6,20), 'r.', ...
    168, linspace(2,6,20), 'r.', ...
    192, linspace(2,6,20), 'r.')
text(t(locs1)+.01,pks1+.05,num2str(t(locs1)),"FontSize",6)
text(t(locs2)+.05,pks2+0.05,num2str(t(locs2)),"FontSize",5)
title("LD 2:2")
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100,200])
ylim([3.2,4.2])
legend(["frq","FRQn"])