clear; close all; clc;

x0 = [0;0;0;0;0;0;0;0;0;0;0];
%tM = linspace(0,200);
tM = 0:0.05:200;
[t,x] = ode45(@modelo_simple,tM,x0);

%Maximos
[pks1,locs1] = findpeaks(x(:,5),"MinPeakDistance",10,"MinPeakHeight",0.165);
[pks2,locs2] = findpeaks(x(:,8));

% Plot ODE results

%ax1=nexttile;
%plot(ax1,t,x(:,5),'-', t(locs1),pks1, "bo",...
%    120, linspace(0,.2,20), 'r.', ...
%    144, linspace(0,0.2,20), 'r.', ...
%    168, linspace(0,0.2,20), 'r.', ...
%    192, linspace(0,0.2,20), 'r.')
%text(ax1,t(locs1)+.01,pks1+.01,num2str(t(locs1)))
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([100 200])
%ylim([0.13,0.2])
%legend("frq")

%ax2=nexttile;
%plot(ax2,t,x(:,8),'-', t(locs2),pks2, "bo",...
%    120, linspace(0.8,1,20), 'r.', ...
%    144, linspace(0.8,1,20), 'r.', ...
%    168, linspace(0.8,1,20), 'r.', ...
%    192, linspace(0.8,1,20), 'r.')
%text(ax2,t(locs2)+.01,pks2+.01,num2str(t(locs2)))
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([100 200])
%ylim([0.82,0.88])
%legend("FRQn")

ax1=nexttile;
plot(ax1,t,x(:,5),'-', t(locs1),pks1, "bo")
text(ax1,t(locs1)+.05,pks1+.05,num2str(t(locs1)))
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100 200])
legend("frq")

ax2=nexttile;
plot(ax2,t,x(:,8),'-', t(locs2),pks2, "bo")
text(ax2,t(locs2)+.01,pks2+.001,num2str(t(locs2)))
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100 200])
legend("FRQn")