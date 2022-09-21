clear; close all; clc;

x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
tM = 0:0.05:228;   
[t,x] = ode45(@modelo_complejo,tM,x0);

x01 = x(end,:);
tM1 = 228:0.05:480;
[t1,x1] = ode45(@modelo_complejoDD,tM1,x01);

% Plot ODE results
X = [x; x1];
y = [t; t1];

%[pks1,locs1] = findpeaks(X(:,8));   %To localize the maximus points.
%[pks2,locs2] = findpeaks(X(:,12));

plot(y,X(:,[8,12]),'-')
%text(y(locs1)+.01,pks1,num2str(y(locs1)),"FontSize",6)
%text(y(locs2)+.05,pks2,num2str(y(locs2)),"FontSize",5)
title("LD DD 2.0")
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100,480])
%ylim([2.2,2.4])
legend(["frq","FRQn"])

%Plot ODE results to compare molecules with different scale of concentration.
%ax1=nexttile;
%plot(ax1,t1,x1(:,8),'-', t1(locs1),pks1, "bo")
%text(ax1,t1(locs1)+.05,pks1+.005,num2str(t1(locs1)))
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([300 480])
%legend("frq")

%ax2=nexttile;
%plot(ax2,t1,x1(:,12),'-', t1(locs2),pks2, "bo")
%text(ax2,t1(locs2)+.01,pks2+.001,num2str(t1(locs2)))
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([300 480])
%legend("FRQn")
