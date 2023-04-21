clear; close all; clc;

x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];   %initial values of concentration of the molecules.
tM = 0:0.05:800;                                      %time to the simulation
[t,x] = ode45(@modelo_complejo,tM,x0);        %Solve the ODE system, x contains the concentration of the molecules at t time.


% Plot ODE results
%[pks1,locs1] = findpeaks(x(:,8));   %To localize the maximus points.
%[pks2,locs2] = findpeaks(x(:,12));

plot(t,x(:,6),'-')
%text(t(locs1)+.01,pks1,num2str(t(locs1)),"FontSize",6)
%text(t(locs2)+.05,pks2,num2str(t(locs2)),"FontSize",5)
title("LD DD 3.0")
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100,800])
%ylim([2.2,2.4])
%legend(["laWCC","frq","FRQn"])

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
