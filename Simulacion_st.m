clear; clc;

% Simulate stochastic dynamics:
T = zeros(1000,1);  % Simulated time points
X = zeros(1000,15);  % Simulated states (size # of species)

t = 0;              % Current (initial) time point
x = ones(1,15)*10;       % Current (initial) state (size # of species)

tMAX = 200;         % Maximum simulation time

i = 0;

while(t<tMAX)
    i = i + 1;
    X(i,:) = x;
    T(i)   = t;
    i = i + 1;
    X(i,:) = x;
    [t,x] = modelo_st(t , x );
    T(i)   = t;
end

clear i t x

% Plot stochastic dynamics

ax1=nexttile;
plot(ax1,T,X(:,8),'-','MarkerSize',3)
title("LD 12:12")
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100, 200])
legend("frq")

ax2=nexttile;
plot(ax2,T,X(:,12),'-','MarkerSize',3)
xlabel("Tiempo")
ylabel("Concentracion")
xlim([100, 200])
legend("FFCn")
