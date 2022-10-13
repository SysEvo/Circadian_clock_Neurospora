% Simulate stochastic dynamics:
t = 0;              % Current (initial) time point
x = ones(1,15)*10;       % Current (initial) state (size # of species)

tMAX = 100;         % Maximum simulation time

y = parfeval(@Gillespie, 2, x, t, tMAX);
[T,X] = fetchOutputs(y);
% Plot stochastic dynamics

% Plot stochastic dynamics

plot(T,X(:,[8 12]),'-','MarkerSize',3)
xlabel("Tiempo")
ylabel("Concentracion")
title("LD n=3")
xlim([100, 200])
legend(["frq","FFCn"])

%ax1=nexttile;
%plot(ax1,T,X(:,8),'-','MarkerSize',3)
%title("LD 12:12")
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([100, 200])
%legend("frq")

%ax2=nexttile;
%plot(ax2,T,X(:,12),'-','MarkerSize',3)
%xlabel("Tiempo")
%ylabel("Concentracion")
%xlim([100, 200])
%legend("FFCn")

%Save
save("Data_st.txt", "X", "T", "-ascii")

%Function of the Gillespie algorithm.
function [T,X] = Gillespie(x,t,tMAX )
    T = zeros(1000,1);  % Simulated time points
    X = zeros(1000,15);  % Simulated states (size # of species)
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
end
