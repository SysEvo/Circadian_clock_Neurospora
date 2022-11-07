% Simulate stochastic dynamics:
t = 0;              % Current (initial) time point
x = [12000,1,19664,4,5,25,12,11,95,8,29,21,2030,3186,6648];       % Current (initial) state (size # of species)

tMAX = 192;         % Maximum simulation time


[T,X] = Gillespie(x,t,tMAX);
% Plot stochastic dynamics

% Plot stochastic dynamics

%plot(T,X(:,[8 12]),'-','MarkerSize',3)
%xlabel("Tiempo")
%ylabel("Concentracion")
%title("LD n=3")
%xlim([100, 200])
%legend(["frq","FFCn"])

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

%Function of the Gillespie algorithm.
function [T,X] = Gillespie(x,t,tMAX )
    T = zeros(5*10^7,1);  % Simulated time points
    X = zeros(5*10^7,15);  % Simulated states (size # of species)
    i = 0;
    while(t<tMAX)
        i = i + 1;
        X(i,:) = x;
        T(i)   = t;
        i = i + 1;
        X(i,:) = x;
        [t,x] = modelo_st(t , x );
        T(i)   = t;
        m=mod(t,24);
        if m>23.999
            %save the data
            save("Data_st1.mat", "X", "T","-v7.3")
        end
    end
end