% Simulate stochastic dynamics:
t = 0;              % Current (initial) time point

x = [1200,2,1966,8,10,50,24,22,95,32,29,42,2030,3186,6648];       % Current (initial) state (size # of species)
%x=ones(1,15);
tMAX = 320;         % Maximum simulation time
prompt="Name file of save:";
file=input(prompt, 's');


[T,X] = Gillespie(x,t,tMAX,file);
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
function [T,X] = Gillespie(x,t,tMAX ,file)
    T = zeros(6*10^7,1);  % Simulated time points
    X = zeros(6*10^7,15);  % Simulated states (size # of species)
    i = 0;
    ranum1=rand(5*10^7,1);
    ranum2=rand(5*10^7,1);    
    while(t<tMAX)
        i = i + 1;
        X(i,:) = x;
        T(i)   = t;
        i = i + 1;
        X(i,:) = x;
	[t,x] = modelo_st(t , x , 144,ranum1(i/2),ranum2(i/2))
        T(i)   = t;
        m=mod(t,12);
        if m>11.999
            %save the data
            save("Data_st1.mat", "X", "T","-v7.3")
        end
    end
end
