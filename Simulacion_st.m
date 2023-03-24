% Simulate stochastic dynamics:
t = 0;              % Current (initial) time point
x = [12700,15819,496,332515,8851,1338,1,99,491,25,103,248,0,12,1607];       % Current (initial) state (size # of species)

tMAX = 100;         % Maximum simulation time
prompt="Name of the file to save:";	%File name
file=input(prompt, 's');

%Simulation with Gillespie algorithm.
[T,X] = Gillespie(x,t,tMAX,file);

%Function of the Gillespie algorithm.
function [T,X] = Gillespie(x,t,tMAX ,file)
    T = zeros(3*10^6,1);  % Simulated time points
    X = zeros(3*10^6,15);  % Simulated states (size # of species)
    i = 0;
    ranum1=rand(2*10^6,1);	%Random numbers to the aleatory time.
    ranum2=rand(2*10^6,39);    %Random numbers to choose wich reaction occur at an event.
    %Random walk
    while(t<tMAX)
        i = i + 1;
        X(i,:) = x;
        T(i)   = t;
	[t,x] = modelo_st(t , x , 120,ranum1(i),ranum2(i,:));
        m=mod(t,24);
        if m>23.999
            %save backups each ~24 hrs.
            save("Data_st.mat", "X", "T","-v7.3");
        end
    end
    %Save the data
    z=find(T,1,"last");
    T=T(1:z);
    X=X(1:z,:);
    save(file, "X", "T","-v7.3");
end
