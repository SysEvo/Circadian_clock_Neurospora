% Simulate stochastic dynamics:
t = 0;              % Current (initial) time point
x = [12700,15819,496,332515,8851,1338,1,99,491,25,103,248,0,12,1607];       % Current (initial) state (size # of species)
%x=ones(1,15);
tMAX = 264;         % Maximum simulation time
prompt="Name file to save:";
file=input(prompt, 's');


[T,X] = Gillespie(x,t,tMAX,file);

%Function of the Gillespie algorithm.
function [T,X] = Gillespie(x,t,tMAX ,file)
    T = zeros(3*10^8,1);  % Simulated time points
    X = zeros(3*10^8,15);  % Simulated states (size # of species)
    i = 0;
    ranum1=rand(2*10^8,1);
    ranum2=rand(2*10^8,1);    
    while(t<tMAX)
        i = i + 1;
        X(i,:) = x;
        T(i)   = t;
        i = i + 1;
        X(i,:) = x;
	[t,x] = modelo_st(t , x , 120,ranum1(i/2),ranum2(i/2));
        T(i)   = t;
        m=mod(t,24);
        if m>23.999
            %save the data
            save("Data_st1.mat", "X", "T","-v7.3")
        end
    end
    z=find(T,1,"last");
    T=T(1:z);
    X=X(1:z,:);
    save(file, "X", "T","-v7.3")
end
