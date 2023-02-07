%Function that contains the propensities to stochastic simulations of the circadian clock of N. crassa.
%t variable input is the actual time of the simulation, x is a vector with the concentration of the molecules at time t,
%tDD is the time that you want to simulate at LD cycles, rnum1 and rnum2 are the random numbers to the next time and reaction event.
function [t,x] = modelo_st(t,x,tDD,rnum1,rnum2)

    %Parameters
    k = [29750, 1.2, 90, 40000, 0.03, 0.226, 2.4, 2, 1.888e-05, 0.3, 0.001, 50, 0.001, 20, 182500, 8000000, 5.4, 0.15, 2, 0.05, 20000000, 0.68, 0.3];
    n1 = 2;
    n2 = 3;
    K = [1.2e-6, 7500, 1250, 125000, 50000, 4500, 500];
    kd = [2.4, 2.5, 0.135, 0.085, 0.05, 0.05, 6, 1, 0.69, 0.34, 0.34, 0.1, 6.2, 0.24, 0.24];
    
    l = [0, 0.2];
    %Periods of LD cycles at 12:12, and after tDD it changes to DD cycles.
    m=mod(t,24);
    if m>=12 || t>tDD
        L = l(1);
    else
        L= l(2);
    end

    %Reactions
    function f = f_M(P_2,P_1,K_m)
        f = (P_1 + P_2 + K_m -  sqrt((P_1 + P_2 + K_m)^2 - 4*P_1*P_2))/2;
    end

    lambda = zeros(39,1);
    lambda(1) = k(1); 
    lambda(2) = k(2)*x(6); 
    lambda(3) = k(3)*x(7); 
    lambda(4) = k(4)/(1 + K(1)*x(6)); 
    lambda(5) = k(5)*x(12); 
    lambda(6) = k(6)*x(1); 
    lambda(7) = k(7)*f_M(x(1),x(11),K(2)); 
    lambda(8) = k(8)*x(2); 
    lambda(9) = k(9)*x(3)*x(4); 
    lambda(10) = k(10)*x(5); 
    lambda(11) = k(11)*x(7); 
    lambda(12) = k(12)*x(6)*(x(12)^n1/(x(12)^n1 + K(3)^n1)); 
    lambda(13) = k(13)*x(6);
    lambda(14) = k(14)*x(7)*(x(15)^n1/(x(15)^n1 + K(4)^n1)); 
    lambda(15) = (k(15)*(K(6)*x(6))^n2)/((K(5)*K(6))^n2 + (K(6)*x(6))^n2 + (K(5)*x(7))^n2); 
    lambda(16) = (k(16)*(K(5)*x(7))^n2)/((K(5)*K(6))^n2 + (K(6)*x(6))^n2 + (K(5)*x(7))^n2); 
    lambda(17) = k(17)*x(8);
    lambda(18) = k(18)*x(9); 
    lambda(19) = k(19)*x(10); 
    lambda(20) = k(20)*x(11); 
    lambda(21) = k(21)*(x(7)^n2/(K(7)^n2 + x(7)^n2)); 
    lambda(22) = k(22)*x(13); 
    lambda(23) = k(23)*x(14); 
    lambda(24) = L*x(6); 
    lambda(25) = kd(1)*x(1); 
    lambda(26) = kd(2)*x(2); 
    lambda(27) = kd(3)*x(3); 
    lambda(28) = kd(4)*x(4); 
    lambda(29) = kd(5)*x(5); 
    lambda(30) = kd(6)*x(6); 
    lambda(31) = kd(7)*x(7); 
    lambda(32) = kd(8)*x(8); 
    lambda(33) = kd(9)*x(9); 
    lambda(34) = kd(10)*x(10); 
    lambda(35) = kd(11)*x(11); 
    lambda(36) = kd(12)*x(12); 
    lambda(37) = kd(13)*x(13); 
    lambda(38) = kd(14)*x(14); 
    lambda(39) = kd(15)*x(15);

    %Stequiometric matrix
    V = [[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,-1,-1,1,0,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,1,-1,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0];
         [0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1];
         [0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0];
         [-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0];
         [0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0];
         [0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0];
         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]];

    %Aleatory time.
    tR = -log(1-rnum1)/sum(lambda);
    t = t + tR;

    %Reaction that occur at the aleatory event.
    rR = sum(rnum2>(cumsum(lambda)/sum(lambda)))+1;
    x = x + V(rR,:);
end
