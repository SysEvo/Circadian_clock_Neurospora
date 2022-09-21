function [t,x] = modelo_st(t,x)

    %Parametros
    k = [1.19, 1.2, 0.5, 1.6, 0.03, 0.226, 2.4, 1, 0.472, 0.3, 0.001, 1, 5, 0.001, 7.3, 6, 5.4, 0.15, 2, 0.05, 800, 0.68, 0.3];
    n = 2;
    K = [0.3, 0.5, 5, 1.8, 0.18, 0.02];
    kd = [2.4, 2.5, 0.135, 0.085, 0.05, 0.05, 0.02, 0.23, 0.69, 0.34, 0.34, 0.1, 6.2, 0.24, 0.24];
    
    l = [0, 0.2];
    
    %Periodos LD
    m=mod(t,24);
    if m>=12
        L = l(1);
    else
        L= l(2);
    end

    %Reacciones
    lambda = zeros(39,1);
    lambda(1) = k(1); 
    lambda(2) = k(2)*x(6); 
    lambda(3) = k(3)*x(7); 
    lambda(4) = k(4)/(1 + x(6)); 
    lambda(5) = k(5)*x(12); 
    lambda(6) = k(6)*x(1); 
    lambda(7) = k(7)*x(1) *(x(11)^n/(K(1)^n + x(11)^n)); 
    lambda(8) = k(8)*x(2); 
    lambda(9) = k(9)*x(3)*x(4); 
    lambda(10) = k(10)*x(5); 
    lambda(11) = k(11)*x(7); 
    lambda(12) = k(12)*x(6)*(x(12)^n/(K(2)^n + x(12)^n)); 
    lambda(13) = k(13)*x(6);
    lambda(14) = k(14)*x(7)*(x(15)^n/(K(3)^n + x(15)^n)); 
    lambda(15) = k(15)*(x(6)^n/(K(4)^n + x(6)^n)); 
    lambda(16) = k(16)*(x(7)^n/(K(5)^n + x(7)^n)); 
    lambda(17) = k(17)*x(8);
    lambda(18) = k(18)*x(9); 
    lambda(19) = k(19)*x(10); 
    lambda(20) = k(20)*x(11); 
    lambda(21) = k(21)*(x(7)^n/(K(6)^n + x(7)^n)); 
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

    %Matriz estequiométrica
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

    %Tiempo
    tR = -log(1-rand())/sum(lambda);
    t = t + tR;

    %Reacción a ocurrir.
    rR = sum(rand()>(cumsum(lambda)/sum(lambda)))+1;
    x = x + V(rR,:);
end