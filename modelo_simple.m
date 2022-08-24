function y = modelo_simple(tm,x)
%Parametros
    %k = [0.16, 0.5, 2.4, 0.2, 0.001 , 1, 5, 0.001, 1, 6, 5.4, 2, 0.4, 800, 0.68, 0.3];
    %K = [1, 0.5, 5, 1.8, 0.18, 0.02];
    %Kd = [0.18, 0.05, 0.02, 0.02, 1, 0.69, 0.34, 0.6, 6.2, 0.24, 0.24];
    
    kr = [1.19, 90, 0.226, 0.472, 0.001, 0.6, 20, 0.001, 7.3, 320, 0.19, 0.1, 0.1, 800, 0.68, 0.3];
    Kr = [1, 0.475, 5, 0.1, 0.18, 0.02];
    Kdr = [2.4, 0.135, 0.05, 100, 2, 0.27, 0.34, 0.27, 6.2, 0.24, 0.24];

    %ka = [0.16, 0.5, 2.4, 40, 0.001 , 50, 5, 0.001, 5.4, 6, 5.4, 2, 0.05, 1800, 10, 0.3];
    %Ka = [1, 0.5, 5, 1.25, 0.18, 0.02];
    %Kda = [0.1, 0.05, 0.02, 0.02, 0.23, 0.69, 0.34, 0.1, 6.2, 0.24, 0.24];
    
    n = 2;
    l = [0, 0.2];
%mRNAs y Proteinas
    wc1 = x(1);
    WC1 = x(2);
    WCCn = x(3);
    laWCC = x(4);
    frq = x(5);
    FRQ = x(6);
    FRQp = x(7);
    FRQn = x(8);
    vvd = x(9);
    VVD = x(10);
    VVDn = x(11);
%Tiempo y periodos LD
    m=mod(tm,24);
    if m>=12
        L = l(1);
    else
        L= l(2);
    end
    
%Ecuaciones
    dwc1 = kr(1) + kr(2)*laWCC^n/(Kr(1)^n+laWCC^n) - Kdr(1)*wc1;
    dWC1 = kr(3)*wc1 - kr(4)*WC1 - Kdr(2)*WC1;
    dWCCn = kr(4)*WC1 + kr(5)*laWCC - Kdr(3)*WCCn - kr(6)*(FRQn^n/(Kr(2)^n + FRQn^n))*WCCn - L*WCCn;
    dlaWCC = L*WCCn - Kdr(4)*laWCC - kr(7)*(VVDn^n/(Kr(3)^n + VVDn^n))*laWCC - kr(8)*laWCC;
    dfrq = kr(9)*WCCn^n/(Kr(4)^n + WCCn^n) + kr(10)*laWCC^n/(Kr(5)^n + laWCC^n) - Kdr(5)*frq;
    dFRQ = kr(11)*frq - kr(12)*FRQ - Kdr(6)*FRQ;
    dFRQp = kr(12)*FRQ - kr(13)*FRQp - Kdr(7)*FRQp;
    dFRQn = kr(13)*FRQp - Kdr(8)*FRQn;
    dvvd = kr(14)*laWCC^n/(Kr(6)^n + laWCC^n) - Kdr(9)*vvd;
    dVVD = kr(15)*vvd - kr(16)*VVD - Kdr(10)*VVD;
    dVVDn = kr(16)*VVD - Kdr(11)*VVDn;


    y= [dwc1; dWC1; dWCCn; dlaWCC; dfrq; dFRQ; dFRQp; dFRQn; dvvd; dVVD; dVVDn];

end