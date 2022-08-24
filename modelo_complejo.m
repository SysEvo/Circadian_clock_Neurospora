function y = modelo_complejo(tm,x)
%Parametros
    k = [1.19, 90, 1.6, 0.03, 0.226, 2.4, 1, 0.472, 0.3, 0.001, 1, 5, 0.001, 7.3, 320, 5.4, 0.15, 2, 0.05, 800, 0.68, 0.3];
    n = 2;
    K = [1, 0.03, 0.3, 0.5, 5, 1.8, 0.18, 0.02];
    kd = [2.4, 2.5, 0.135, 0.085, 0.05, 0.05, 6, 0.23, 0.69, 0.34, 0.34, 0.1, 6.2, 0.24, 0.24];
    l = [0, 0.2];
%mRNAs y Proteinas
    wc1 = x(1);
    wc2 = x(2);
    WC1 = x(3);
    WC2 = x(4);
    WCCc = x(5);
    WCCn = x(6);
    laWCC = x(7);
    frq = x(8);
    FRQ = x(9);
    FRQ_FRH = x(10);
    FRQp = x(11);
    FRQn = x(12);
    vvd = x(13);
    VVD = x(14);
    VVDn = x(15);
%Tiempo y periodos LD
    m=mod(tm,4);
    if m>=2
        L = l(1);
    else
        L= l(2);
    end
    
%Ecuaciones
    dwc1 = k(1) + k(2)*laWCC^n/(K(1)^n+laWCC^n) - kd(1)*wc1;
    dwc2 = k(3)/(K(2) + WCCn) + k(4)*FRQn - kd(2)*wc2;
    dWC1 = k(5)*wc1 + k(6)*(FRQp^n/(K(3)^n + FRQp^n))*wc1 - kd(3)*WC1 - k(7)*WC1*WC2;
    dWC2 = k(8)*wc2 - kd(4)*WC2 - k(7)*WC1*WC2;
    dWCCc = k(7)*WC1*WC2 - kd(5)*WCCc - k(9)*WCCc;
    dWCCn = k(9)*WCCc + k(10)*laWCC - kd(6)*WCCn - k(11)*(FRQn^n/(K(4)^n + FRQn^n))*WCCn - L*WCCn;
    dlaWCC = L*WCCn - kd(7)*laWCC - k(12)*(VVDn^n/(K(5)^n + VVDn^n))*laWCC - k(13)*laWCC;
    dfrq = k(14)*WCCn^n/(K(6)^n + WCCn^n) + k(15)*laWCC^n/(K(7)^n + laWCC^n) - kd(8)*frq;
    dFRQ = k(16)*frq - kd(9)*FRQ - k(17)*FRQ;
    dFRQ_FRH = k(17)*FRQ - kd(10)*FRQ_FRH - k(18)*FRQ_FRH;
    dFRQp = k(18)*FRQ_FRH - kd(11)*FRQp - k(19)*FRQp;
    dFRQn = k(19)*FRQp - kd(12)*FRQn;
    dvvd = k(20)*laWCC^n/(K(8)^n + laWCC^n) - kd(13)*vvd;
    dVVD = k(21)*vvd - kd(14)*VVD - k(22)*VVD;
    dVVDn = k(22)*VVD - kd(15)*VVDn;



    y= [dwc1; dwc2; dWC1; dWC2; dWCCc; dWCCn; dlaWCC; dfrq; dFRQ; dFRQ_FRH; dFRQp; dFRQn; dvvd; dVVD; dVVDn];

end