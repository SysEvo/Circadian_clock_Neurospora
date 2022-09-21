function y = modelo_complejo(tm,x)
%Parametros    
    k = [1.19, 1.2, 90, 1.6, 0.03, 0.226, 2.4, 1, 0.472, 0.3, 0.001, 50, 5, 0.001, 7.3, 320, 5.4, 0.15, 2, 0.05, 800, 0.68, 0.3];
    n = 2;
    K = [0.3, 0.5, 5, 1.8, 0.18, 0.02];
    kd = [2.4, 2.5, 0.135, 0.085, 0.05, 0.05, 0.02, 0.23, 0.69, 0.34, 0.34, 0.1, 6.2, 0.24, 0.24];
    
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
    FFC = x(10);
    FFCp = x(11);
    FFCn = x(12);
    vvd = x(13);
    VVDc = x(14);
    VVDn = x(15);
%Tiempo y periodos LD
    m=mod(tm,24);
    if m>=12
        L = l(1);
    else
        L= l(2);
    end
    
%Ecuaciones
    dwc1 = k(1) + k(2)*WCCn + k(3)*laWCC - kd(1)*wc1;
    dwc2 = k(4)/(1 + WCCn) + k(5)*FFCn - kd(2)*wc2;
    dWC1 = k(6)*wc1 + k(7)*(FFCp^n/(K(1)^n + FFCp^n))*wc1 - kd(3)*WC1 - k(9)*WC1*WC2;
    dWC2 = k(8)*wc2 - kd(4)*WC2 - k(9)*WC1*WC2;
    dWCCc = k(9)*WC1*WC2 - kd(5)*WCCc - k(10)*WCCc;
    dWCCn = k(10)*WCCc + k(11)*laWCC - kd(6)*WCCn - k(12)*(FFCn^n/(K(2)^n + FFCn^n))*WCCn - L*WCCn;
    dlaWCC = L*WCCn - kd(7)*laWCC - k(13)*(VVDn^n/(K(3)^n + VVDn^n))*laWCC - k(14)*laWCC;
    dfrq = k(15)*WCCn^n/(K(4)^n + WCCn^n) + k(16)*laWCC^n/(K(5)^n + laWCC^n) - kd(8)*frq;
    dFRQ = k(17)*frq - kd(9)*FRQ - k(18)*FRQ;
    dFFC = k(18)*FRQ - kd(10)*FFC - k(19)*FFC;
    dFFCp = k(19)*FFC - kd(11)*FFCp - k(20)*FFCp;
    dFFCn = k(20)*FFCp - kd(12)*FFCn;
    dvvd = k(21)*laWCC^n/(K(6)^n + laWCC^n) - kd(13)*vvd;
    dVVDc = k(22)*vvd - kd(14)*VVD - k(23)*VVD;
    dVVDn = k(23)*VVD - kd(15)*VVDn;



    y= [dwc1; dwc2; dWC1; dWC2; dWCCc; dWCCn; dlaWCC; dfrq; dFRQ; dFFC; dFFCp; dFFCn; dvvd; dVVD; dVVDn];

end
