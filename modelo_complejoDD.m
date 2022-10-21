function y = modelo_complejoDD(~,x)
%Parametros    
    k = [1.19, 1.2, 90, 1.6, 0.03, 0.226, 2.4, 2, 0.472, 0.3, 0.001, 0.6, 0.001, 5, 1, 6, 5.4, 0.15, 2, 0.05, 800, 0.68, 0.3];
    n1 = 2;
    n2 = 2;
    K = [0.03, 0.3, 0.05, 5, 2, 0.18, 0.02];
    kd = [2.4, 2.5, 0.135, 0.085, 0.05, 0.05, 0.02, 0.23, 0.69, 0.34, 0.34, 0.1, 6.2, 0.24, 0.24];
    
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
    L = 0;

%Ecuaciones
    function f = f_M(P_2,P_1,K_m)
        f = (P_1 + P_2 + K_m -  sqrt((P_1 + P_2 + K_m)^2 - 4*P_1*P_2))/2;
        %f = P_1*((P_2^n)/(P_2^n + K_m^n));
    end

    dwc1 = k(1) + k(2)*WCCn + k(3)*laWCC - kd(1)*wc1;
    dwc2 = k(4)/(1 + K(1)*WCCn) + k(5)*FFCn - kd(2)*wc2;
    dWC1 = k(6)*wc1 + k(7)*f_M(wc1,FFCp,K(2)) - k(9)*WC1*WC2 - kd(3)*WC1;
    dWC2 = k(8)*wc2 - kd(4)*WC2 - k(9)*WC1*WC2;
    dWCCc = k(9)*WC1*WC2 - kd(5)*WCCc -  k(10)*WCCc;
    dWCCn = k(10)*WCCc + k(11)*laWCC - kd(6)*WCCn - L*WCCn - k(13)*WCCn - k(12)*WCCn*(FFCn^n1/(K(3)^n1 + FFCn^n1));
    dlaWCC = L*WCCn + k(13)*WCCn - kd(7)*laWCC - k(11)*laWCC - k(14)*laWCC*(VVDn^n1/(K(4)^n1 + VVDn^n1));
    dfrq =  (k(15)*(K(6)*WCCn)^n2 + k(16)*(K(5)*laWCC)^n2)/((K(5)*K(6))^n2 + (K(6)*WCCn)^n2 + (K(5)*laWCC)^n2) - kd(8)*frq;
    dFRQ =  k(17)*frq - kd(9)*FRQ - k(18)*FRQ;
    dFFC = k(18)*FRQ - kd(10)*FFC - k(19)*FFC;
    dFFCp = k(19)*FFC - kd(11)*FFCp - k(20)*FFCp;
    dFFCn = k(20)*FFCp - kd(12)*FFCn;
    dvvd = (k(21)*laWCC^n2)/(K(7)^n2 + laWCC^n2) - kd(13)*vvd;
    dVVDc = k(22)*vvd - kd(14)*VVDc -  k(23)*VVDc;
    dVVDn = k(23)*VVDc - kd(15)*VVDn;

    y= [dwc1; dwc2; dWC1; dWC2; dWCCc; dWCCn; dlaWCC; dfrq; dFRQ; dFFC; dFFCp; dFFCn; dvvd; dVVDc; dVVDn];

end
