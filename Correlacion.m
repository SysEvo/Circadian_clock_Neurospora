x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
tM = 0:0.05:450;   
[t,x] = ode45(@modelo_complejo,tM,x0);
t = t(1001:end);
x = x(1001:end,:);
l = 3560;
period1 = 24;

[pks1,locs1] = findpeaks(x(:,8),"MinPeakDistance",20);
period2 = diff(t(locs1([end-1, end])));
%LD period1
d_ind=ceil(period1/(t(2)-t(1)));
for indd=1:15
     me=mean(x(1:l,indd));
     acf=mean((x(1:(l-d_ind),indd)-me).*(x((1+d_ind):l,indd)-me));
     s=mean((x(1:l,indd)-me).*(x(1:l,indd)-me));
     
     if ~isempty(find(autocorr(x(1:l,indd),'NumLags',d_ind-1)<0, 1))
         if acf>0
             tau=-1/log((acf/s));
         else
             tau=0;
         end
     else
         tau=0;
     end
     
     
    eval(['tau',num2str(indd),'=tau;']);
    eval(['acf',num2str(indd),'=acf;']);
end

d_ind=ceil(3*period1/(t(2)-t(1)));
t_cal=t(1:d_ind);
%for indd=1:15
%    acf=autocorr(x(1:l,indd),'NumLags',d_ind-1);
%    acf=acf.*length(x(1:l,:))./(length(x(1:l,:))-(0:(d_ind-1)))';
%    warning('off');[acf_peak,loc]=findpeaks(acf,t_cal,'MinPeakHeight',0,'MinPeakDistance',period1/2);
%    eval(['acf',num2str(indd),'=acf;']);
%end

%DD period2
d_ind=ceil(period2/(t(end)-t(end-1)));
for indd=1:15
     me=mean(x(l:end,indd));
     acf=mean((x(l:(end-d_ind),indd)-me).*(x((l+d_ind):end,indd)-me));
     s=mean((x(l:end,indd)-me).*(x(l:end,indd)-me));
     if ~isempty(find(autocorr(x(l:end,indd),'NumLags',d_ind-1)<0, 1))
         if acf>0
             tau=-1/log((acf/s));
         else
             tau=0;
         end
     else
         tau=0;
     end
     
     
    eval(['tau_',num2str(indd),'=tau;']);
    eval(['acf_',num2str(indd),'=acf;']);
end

d_ind=ceil(3*period2/(t(end)-t(end-1)));
t_cal1=t(l:l+d_ind-1);
%for indd=1:15
%    acf=autocorr(x(l:end,indd),'NumLags',d_ind-1);
%    acf=acf.*length(x(l:end,:))./(length(x(l:end,:))-(0:(d_ind-1)))';
%    warning('off');[acf_peak2,loc2]=findpeaks(acf,t_cal1,'MinPeakHeight',0,'MinPeakDistance',period2/2);
%    eval(['acf_',num2str(indd),'=acf;']);
%end

figure;

subplot(2,2,1);  % autocorrelation function
plot(t_cal,[c_t(period1,tau8,t_cal), c_t(period1, tau12, t_cal), ...
    exp(-t_cal./(tau8*period1)), exp(-t_cal./(tau12*period1))],'-','linewidth',0.5);
title("Correlation")
legend('frq','FFCn');


subplot(2,2,2);plot(t(1:l),x(1:l,[8,12]));
title('Deterministic model LD');legend(["frq","FFCn"]);xlim([50, 228])


subplot(2,2,3);  % autocorrelation function
plot(t_cal1,[c_t(period2,tau_8,t_cal1), c_t(period2, tau_12, t_cal1), ...
    exp(-t_cal1./(tau_8*period2)), exp(-t_cal1./(tau_12*period2))],'-','linewidth',0.5);
title("Correlation")
legend('frq','FFCn');


subplot(2,2,4);plot(t(l:end),x(l:end,[8,12]));
title('Deterministic model DD');xlim([228,450]);


function eccor = c_t(period,tau,tcal)
    eccor = exp(-tcal./(tau*period)).*cos(2*pi*tcal/period);
end
