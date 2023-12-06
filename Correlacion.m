x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; 
tM = 0:0.05:600;   
[t,x] = ode45(@modelo_complejo,tM,x0);

l = find(t==260,1);
y= find(t==100,1);
z= find(t==350,1);
period1 = 24;

[pks1,locs1] = findpeaks(x(:,8),"MinPeakDistance",20);
period2 = diff(t(locs1([end-1, end])));
%LD period1
d_ind=ceil(period1/(t(2)-t(1)));
for indd=1:15
     me=mean(x(y:l,indd));
     acf=mean((x(y:(l-d_ind),indd)-me).*(x((y+d_ind):l,indd)-me));
     s=mean((x(y:l,indd)-me).*(x(y:l,indd)-me));
     
     if ~isempty(find(autocorr(x(y:l,indd),'NumLags',d_ind-1)<0, 1))
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
t_cal=t(y:d_ind);
%for indd=1:15
%    acf=autocorr(x(1:l,indd),'NumLags',d_ind-1);
%    acf=acf.*length(x(1:l,:))./(length(x(1:l,:))-(0:(d_ind-1)))';
%    warning('off');[acf_peak,loc]=findpeaks(acf,t_cal,'MinPeakHeight',0,'MinPeakDistance',period1/2);
%    eval(['acf',num2str(indd),'=acf;']);
%end

%DD period2
d_ind=ceil(period2/(t(end)-t(end-1)));
for indd=1:15
     me=mean(x(z:end,indd));
     acf=mean((x(z:(end-d_ind),indd)-me).*(x((z+d_ind):end,indd)-me));
     s=mean((x(z:end,indd)-me).*(x(z:end,indd)-me));
     if ~isempty(find(autocorr(x(z:end,indd),'NumLags',d_ind-1)<0, 1))
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
t_cal1=t(z:z+d_ind-1);
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
title("Correlation");
legend('frq','FFCn','decay frq','decay FFCn');


subplot(2,2,2);plot(t(y:l),x(y:l,[8,12]));
title('Deterministic model LD');legend(["frq","FFCn"]);

subplot(2,2,3);  % autocorrelation function
plot(t_cal1,[c_t(period2,tau_8,t_cal1), c_t(period2, tau_12, t_cal1), ...
    exp(-t_cal1./(tau_8*period2)), exp(-t_cal1./(tau_12*period2))],'-','linewidth',0.5);
title("Correlation");
legend('frq','FFCn','decay frq','decay FFCn');


subplot(2,2,4);plot(t(z:end),x(z:end,[8,12]));
title('Deterministic model DD');


function eccor = c_t(period,tau,tcal)
    eccor = exp(-tcal./(tau*period)).*cos(2*pi*tcal/period);
end
