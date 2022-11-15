
prompt="Name file of stochastic data:";
file=input(prompt, 's');
load(file)

y = find(T>120,1);
period1 = 24;
x = zeros(length(T),15);
for i = 1:15
    x(:,i)=smooth(X(:,i));
end

[pks1,locs1] = findpeaks(x(y:end,8),"MinPeakDistance",y/6,"MinPeakHeight",220);
period2 = diff(T(locs1([end-1, end])));
%LD period1
d_ind=find(T>period1,1);
for indd=1:15
     me=mean(X(1:y,indd));
     acf=mean((X(1:(y-d_ind),indd)-me).*(X((1+d_ind):y,indd)-me));
     s=mean((X(1:y,indd)-me).*(X(1:y,indd)-me));
     
     if ~isempty(find(autocorr(X(1:y,indd),'NumLags',d_ind-1)<0, 1))
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

d_ind=find(T>3*period1,1);
t_cal=T(1:d_ind);
%for indd=1:15
%    acf=autocorr(x(1:l,indd),'NumLags',d_ind-1);
%    acf=acf.*length(x(1:l,:))./(length(x(1:l,:))-(0:(d_ind-1)))';
%    warning('off');[acf_peak,loc]=findpeaks(acf,t_cal,'MinPeakHeight',0,'MinPeakDistance',period1/2);
%    eval(['acf',num2str(indd),'=acf;']);
%end

%DD period2
d_ind=find(T(y:end)>period2,1);
for indd=1:15
     me=mean(x(y:end,indd));
     acf=mean((x(y:(end-d_ind),indd)-me).*(x((y+d_ind):end,indd)-me));
     s=mean((x(y:end,indd)-me).*(x(y:end,indd)-me));
     if ~isempty(find(autocorr(x(y:end,indd),'NumLags',d_ind-1)<0, 1))
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

d_ind=find(T(y:end)>3*period2,1);
t_cal1=T(y:y+d_ind-1);
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


subplot(2,2,2);plot(T(1:y),X(1:y,[8,12]));
title('Deterministic model LD');legend(["frq","FFCn"]);xlim([0, 120])


subplot(2,2,3);  % autocorrelation function
plot(t_cal1,[c_t(period2,tau_8,t_cal1), c_t(period2, tau_12, t_cal1), ...
    exp(-t_cal1./(tau_8*period2)), exp(-t_cal1./(tau_12*period2))],'-','linewidth',0.5);
title("Correlation")
legend('frq','FFCn');


subplot(2,2,4);plot(T(y:end),X(y:end,[8,12]));
title('Deterministic model DD');xlim([120,264]);


function eccor = c_t(period,tau,tcal)
    eccor = exp(-tcal./(tau*period)).*cos(2*pi*tcal/period);
end