%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% Estimates regulation requirements for a wind penetration 
% Change data for new wind penetration and re-run!

clear all;
load CAISOdata
%% 1 year sim
% load YData1min
% LOAD_CAISO = load_1min;
% WIND_CAISO = wind_1min;
% %%

% load loadhourly
%LWcor = corr(LOAD_CAISO,WIND_CAISO);
% LWcov = cov(LOAD_CAISO,WIND_CAISO);
% % % wp1 = 1.2*max(WIND_CAISO);
wp1 = 1530;
wp2 = 2040;
wp3 = 1530;

pen = 100*(wp1+wp2+wp3)/(3405+(wp1+wp2+wp3));

LOAD_CAISO = LOAD_CAISO*3000/(max(LOAD_CAISO)); % max 24 bus syst load w/o load shedding ~ 3000MW
N = 1.2*max(WIND_CAISO)/100;% 100 MW base plant
%% Below changes with wind plant rating/penetration
NL = LOAD_CAISO-(WIND_CAISO*(wp1/100)/N)-(WIND_CAISO*(wp2/100)/N)-(WIND_CAISO*(wp3/100)/N);
CAISO_wind_proj = WIND_CAISO*(wp1+wp2+wp3)/(100*N);
%Int = 60;
Int = 5;
%T60min = length(LOAD_CAISO)/60; % No. of 1 hour intervals
T5min = length(LOAD_CAISO)/Int; % No. of 5 min. intervals
for a = 1:T5min
    Load_hour(:,a) = LOAD_CAISO(((a-1)*Int)+1:a*Int);% data grouped for each hour for total 744 hours (Jan)
    Mean_hr(a,1) = mean(Load_hour(:,a));% 744 hours mean obtained
    Var_hr(a,1) = std(Load_hour(:,a));
    Wind_hour(:,a) = CAISO_wind_proj(((a-1)*Int)+1:a*Int);% data grouped for each hour for total 744 hours (Jan)
    WindMean_hr(a,1) = mean(Wind_hour(:,a));% 744 hours mean obtained
    WindVar_hr(a,1) = std(Wind_hour(:,a))/sqrt(N/(sum(wp1+wp2+wp3)/100));  
    NetloadVar_hr(a,1) = sqrt(WindVar_hr(a,1)^2+Var_hr(a,1)^2);
end

for a = 1:T5min
    %% net load actual
    NL5(:,a) = NL(((a-1)*Int)+1:a*Int);% data grouped for each hour for total 744 hours (Jan)
    NLMean5(a,1) = mean(NL5(:,a));% 744 hours mean obtained
    NLVar5(a,1) = std(NL5(:,a)); 
end

plot(LOAD_CAISO)
hold on
plot(CAISO_wind_proj,'r')
hold on
plot(NL,'g.')
xlabel('Minute');
ylabel('MW');
title('NET Load');
legend('Load','Wind','Net Load');
% 
% figure
% % plot(LOAD_CAISO(1:60,1))
% % hold on
% % plot(CAISO_wind_proj(1:60,1),'r')
% % hold on
% plot(NL(1:60,1),'g')
% % hold on
% % plot(NL(301:360,1),'r')% hour 7
% % hold on
% % plot(NL(1201:1260,1))% hour 20
% xlabel('minute ');
% ylabel('MW');
% title('Hour 1 NET Load CAISO');
% % legend('Hour 1','Hour 7','Hour 20');
% 
% figure
% %  plot(Var_hr)
% %  hold on
% %  plot(NetloadVar_hr,'r')
% %   hold on
% plot(3*NLVar5,'k')
% xlabel('Minute');
% ylabel('Estimated Regulation Requirement in MW');
% title('5-minute Estimated Regulation Requirement');
% % plot(NetloadVar_hr,'g')
% 
 
%% Hourly regulation from 5-min Net load data

h = length(LOAD_CAISO)/60; % hours
Int5 = T5min/h; % 5 min intervals in 1 hour

for b = 1:h
    %% net load actual
    NL_hr(:,b) = NLMean5(((b-1)*Int5)+1:b*Int5);% data grouped for each hour for total 744 hours (Jan)
    NLMean_hr(b,1) = mean(NL_hr(:,b));% 744 hours mean obtained
    NLVar_hr(b,1) = std(NL_hr(:,b)); 
    NLVhr(b,1) = max(NLVar5((b-1)*12+1:(b-1)*12+12,1));
end

figure
 plot(3*NLVhr)
%  hold on
% plot(WindVar_hr,'r')
% %  hold on
% plot(3*NLVar_hr,'r')
 xlabel('hour');
ylabel('3*sigma in MW');
title('Hourly Regulation');

save Reg_req5minED NLVhr NLVar5