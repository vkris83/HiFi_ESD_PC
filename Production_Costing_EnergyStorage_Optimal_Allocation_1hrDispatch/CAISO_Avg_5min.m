%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% Estimates regulation requirements for a wind penetration 
% Change data for new wind penetration and re-run!

clear all;
load CAISOdata

%% change here for different wind penetration
wp1 = 300; % Bus 17
wp2 = 400; % Bus 21
wp3 = 300; % Bus 22
pen = 100*(wp1+wp2+wp3)/(3405+(wp1+wp2+wp3)); % penetration percentage

LOAD_CAISO = LOAD_CAISO*3000/(max(LOAD_CAISO)); % max 24 bus syst load w/o load shedding ~ 3000MW
N = 1.2*max(WIND_CAISO)/100;% 100 MW base plant

%% Below changes with wind plant rating/penetration
NL = LOAD_CAISO-(WIND_CAISO*(wp1/100)/N)-(WIND_CAISO*(wp2/100)/N)-(WIND_CAISO*(wp3/100)/N);
CAISO_wind_proj = WIND_CAISO*(wp1+wp2+wp3)/(100*N);

Int = 5;

T5min = length(LOAD_CAISO)/Int; % No. of 5 min. intervals

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
 
%% Hourly regulation from 5-min Net load data

h = length(LOAD_CAISO)/60; % hours
Int5 = T5min/h; % 5 min intervals in 1 hour

for b = 1:h
    NLVhr(b,1) = max(NLVar5((b-1)*12+1:(b-1)*12+12,1)); % for slave_UC
    NLVhravg(b,1) = mean(NLVar5((b-1)*12+1:(b-1)*12+12,1)); % for slave_ED
end

figure
plot(3*NLVhr)
xlabel('hour');
ylabel('3*sigma in MW');
title('Hourly Regulation');

save Reg_req NLVhr NLVhravg NLVar5