% Trishna Das & Venkat Krishnan, Iowa State University
% used to plot the ancillary service results (Regulation, spinning and
% non-spinning reserves etc through both charging & discharging operations
subplot(4,1,1)
plot(CAES_Comspin(49:end),'r');
hold on
plot(CAES_spin(49:end),'DisplayName','CAES_Comdownreg(49:96,1)','YDataSource','CAES_Comdownreg(49:96,1)');figure(gcf)
hold on
plot(CAES_nonspin2(49:end),'g');
xlabel('hours')
ylabel('MW-hr')
legend('Charge SR','DisCharge SR','DisCharge NSR')
subplot(4,1,2)
plot(CAES_Comupreg(49:end),'r');
hold on
plot(CAES_Comdownreg(49:end),'DisplayName','CAES_Comdownreg(49:96,1)','YDataSource','CAES_Comdownreg(49:96,1)');figure(gcf)
xlabel('hours')
ylabel('MW-hr')
legend('Charge RU','Charge RD')
% title('Compressor Regulation')
subplot(4,1,3)
plot(CAES_upreg(49:end),'r');
hold on
plot(CAES_downreg(49:end),'DisplayName','CAES_downreg(49:96,1)','YDataSource','CAES_downreg(49:96,1)');figure(gcf)
xlabel('hours')
ylabel('MW-hr')
legend('DisCharge RU','DisCharge RD')
% title('Turbine Regulation')
% subplot(5,2,5)
% xlabel('hours')
% ylabel('MW')
% title('Compressor Down Regulation')
% subplot(5,2,6)
% 
% xlabel('hours')
% ylabel('MW')
% title('Turbine Down Regulation')
subplot(4,1,4)
plot(CAES_strlevel(49:end),'DisplayName','CAES_strlevel(49:96,1)','YDataSource','CAES_strlevel(49:96,1)');figure(gcf)
xlabel('hours')
ylabel('MW-hr')
legend('Storage Level')

figure
plot(1:1:dpi,CAES_charge(49:end),'r');
hold;
[AX,H1,H2]=plotyy(1:1:dpi,-CAES_dischar(49:end),1:1:dpi,LMP_21);
title('CAES output');
xlabel('Hours')
ylabel('MW')
legend('Charging','discharging','LMP');

% subplot(5,2,1)
% plot(CAES_Comspin(49:96,1),'DisplayName','CAES_Comspin(49:96,1)','YDataSource','CAES_Comspin(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Compressor Spinning Reserve')
% subplot(5,2,2)
% plot(CAES_spin(49:96,1),'DisplayName','CAES_spin(49:96,1)','YDataSource','CAES_spin(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Turbine Spinning Reserve')
% subplot(5,2,3)
% plot(CAES_Comupreg(49:96,1),'DisplayName','CAES_Comupreg(49:96,1)','YDataSource','CAES_Comupreg(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Compressor Up Regulation')
% subplot(5,2,4)
% plot(CAES_upreg(49:96,1),'DisplayName','CAES_upreg(49:96,1)','YDataSource','CAES_upreg(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Turbine Up Regulation')
% subplot(5,2,5)
% plot(CAES_Comdownreg(49:96,1),'DisplayName','CAES_Comdownreg(49:96,1)','YDataSource','CAES_Comdownreg(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Compressor Down Regulation')
% subplot(5,2,6)
% plot(CAES_downreg(49:96,1),'DisplayName','CAES_downreg(49:96,1)','YDataSource','CAES_downreg(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Turbine Down Regulation')
% subplot(5,2,7)
% plot(CAES_strlevel(49:96,1),'DisplayName','CAES_strlevel(49:96,1)','YDataSource','CAES_strlevel(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('CAES_strlevel')
% subplot(5,2,8)
% plot(CAES_nonspin2(49:96,1),'DisplayName','CAES_nonspin2(49:96,1)','YDataSource','CAES_nonspin2(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Turbine Non spin2')
% subplot(5,2,9)
% plot(CAES_charge(49:96,1),'DisplayName','CAES_charge(49:96,1)','YDataSource','CAES_charge(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Compressor Charge')
% subplot(5,2,10)
% plot(CAES_dischar(49:96,1),'DisplayName','CAES_dischar(49:96,1)','YDataSource','CAES_dischar(49:96,1)');figure(gcf)
% xlabel('hours')
% ylabel('MW')
% title('Turbine Discharge')