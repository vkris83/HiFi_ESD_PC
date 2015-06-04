%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% Main program that initiates monte carlo simulation (changes in random
% gen. outages, prices...), and calls for programs Slave_UC.m (SCUC) and Slave_ED.m (SCED), and gets the output from SCED for plotting purposes 

clc
clear;
tic;
n=70;  % Number of monte carlo simulation
tolerance=0.001;
count=0;
wc=0;
lead_time=1;
carbon_tax=0; % $/ton of carbob dioxide emission
h = 24;
d = 2;
a = 2;
b = 24;
E2P = 4; % storage energy to power ratio (usually 4)

dpi = h*d; %data per item for simulation period with 10 min intervals
STOR1_n = 2;
STOR2_n = 21;% storage location; can be changed
load UpDn % minimum up and down time for generators

Slave_UC; % run the SCUC slave problem

%%%%%%%%%%%%%%%% update the wind generation information%%%%%%%%%%%%%%%%%%%%
% the upper bound of wind should be determined by wind speed, so it needs to be updated
load Windfc % wind hourly data
wind_gen = [wind_gen;wind_gen;wind_gen]; % same wind shape assumed for 3 wind farms at bus 17, 21, and 22
wo=wind_gen(1:num_wind+num_wind_new);
upbd=x_U_UC(st_wind:st_wind+num_wind+num_wind_new-1);
for i=1:length(upbd)
    if upbd(i)>10^6
        upbd(i)=10^6;
    end
    if upbd(i)==0
        upbd(i)=wc;
    end
end
wind_num_all=zeros(num_wind+num_wind_new,1);
wind_num_all(1:num_wind)=1;

Slave_ED; % run the SCED slave problem

% Monte Carlo simulation is executed (variables changed are fuel prices and gen./trans. outages)
solu_primal=0;
solu_dual=0;
B_rating=0;  % rating of branch if there is outage (hourly outage decided by random sampling)
G_rating=0; % rating of gen if there is outage
num_each=a*b;
dual_tmp_wind=zeros(num_wind_new,1); % the dual of wind generation node
tt_tmp=num_wind_new/b/a; % num. of new wind gen

% Deal with the oil prices
price1(1,1:n)= randn(1,n); % random number generation for Monte Carlo simulation

% Deal with the natural gas prices
price2(1,1:n)= randn(1,n);

x_U_UC(st_wind:st_wind+num_wind+num_wind_new-1)=wind_num_all.*upbd.*wind_gen(1:num_wind+num_wind_new); % update the upper bound of wind capacity at each hour
x_U_ED(st_wind:st_wind+num_wind+num_wind_new-1)=wind_num_all.*upbd.*wind_gen(1:num_wind+num_wind_new); % update the upper bound of wind capacity at each hour
% Remember: x (optimization solution) is determined, then everything evolves around x.

flag_mt=0;
Profit_mean=0;
for kk=1:length(x_U_UC)
    if x_U_UC(kk)==Inf
        x_U_UC(kk)=10^6;
    end
end
for kk=1:length(x_U_ED)
    if x_U_ED(kk)==Inf
        x_U_ED(kk)=10^6;
    end
end

for i=1:n  % Monte Carlo simulation for n times
    % First run the UC problem.
    fot_rate=zeros(num_gen,1);
    Branch_rate=zeros(num_br+num_new,1);
    for k1=st_gen:st_gen + num_gen -1 % Forced-outage rate
        fot_rate(k1-st_gen+1,1)=str2double(cell2mat(arcstmp{k1,29}));
    end
    
   Gen_udf=rand(num_gen,1); % For generators
   Branch_udf=rand(num_br+num_new,1); % for branches
    for j=1:num_gen
        if Gen_udf(j,1)>=fot_rate(num_gen,1)
            Gen_udf(j,1)=1;
        else
            Gen_udf(j,1)=G_rating; % Generator outage
        end
    end
    for k1=st_br:(st_br + num_br+num_new -1)
        Branch_rate(k1-st_br+1,1)=str2double(cell2mat(arcstmp{k1,29}));
    end
    for j=1:num_br+num_new
        if Branch_udf(j,1)>=Branch_rate(j,1)
            Branch_udf(j,1)=1;
        else
            Branch_udf(j,1)=B_rating; % Branch outage
        end
    end
    
    % the decomposition function changes x_L_UC and x_U_UC
    mt_c_UC=c_UC; % Prices
    mt_c_UC(1:num_each,1)=mt_c_UC(1:num_each,1)+price1(1,i)*mt_c_UC(1:num_each,1);
    mt_c_UC(num_each*2+1:num_each*3,1)=mt_c_UC(num_each*2+1:num_each*3,1)+price1(1,i)*mt_c_UC(num_each*2+1:num_each*3,1);
    mt_c_UC(num_each*4+1:num_each*5,1)=mt_c_UC(num_each*4+1:num_each*5,1)+price2(1,i)*mt_c_UC(num_each*4+1:num_each*5,1);
    mt_c_UC(num_each*5+1:num_each*6,1)=mt_c_UC(num_each*5+1:num_each*6,1)+price2(1,i)*mt_c_UC(num_each*5+1:num_each*6,1);
    mt_c_UC(num_each*6+1:num_each*7,1)=mt_c_UC(num_each*6+1:num_each*7,1)+price2(1,i)*mt_c_UC(num_each*6+1:num_each*7,1);
    mt_x_L_UC=x_L_UC;
    mt_x_L_UC(st_gen:st_gen+num_gen-1,1)=mt_x_L_UC(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_L_UC(st_br:st_br+num_br+num_new-1,1)=mt_x_L_UC(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_L_UC(obj_A:obj_A+num_gen-1,1)=mt_x_L_UC(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_L_UC(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_RU:obj_RU+num_gen-1,1)=mt_x_L_UC(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_RD:obj_RD+num_gen-1,1)=mt_x_L_UC(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC=x_U_UC;
    mt_x_U_UC(st_gen:st_gen+num_gen-1,1)=mt_x_U_UC(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_U_UC(st_br:st_br+num_br+num_new-1,1)=mt_x_U_UC(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_U_UC(obj_A:obj_A+num_gen-1,1)=mt_x_U_UC(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_U_UC(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC(obj_RU:obj_RU+num_gen-1,1)=mt_x_U_UC(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC(obj_RD:obj_RD+num_gen-1,1)=mt_x_U_UC(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC(obj_UC:obj_UC+num_gen-1,1)=mt_x_U_UC(obj_UC:obj_UC+num_gen-1,1).*Gen_udf(:,1); % Oct 10 2011
    mt_x_U_UC(obj_UCNsp:obj_UCNsp+num_gen-1,1)=mt_x_U_UC(obj_UCNsp:obj_UCNsp+num_gen-1,1).*Gen_udf(:,1); % Oct 10 2011
    mt_b_L_UC=b_L_UC; %% Oct 20 2011
    mt_b_L_UC(obj_row_S_pter:obj_row_S_pter+num_gen-1,1)=mt_b_L_UC(obj_row_S_pter:obj_row_S_pter+num_gen-1,1).*Gen_udf(:,1);
    
    [x_uc,slack_uc,v_uc,nouse_uc,objv_uc,nouse2_uc,nouse3_uc,inform_uc] = cplex(mt_c_UC,A_UC,mt_x_L_UC,mt_x_U_UC,mt_b_L_UC,b_U_UC,[],[],[],[],IntVars_UC);
    
    UC_status = zeros(length(obj_UC_num),1);
    UC_status = x_uc(obj_UC:obj_UC+obj_UC_num-1,1);
    
    Xst = x_uc(obj_X:obj_Y-1,1); % Start up status
    Yst = x_uc(obj_Y:obj_Y+(obj_Y-obj_X)-1,1); % Shut down status
    
    Xnst = x_uc(obj_Xn:obj_Yn-1,1); % Start up status NSP Oct 21
    Ynst = x_uc(obj_Yn:obj_Yn+(obj_Yn-obj_Xn)-1,1); % Shut down status NSP Oct 21
    
    Nsp_status = zeros(length(obj_UCNsp_num),1); %% Nsp status
    Nsp_status = x_uc(obj_UCNsp:obj_UCNsp+obj_UCNsp_num-1,1);  
    
    Unsp = UC_status + Nsp_status;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the decomposition function changes x_L_ED and x_U_ED
    mt_c_ED=c_ED; % Prices
    mt_c_ED(1:num_each,1)=mt_c_ED(1:num_each,1)+price1(1,i)*mt_c_ED(1:num_each,1);
    mt_c_ED(num_each*2+1:num_each*3,1)=mt_c_ED(num_each*2+1:num_each*3,1)+price1(1,i)*mt_c_ED(num_each*2+1:num_each*3,1);
    mt_c_ED(num_each*4+1:num_each*5,1)=mt_c_ED(num_each*4+1:num_each*5,1)+price2(1,i)*mt_c_ED(num_each*4+1:num_each*5,1);
    mt_c_ED(num_each*5+1:num_each*6,1)=mt_c_ED(num_each*5+1:num_each*6,1)+price2(1,i)*mt_c_ED(num_each*5+1:num_each*6,1);
    mt_c_ED(num_each*6+1:num_each*7,1)=mt_c_ED(num_each*6+1:num_each*7,1)+price2(1,i)*mt_c_ED(num_each*6+1:num_each*7,1);
    mt_x_L_ED=x_L_ED;
    mt_x_L_ED(st_gen:st_gen+num_gen-1,1)=mt_x_L_ED(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_L_ED(st_br:st_br+num_br+num_new-1,1)=mt_x_L_ED(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_L_ED(obj_A:obj_A+num_gen-1,1)=mt_x_L_ED(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_L_ED(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_RU:obj_RU+num_gen-1,1)=mt_x_L_ED(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_RD:obj_RD+num_gen-1,1)=mt_x_L_ED(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED=x_U_ED;
    mt_x_U_ED(st_gen:st_gen+num_gen-1,1)=mt_x_U_ED(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_U_ED(st_br:st_br+num_br+num_new-1,1)=mt_x_U_ED(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_U_ED(obj_A:obj_A+num_gen-1,1)=mt_x_U_ED(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_U_ED(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_RU:obj_RU+num_gen-1,1)=mt_x_U_ED(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_RD:obj_RD+num_gen-1,1)=mt_x_U_ED(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_b_L_ED=b_L_ED; %% Oct 20 2011
    mt_b_L_ED(obj_row_S_pter_ED:obj_row_S_pter_ED+num_gen-1,1)=mt_b_L_ED(obj_row_S_pter_ED:obj_row_S_pter_ED+num_gen-1,1).*Gen_udf(:,1);
    mt_b_U_ED=b_U_ED; %% Oct 20 2012
    
    % After updating the price and gen/tran availability, update the unit commitment decision.
    kng = 0;
    for kn=0:obj_UC_num-1  % To deal with the Inf*0 problem
        if (x_uc(obj_UC+kn))<0.5
            kng = kn;
            if kn>=num_gen %% wind gen. starts, skip until storage gen.
                kng = kn+num_wind;
            end               
                
            mt_x_U_ED(obj_A+kn)=0; 
            mt_x_U_ED(st_gen+kng)=0;
            mt_x_U_ED(obj_RU+kn)=0;
            mt_x_U_ED(obj_RD+kn)=0;

            mt_x_L_ED(obj_A+kn)=0; 
            mt_x_L_ED(st_gen+kng)=0;
            mt_x_L_ED(obj_RU+kn)=0;
            mt_x_L_ED(obj_RD+kn)=0;
            
            mt_b_L_ED(obj_row_S_pter_ED+kn)=0;           
        end
        if (x_uc(obj_UCNsp+kn))<0.5
            mt_x_U_ED(obj_Bn+kn)=0;                   
            mt_x_L_ED(obj_Bn+kn)=0;           
        end
    end

    %% June 27 2012 - Disjoint Hourly Charge/Discharge Operation in Bulk Storage
    %% (this section commented for short-term storage simulation)
    for kn=0:num_tur-1
        if (x_uc(obj_UC+obj_UC_num-num_tur+kn)+x_uc(obj_UCNsp+obj_UC_num-num_tur+kn))>0.5
            mt_x_U_ED(st_com+kn)=0;
            mt_x_L_ED(st_com+kn)=0;
            mt_b_U_ED(obj_row_StorComLDJ+kn)=0;            
        end
    end

   
    [x_ed,slack_ed,v_ed,nouse_ed,objv_ed,nouse2_ed,nouse3_ed,inform_ed] = cplex(mt_c_ED, A_ED, mt_x_L_ED, mt_x_U_ED, mt_b_L_ED, mt_b_U_ED,[],[],[],[],[],[],[],[],[],[],F_ED);
    if flag_mt==0
        flag_mt=1;
        xx_UC=zeros(length(x_uc),1);
        v_k_UC=zeros(length(v_uc),1);
        op_c_UC=0;
        xx_ED=zeros(length(x_ed),1);
        v_k_ED=zeros(length(v_ed),1);
        op_c_ED=0;
    end
    xx_UC=x_uc+xx_UC;
    op_c_UC=objv_uc+op_c_UC;
    xx_ED=x_ed+xx_ED;
    v_k_ED=v_ed+v_k_ED;
    op_c_ED=objv_ed+op_c_ED;
  
    STOR_charge=x_ed(st_com:st_com+num_com-1);
    STOR_dischar=x_ed(st_tur:st_tur+num_tur-1);
    STOR_strlevel=[x_ed(st_res:st_res+dpi-2);x_ed(obj_final_st);x_ed(st_res+dpi-1:st_res+num_res-1);x_ed(obj_final_st+1)];
    Wind_output=x_ed(st_wind:st_wind+num_wind+num_wind_new-1);
    gen=x_ed(st_res+num_res:st_res+num_res-1+dpi*14);
    spin=x_ed(st_spin:st_spin+num_res1-1);
    nonspin2=x_ed(obj_Bn:obj_Bn+num_res1-1); % oCT 09 2011
    Upreg=x_ed(obj_RU:obj_RU+num_reg-1);
    Downreg=x_ed(obj_RD:obj_RD+num_reg-1);
    
    LMP_2=v_ed(nst_tran+(STOR1_n - 1)*dpi:STOR1_n * dpi);
    LMP_21=v_ed(nst_tran+(STOR2_n - 1)*dpi:STOR2_n * dpi);
    MCP_SR2=v_ed(obj_row_2:obj_row_2 + dpi-1);
    MCP_NSR2=v_ed(obj_row_4:obj_row_4 + dpi-1);
    MCP_SR1=v_ed(obj_row_1:obj_row_1 + dpi-1);
    MCP_NSR1=v_ed(obj_row_3:obj_row_3 + dpi-1);
    MCP_R=v_ed(obj_row_ru:obj_row_ru + dpi-1);
    MCP_D=v_ed(obj_row_rd:obj_row_rd + dpi-1);

    MCP_ru = MCP_R+MCP_SR1+MCP_NSR1;
    MCP_rd = MCP_D;
    MCP_1 = MCP_SR1+MCP_NSR1; % Spinning reserve MCPs
    MCP_3 = MCP_NSR1; % Non-spinning reserve MCPs

    STOR_spin=x_ed(st_STOR_spin:st_STOR_spin+num_res1_STOR-1);
    STOR_nonspin2=x_ed(obj_Bn+num_gen:obj_Bn+num_gen+num_res1_STOR-1);% oCT 09 2011
    STOR_upreg=Upreg((st_STOR_reg-obj_RU)+1:end);
    STOR_downreg=Downreg((st_STOR_reg-obj_RU)+1:end);
    STOR_Comspin=x_ed(obj_comA:obj_comA+num_comA-1);
    STOR_Comupreg=x_ed(obj_comRU:obj_comRU+num_comRU-1);
    STOR_Comdownreg=x_ed(obj_comRD:obj_comRD+num_comRD-1);


    STOR_energy_profit_2(i)=(STOR_dischar(1:dpi)'*LMP_2-STOR_charge(1:dpi)'*LMP_2);
    STOR_energy_profit_21(i)=(STOR_dischar(dpi+1:end)'*LMP_21-STOR_charge(dpi+1:end)'*LMP_21);
    STOR_spin_profit_2(i)=STOR_spin(1:dpi)'*MCP_1;
    STOR_spin_profit_21(i)=STOR_spin(dpi+1:end)'*MCP_1;
    STOR_nonspin2_profit_2(i)=STOR_nonspin2(1:dpi)'*MCP_3;
    STOR_nonspin2_profit_21(i)=STOR_nonspin2(dpi+1:end)'*MCP_3; 
    STOR_upreg_profit_21(i)=STOR_upreg(dpi+1:end)'*MCP_ru;
    STOR_downreg_profit_21(i)=STOR_downreg(dpi+1:end)'*MCP_rd;    
    STOR_Comspin_profit_21(i)=STOR_Comspin(dpi+1:end)'*MCP_1;    
    STOR_Comupreg_profit_21(i)=STOR_Comupreg(dpi+1:end)'*MCP_ru;
    STOR_Comdownreg_profit_21(i)=STOR_Comdownreg(dpi+1:end)'*MCP_rd;
    STOR_ancillary_profit_21(i)=STOR_spin_profit_21(i)+STOR_nonspin2_profit_21(i)+STOR_upreg_profit_21(i)+...
        STOR_downreg_profit_21(i)+STOR_Comupreg_profit_21(i)+STOR_Comdownreg_profit_21(i)+STOR_Comspin_profit_21(i);


    for kt=1:tt_tmp
        dual_tmp_wind((kt-1)*expand_num+1:kt*expand_num)=v_ed(nst_tran+(Wind_tran(kt,1)-1)*expand_num:nst_tran+Wind_tran(kt,1)*expand_num-1);
    end
end

primal_mean_UC=xx_UC/n;
primal_mean_ED=xx_ED/n;
dual_mean_ED=v_k_ED/n;
optimalv_mean_ED=op_c_ED/n;
optimalv_mean_UC=op_c_UC/n;

STOR_charge=primal_mean_ED(st_com:st_com+num_com-1);
STOR_dischar=primal_mean_ED(st_tur:st_tur+num_tur-1);
STOR_strlevel=[primal_mean_ED(st_res:st_res+dpi-2);primal_mean_ED(obj_final_st);primal_mean_ED(st_res+dpi-1:st_res+num_res-1);primal_mean_ED(obj_final_st+1)];
Wind_output=primal_mean_ED(st_wind:st_wind+num_wind+num_wind_new-1);
LMP=dual_mean_ED(nst_tran:24 * dpi); % 48 hour LMPs from bus 1 to n in 'n' bus system
LMP_2=dual_mean_ED(nst_tran+(STOR1_n - 1)*dpi:STOR1_n * dpi);
LMP_21=dual_mean_ED(nst_tran+(STOR2_n - 1)*dpi:STOR2_n * dpi);

    MCP_SR=dual_mean_ED(obj_row_1:obj_row_1 + dpi-1);
    MCP_NSR=dual_mean_ED(obj_row_3:obj_row_3 + dpi-1);
    MCP_R=dual_mean_ED(obj_row_ru:obj_row_ru + dpi-1);
    MCP_D=dual_mean_ED(obj_row_rd:obj_row_rd + dpi-1);

    MCP_ru = MCP_R+MCP_SR+MCP_NSR;
    MCP_rd = MCP_D;
    MCP_1 = MCP_SR+MCP_NSR;
    MCP_3 = MCP_NSR;
    
gen_avg=primal_mean_ED(st_res+num_res:st_res+num_res-1+dpi*14);

Wactual = wind_num_all.*upbd.*wind_gen(1:num_wind+num_wind_new);
Wspillage = Wactual - Wind_output;
Wspillagep = sum(Wspillage)*100/sum(Wactual);

figure
plot(1:1:dpi,STOR_charge(dpi+1:end),'r');
hold;
[AX,H1,H2]=plotyy(1:1:dpi,-STOR_dischar(dpi+1:end),1:1:dpi,LMP_21);
title('STOR output');
xlabel('Hours')
ylabel('MW')
legend('Charging','discharging','LMP');


STOR_energy_profit_2_avg=mean(STOR_energy_profit_2);
STOR_energy_profit_21_avg=mean(STOR_energy_profit_21);
STOR_ancillary_profit_21_avg=mean(STOR_ancillary_profit_21);

Load_unmet = primal_mean_ED(obj_loadunmet:obj_loadunmet+num_loadunmet-1);
for lc = 1:a*b %hours
    L_unmet(lc,1)=0;
    for lci = 1:24 % nodes
        L_unmet(lc,1) = L_unmet(lc,1)+Load_unmet(((lci-1)*a*b)+lc,1); %node 1 hour 1
    end
end

for windc = 1:a*b %hours
    Wind_spill(windc,1)=0;
    for wci = 1:3 % wind plants
        Wind_spill(windc,1) = Wind_spill(windc,1)+Wspillage(((wci-1)*a*b)+windc,1); %total wind 1 hour 1
    end
end

for windc = 1:a*b %hours
    Wind_sup(windc,1)=0;
    for wci = 1:3 % wind plants
        Wind_sup(windc,1) = Wind_sup(windc,1)+Wind_output(((wci-1)*a*b)+windc,1); %total wind 1 hour 1
    end
end

% houw much of each service storaeg supplies
STORspinp = sum(STOR_spin)*100/sum(spin);
STORupregp = (sum(STOR_upreg)+sum(STOR_Comupreg))*100/(sum(Upreg)+sum(STOR_Comupreg));
STORdownregp = (sum(STOR_downreg)+sum(STOR_Comdownreg))*100/(sum(Downreg)+sum(STOR_Comdownreg));
STORnonspinp = (sum(STOR_nonspin2))*100/(sum(nonspin2));

% Hourly total MW generation from each generation type
for gi = 1:a*b,
    Gen_Coal(gi,1) = gen_avg(a*b*(2-1)+gi,1)+gen_avg(a*b*(4-1)+gi,1)+gen_avg(a*b*(8-1)+gi,1)+gen_avg(a*b*(9-1)+gi,1)+gen_avg(a*b*(12-1)+gi,1)+gen_avg(a*b*(13-1)+gi,1)+gen_avg(a*b*(14-1)+gi,1);
    Gen_Nuc(gi,1) = gen_avg(a*b*(10-1)+gi,1)+gen_avg(a*b*(11-1)+gi,1);
    Gen_NG(gi,1) = gen_avg(a*b*(5-1)+gi,1)+gen_avg(a*b*(6-1)+gi,1)+gen_avg(a*b*(7-1)+gi,1);
    Gen_Oil(gi,1) = gen_avg(a*b*(1-1)+gi,1)+gen_avg(a*b*(3-1)+gi,1);
end

toc