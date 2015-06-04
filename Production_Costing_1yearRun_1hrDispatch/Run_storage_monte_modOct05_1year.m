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
n=1;  % Number of monte carlo simulation
jmp= 0; % round(30*rand(1,1))
kf = 1; %rand(1,1); % ancillary deployment
ke = 0.000001 ; %0.2*rand(1,1);
tolerance=0.001;
count=0;
wc=0;
lead_time=1;
carbon_tax=0; % $/ton of carbob dioxide emission
w=1;
h = 24;
d = 2;
a = 2;
b = 24;

dpi = w*d*h; %data per item for simulation period with 10 min intervals
Caes1_n = 2;
Caes2_n = 21;% different CAES positions
%        UT = 1;
%        DT = 1;
load UpDn
% UT=ones(dpi*16,1);
% DT=ones(dpi*16,1);
spin_bid =[4.2000,4.2000,4.2000,4.2000,4.2000,4.2000,4.2000,5.6000,5.6000,5.6000,5.6000,7.,7.0000,7.0000,...
    7.0000,7.000,5.6000,5.6000,6.3000,6.3000,6.3000,6.3000,4.9000,4.2000,4.2000,4.2000,4.2000,4.2000,4.2000,...
    4.2000,4.2000,5.6000,5.6000,5.6000,5.6000,7.0000,7.0000,7.0000,7.0000,7.0000,5.6000,5.6000,6.3000,6.3000,6.3000,6.3000,4.9000,4.2000];
nonspin_bid=[2.4;2.4;2.4;2.4;2.4;2.4;2.4;3.2;3.20;3.20;3.2;4;4;4;4;4;3.2;3.2;3.60;3.6;3.60;3.60;2.80;2.40;2.40;2.4;2.40;2.4;2.4;2.4;...
    2.4;3.2;3.2;3.20;3.200;4;4;4;4;4;3.2;3.2;3.6;3.6;3.60;3.60;2.800;2.40]';
reg_bid = [7.87,9.50,9.58,9.50,9.50,10.33,10.13,10.29,9.51,7.38,7.17,7.30,8.23,7.54,9.32,9.06,10.21,7.72,...
           6.99,6.99,6.81,8.88,9.83,10.21,7.49,8.45,9.50,9.50,9.50,9.50,16.49,20,11.68,9.19,13.68,9.51,9.31,...
           9.31,9.31,7.54,7.20,9.53,9.43,9.88,7.49,6.99,8.92,9.31];
% Wind_tran=[1,2]';
% cap_factor=[0.3,0.35,0.3]';  % one cap. factor for each wind farm
% Davg=[1250,1250,1250,1250,1250]'; % Average demand for each year
% Pene_target=[0,0.05,0.10,0.15,0.2]'; % Wind penetration level for each year

% Monte Carlo simulation is executed
solu_primal=0;
solu_dual=0;
B_rating=0;  % rating of branch if there is outage
G_rating=0;
num_each=w*b*a;

% Deal with the oil prices
price1(1,1:n)= 0.1*0.3480; % randn(1,n);% 
% Deal with the natural gas prices
price2(1,1:n)= 0.05*-0.7197; %randn(1,n);% 


% Remember: X is determined, then everything evolves around X.
flag_mt=0;
Profit_mean=0;

% Feb-06-2013
ii=1;
load Yearly 
load Reg_req
Yregh = NLVhr;

Slave_UC_1year; %_wNSp_NominON;  % run the UC slave problem
% Slave_UC_No_ReservoirAnc_1year; % Existing model
%%%%%%%%%%%%%%%% update the wind generation information%%%%%%%%%%%%%%%%%%%%
% the upper bound of wind should be determined by wind speed, so it need to be updated

    % First run the UC problem.
    
    %% new Aug 05 2012
    %load Windfc % Apr 5 - wind hourly
load Yearlywind
wind_gen_y=wind_gen;

Slave_ED_1year; %_wNSp_NominON;
% Slave_ED_No_ReservoirAnc_1year; % Existing model

for ii=1:n  % monte carlo simulation for n times
F_ED=F_ED1;c_ED=c_ED1;A_ED=A_ED1;b_L_ED=b_L_ED1;b_U_ED=b_U_ED1;x_L_ED=x_L_ED1;x_U_ED=x_U_ED1;IntVars_ED=IntVars_ED1;
F_UC=F_UC1;c_UC=c_UC1;A_UC=A_UC1;b_L_UC=b_L_UC1;b_U_UC=b_U_UC1;x_L_UC=x_L_UC1;x_U_UC=x_U_UC1;IntVars_UC=IntVars_UC1;

% load update
loadhr = Yloadh(1+(jmp+ii-1)*dpi:dpi*(jmp+ii),1)*3000/(max(Yloadh));% new
for numLd=1:length(stressf)
    x_L_UC(st_dm+(numLd-1)*dpi:st_dm+numLd*dpi-1,1)=0.6*loadhr*stressf(numLd,1);
    x_U_UC(st_dm+(numLd-1)*dpi:st_dm+numLd*dpi-1,1)=0.6*loadhr*stressf(numLd,1);
    x_L_ED(st_dm+(numLd-1)*dpi:st_dm+numLd*dpi-1,1)=0.6*loadhr*stressf(numLd,1);
    x_U_ED(st_dm+(numLd-1)*dpi:st_dm+numLd*dpi-1,1)=0.6*loadhr*stressf(numLd,1);
end

% regulation update
% NLVhr = NLVhravg(1+(ii-1)*dpi:dpi*ii,1);%% new
NLVhr = Yregh(1+(jmp+ii-1)*dpi:dpi*(jmp+ii),1);%% new
b_L_ED(obj_row_1:obj_row_1+dpi-1,1)=0.5*Max_gen+3*NLVhr;%'1'-LB
b_L_ED(obj_row_2:obj_row_2+dpi-1,1)=3*NLVhr;%'2'-LB
b_L_ED(obj_row_3:obj_row_3+dpi-1,1)=Max_gen+3*NLVhr;%'3'-LB
b_L_ED(obj_row_4:obj_row_4+dpi-1,1)=3*NLVhr;%'4'-LB
b_L_ED(obj_row_ru:obj_row_ru+dpi-1,1)=3*NLVhr;%'RU'-LB
b_L_ED(obj_row_rd:obj_row_rd+dpi-1,1)=3*NLVhr;%'RD'-LB - 48 values of hours
b_U_ED(obj_row_DRlim:obj_row_DRlim+dpi-1,1)=3*NLVhr;%'5'-UB(-inf to 3*NLVhHr)-1st storage*hours
b_U_ED(obj_row_DRlim+dpi:obj_row_DRlim+dpi*2-1,1)=3*NLVhr;%'5'-UB(-inf to 3*NLVhHr)-2nd storage*hours

UC_ED_diff = row_count-row_count_ED+dpi*2;%dpi*2 extra in ED for '6' Oct 10 addition (Comp DR+Tur UR<eff*rating, if simultaneous)

b_L_UC(UC_ED_diff+obj_row_1:UC_ED_diff+obj_row_1+dpi-1,1)=0.5*Max_gen+3*NLVhr;%'1'-LB
b_L_UC(UC_ED_diff+obj_row_2:UC_ED_diff+obj_row_2+dpi-1,1)=3*NLVhr;%'2'-LB
b_L_UC(UC_ED_diff+obj_row_3:UC_ED_diff+obj_row_3+dpi-1,1)=Max_gen+3*NLVhr;%'3'-LB
b_L_UC(UC_ED_diff+obj_row_4:UC_ED_diff+obj_row_4+dpi-1,1)=3*NLVhr;%'4'-LB
b_L_UC(UC_ED_diff+obj_row_ru:UC_ED_diff+obj_row_ru+dpi-1,1)=3*NLVhr;%'RU'-LB
b_L_UC(UC_ED_diff+obj_row_rd:UC_ED_diff+obj_row_rd+dpi-1,1)=3*NLVhr;%'RD'-LB - 48 values of hours
b_U_UC(UC_ED_diff+obj_row_DRlim:UC_ED_diff+obj_row_DRlim+dpi-1,1)=3*NLVhr;%'5'-UB(-inf to 3*NLVhHr)-1st storage*hours
b_U_UC(UC_ED_diff+obj_row_DRlim+dpi:UC_ED_diff+obj_row_DRlim+dpi*2-1,1)=3*NLVhr;%'5'-UB(-inf to 3*NLVhHr)-2nd storage*hours

wind_gen = [wind_gen_y(1+(jmp+ii-1)*dpi:dpi*(jmp+ii),1);wind_gen_y(1+(jmp+ii-1)*dpi:dpi*(jmp+ii),1);wind_gen_y(1+(jmp+ii-1)*dpi:dpi*(jmp+ii),1)];
wo=wind_gen(1:num_wind+num_wind_new);
upbd=x_U_UC(st_wind:st_wind+num_wind+num_wind_new-1);
for i2=1:length(upbd)
    if upbd(i2)>10^6
        upbd(i2)=10^6;
    end
    if upbd(i2)==0
        upbd(i2)=wc;
    end
end

wind_num_all=zeros(num_wind+num_wind_new,1);
wind_num_all(1:num_wind)=1;
%

    dual_tmp_wind=zeros(num_wind_new,1); % the duel of wind generation node
    tt_tmp=num_wind_new/b/a; % num. of new wind gen

    
x_U_UC(st_wind:st_wind+num_wind+num_wind_new-1)=wind_num_all.*upbd.*wind_gen(1:num_wind+num_wind_new); % update the upper bound of wind capacity at each hour
x_U_ED(st_wind:st_wind+num_wind+num_wind_new-1)=wind_num_all.*upbd.*wind_gen(1:num_wind+num_wind_new); % update the upper bound of wind capacity at each hour

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

%%
    fot_rate=zeros(num_gen,1);
    Branch_rate=zeros(num_br+num_new,1);
    for k1=st_gen:st_gen + num_gen -1 % Forced-outage rate
        fot_rate(k1-st_gen+1,1)=str2double(cell2mat(arcstmp{k1,29}));
    end
    
   load udf % Gen_udf and Branch_udf conatant
%   Gen_udf=rand(num_gen,1); % For generators
   Gen_udf(433:528,1)=1; %% Nov 5 2011 - Nuclear never outage
%   Branch_udf=rand(num_br+num_new,1); % for branches
%     for j=1:num_gen
%         if Gen_udf(j,1)>=fot_rate(num_gen,1)
%             Gen_udf(j,1)=1;
%         else
%             Gen_udf(j,1)=G_rating;
%         end
%     end
    for k1=st_br:(st_br + num_br+num_new -1)
        Branch_rate(k1-st_br+1,1)=str2double(cell2mat(arcstmp{k1,29}));
    end
%     for j=1:num_br+num_new
%         if Branch_udf(j,1)>=Branch_rate(j,1)
%             Branch_udf(j,1)=1;
%         else
%             Branch_udf(j,1)=B_rating;
%         end
%     end
    
    % the decomposition function changes x_L_UC and x_U_UC
    mt_c_UC=c_UC; % Prices
    mt_c_UC(1:num_each,1)=mt_c_UC(1:num_each,1)+price1(1,ii)*mt_c_UC(1:num_each,1);
    mt_c_UC(num_each*2+1:num_each*3,1)=mt_c_UC(num_each*2+1:num_each*3,1)+price1(1,ii)*mt_c_UC(num_each*2+1:num_each*3,1);
    mt_c_UC(num_each*4+1:num_each*5,1)=mt_c_UC(num_each*4+1:num_each*5,1)+price2(1,ii)*mt_c_UC(num_each*4+1:num_each*5,1);
    mt_c_UC(num_each*5+1:num_each*6,1)=mt_c_UC(num_each*5+1:num_each*6,1)+price2(1,ii)*mt_c_UC(num_each*5+1:num_each*6,1);
    mt_c_UC(num_each*6+1:num_each*7,1)=mt_c_UC(num_each*6+1:num_each*7,1)+price2(1,ii)*mt_c_UC(num_each*6+1:num_each*7,1);
    mt_x_L_UC=x_L_UC;
    mt_x_L_UC(st_gen:st_gen+num_gen-1,1)=mt_x_L_UC(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_L_UC(st_br:st_br+num_br+num_new-1,1)=mt_x_L_UC(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_L_UC(obj_A:obj_A+num_gen-1,1)=mt_x_L_UC(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
%     mt_x_L_UC(obj_B:obj_B+num_gen-1,1)=mt_x_L_UC(obj_B:obj_B+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_L_UC(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_RU:obj_RU+num_gen-1,1)=mt_x_L_UC(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_UC(obj_RD:obj_RD+num_gen-1,1)=mt_x_L_UC(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_UC=x_U_UC;
    mt_x_U_UC(st_gen:st_gen+num_gen-1,1)=mt_x_U_UC(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_U_UC(st_br:st_br+num_br+num_new-1,1)=mt_x_U_UC(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_U_UC(obj_A:obj_A+num_gen-1,1)=mt_x_U_UC(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
%     mt_x_U_UC(obj_B:obj_B+num_gen-1,1)=mt_x_U_UC(obj_B:obj_B+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
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
    
    Xst(:,ii) = x_uc(obj_X:obj_Y-1,1); % Start up status
    Yst(:,ii) = x_uc(obj_Y:obj_Y+(obj_Y-obj_X)-1,1); % Shut down status
    
    Xnst(:,ii) = x_uc(obj_Xn:obj_Yn-1,1); % Start up status NSP Oct 21
    Ynst(:,ii) = x_uc(obj_Yn:obj_Yn+(obj_Yn-obj_Xn)-1,1); % Shut down status NSP Oct 21
    
    Nsp_status = zeros(length(obj_UCNsp_num),1); %% Nsp status
    Nsp_status = x_uc(obj_UCNsp:obj_UCNsp+obj_UCNsp_num-1,1);  
    
    Unsp = UC_status + Nsp_status;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % the decomposition function changes x_L_ED and x_U_ED
    mt_c_ED=c_ED; % Prices
    mt_c_ED(1:num_each,1)=mt_c_ED(1:num_each,1)+price1(1,ii)*mt_c_ED(1:num_each,1);
    mt_c_ED(num_each*2+1:num_each*3,1)=mt_c_ED(num_each*2+1:num_each*3,1)+price1(1,ii)*mt_c_ED(num_each*2+1:num_each*3,1);
    mt_c_ED(num_each*4+1:num_each*5,1)=mt_c_ED(num_each*4+1:num_each*5,1)+price2(1,ii)*mt_c_ED(num_each*4+1:num_each*5,1);
    mt_c_ED(num_each*5+1:num_each*6,1)=mt_c_ED(num_each*5+1:num_each*6,1)+price2(1,ii)*mt_c_ED(num_each*5+1:num_each*6,1);
    mt_c_ED(num_each*6+1:num_each*7,1)=mt_c_ED(num_each*6+1:num_each*7,1)+price2(1,ii)*mt_c_ED(num_each*6+1:num_each*7,1);
    mt_x_L_ED=x_L_ED;
    mt_x_L_ED(st_gen:st_gen+num_gen-1,1)=mt_x_L_ED(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_L_ED(st_br:st_br+num_br+num_new-1,1)=mt_x_L_ED(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_L_ED(obj_A:obj_A+num_gen-1,1)=mt_x_L_ED(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
%     mt_x_L_ED(obj_B:obj_B+num_gen-1,1)=mt_x_L_ED(obj_B:obj_B+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_L_ED(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_RU:obj_RU+num_gen-1,1)=mt_x_L_ED(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_L_ED(obj_RD:obj_RD+num_gen-1,1)=mt_x_L_ED(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED=x_U_ED;
    mt_x_U_ED(st_gen:st_gen+num_gen-1,1)=mt_x_U_ED(st_gen:st_gen+num_gen-1,1).*Gen_udf(:,1);
    mt_x_U_ED(st_br:st_br+num_br+num_new-1,1)=mt_x_U_ED(st_br:st_br+num_br+num_new-1,1).*Branch_udf(:,1);
    mt_x_U_ED(obj_A:obj_A+num_gen-1,1)=mt_x_U_ED(obj_A:obj_A+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
%     mt_x_U_ED(obj_B:obj_B+num_gen-1,1)=mt_x_U_ED(obj_B:obj_B+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_Bn:obj_Bn+num_gen-1,1)=mt_x_U_ED(obj_Bn:obj_Bn+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_RU:obj_RU+num_gen-1,1)=mt_x_U_ED(obj_RU:obj_RU+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_x_U_ED(obj_RD:obj_RD+num_gen-1,1)=mt_x_U_ED(obj_RD:obj_RD+num_gen-1,1).*Gen_udf(:,1); % Oct 09 2011
    mt_b_L_ED=b_L_ED; %% Oct 20 2011
    mt_b_L_ED(obj_row_S_pter_ED:obj_row_S_pter_ED+num_gen-1,1)=mt_b_L_ED(obj_row_S_pter_ED:obj_row_S_pter_ED+num_gen-1,1).*Gen_udf(:,1);
    mt_b_U_ED=b_U_ED; %% Oct 20 2011

    
    % After updating the price and gen/tran availability, update the unit commitment decision.
    % Didn't update the UC decisions about storage
    kng = 0;
    for kn=0:obj_UC_num-1  % To deal with the Inf*0 problem
        %if kn<num_gen % short term
        if (x_uc(obj_UC+kn))<0.5
            kng = kn;
            if kn>=num_gen %% wind gen. starts, skip until storage gen.
                kng = kn+num_wind;
            end               
                
%            mt_x_U_ED(obj_A-2*obj_UC_num+kn)=0; %% changed... X and Y at the end
            mt_x_U_ED(obj_A+kn)=0; %% changed... X and Y at the end
            mt_x_U_ED(st_gen+kng)=0;
%             mt_x_U_ED(obj_B+kn)=0; 
            mt_x_U_ED(obj_RU+kn)=0;
            mt_x_U_ED(obj_RD+kn)=0;
%            mt_x_U_ED(obj_Bn+kn)=0;           

            mt_x_L_ED(obj_A+kn)=0; %% changed... X and Y at the end
            mt_x_L_ED(st_gen+kng)=0;
%             mt_x_L_ED(obj_B+kn)=0; 
            mt_x_L_ED(obj_RU+kn)=0;
            mt_x_L_ED(obj_RD+kn)=0;
 %           mt_x_L_ED(obj_Bn+kn)=0;           
            
            mt_b_L_ED(obj_row_S_pter_ED+kn)=0;           
        end
        if (x_uc(obj_UCNsp+kn))<0.5
%            mt_x_U_ED(obj_A-2*obj_UC_num+kn)=0; %% changed... X and Y at the end
            mt_x_U_ED(obj_Bn+kn)=0;                   
            mt_x_L_ED(obj_Bn+kn)=0;           
        end
        %end
    end
    
       
    %% June 27 2012 - Disjoint Charge/Discharge Operation in Storage
    for kn=0:num_tur-1
        if (x_uc(obj_UC+obj_UC_num-num_tur+kn)+x_uc(obj_UCNsp+obj_UC_num-num_tur+kn))>0.5
            mt_x_U_ED(st_com+kn)=0;
            mt_x_L_ED(st_com+kn)=0;
%             mt_b_L_ED(obj_row_StorComLDJ+kn)=0; %OCT 12 2012
            mt_b_U_ED(obj_row_StorComLDJ+kn)=0;            
        end
    end
    
%      for kn=0:num_gen-1  % To deal with the Inf*0 problem
%         if (x_uc(obj_UC+kn))==0
%             mt_x_U_ED(obj_A-2*obj_UC_num+kn)=0;
%             mt_x_U_ED(st_gen+kn)=0;
%             mt_x_U_ED(obj_RU-2*obj_UC_num+kn)=0; %% Oct 03 2011 
%             mt_x_U_ED(obj_RD-2*obj_UC_num+kn)=0; %% Oct 03 2011
%             mt_x_L_ED(obj_A-2*obj_UC_num+kn)=0;
%             mt_x_L_ED(st_gen+kn)=0;
%             mt_x_L_ED(obj_RU-2*obj_UC_num+kn)=0; %% Oct 03 2011 
%             mt_x_L_ED(obj_RD-2*obj_UC_num+kn)=0; %% Oct 03 2011
%         end
%         if(x_uc(obj_UCNsp+kn))==0 %% Oct 03 2011 Nsp
%             mt_x_U_ED(obj_B-2*obj_UC_num+kn)=0;
%             mt_x_L_ED(obj_B-2*obj_UC_num+kn)=0;            
%         end
%      end
    
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
%     CAES_charge=x_ed(817:912);
%     CAES_dischar=x_ed(1823:1918);
%     CAES_strlevel=x_ed(913:1006);
%     Wind_output=x_ed(1679:1822);
    PrCosti(ii,1)=objv_ed;
    CAES_chargei(:,ii)=x_ed(st_com:st_com+num_com-1);
    CAES_dischari(:,ii)=x_ed(st_tur:st_tur+num_tur-1);
    CAES_strleveli(:,ii)=[x_ed(st_res:st_res+dpi-2);x_ed(obj_final_st);x_ed(st_res+dpi-1:st_res+num_res-1);x_ed(obj_final_st+1)];
    
    CAES_charge=x_ed(st_com:st_com+num_com-1);
    CAES_dischar=x_ed(st_tur:st_tur+num_tur-1);
    CAES_strlevel=[x_ed(st_res:st_res+dpi-2);x_ed(obj_final_st);x_ed(st_res+dpi-1:st_res+num_res-1);x_ed(obj_final_st+1)];
    
    Wind_output=x_ed(st_wind:st_wind+num_wind+num_wind_new-1);
    gen=x_ed(st_res+num_res:st_res+num_res-1+dpi*14);
    spin=x_ed(st_spin:st_spin+num_res1-1);
%     nonspin=x_ed(st_spin+num_res1:st_spin+2*num_res1-1);
    nonspin2=x_ed(obj_Bn:obj_Bn+num_res1-1); % oCT 09 2011
    Upreg=x_ed(obj_RU:obj_RU+num_reg-1);
    Downreg=x_ed(obj_RD:obj_RD+num_reg-1);

%     LMP_2=v_ed(49:96);
%     LMP_21=v_ed(961:1008);
    LMP_2=v_ed(nst_tran+(Caes1_n - 1)*dpi:Caes1_n * dpi);
    LMP_21=v_ed(nst_tran+(Caes2_n - 1)*dpi:Caes2_n * dpi);
    LMP_21i(:,ii)=v_ed(nst_tran+(Caes2_n - 1)*dpi:Caes2_n * dpi);
    MCP_SR2=v_ed(obj_row_2:obj_row_2 + dpi-1);
    MCP_NSR2=v_ed(obj_row_4:obj_row_4 + dpi-1);
    MCP_SR1=v_ed(obj_row_1:obj_row_1 + dpi-1);
    MCP_NSR1=v_ed(obj_row_3:obj_row_3 + dpi-1);
    MCP_R=v_ed(obj_row_ru:obj_row_ru + dpi-1);
    MCP_D=v_ed(obj_row_rd:obj_row_rd + dpi-1);

    MCP_rui(:,ii) = MCP_R+MCP_SR1+MCP_NSR1;
    MCP_rdi(:,ii) = MCP_D;
    MCP_1i(:,ii) = MCP_SR1+MCP_NSR1;
    MCP_3i(:,ii) = MCP_NSR1;
    
    MCP_ru = MCP_R+MCP_SR1+MCP_NSR1;
    MCP_rd = MCP_D;
    MCP_1 = MCP_SR1+MCP_NSR1;
    MCP_3 = MCP_NSR1;

%     CAES_spin=x_ed(7343:7438);
%     CAES_nonspin=x_ed(8111:8206);

    CAES_spini(:,ii)=x_ed(st_caes_spin:st_caes_spin+num_res1_caes-1);
%     CAES_nonspin=x_ed(st_caes_nonspin:st_caes_nonspin+num_res2_caes-1);
    CAES_nonspin2i(:,ii)=x_ed(obj_Bn+num_gen:obj_Bn+num_gen+num_res1_caes-1);% oCT 09 2011
    CAES_upregi(:,ii)=Upreg((st_caes_reg-obj_RU)+1:end);
    CAES_downregi(:,ii)=Downreg((st_caes_reg-obj_RU)+1:end);
    CAES_Comspini(:,ii)=x_ed(obj_comA:obj_comA+num_comA-1);
    CAES_Comupregi(:,ii)=x_ed(obj_comRU:obj_comRU+num_comRU-1);
    CAES_Comdownregi(:,ii)=x_ed(obj_comRD:obj_comRD+num_comRD-1);

    CAES_spin=x_ed(st_caes_spin:st_caes_spin+num_res1_caes-1);
%     CAES_nonspin=x_ed(st_caes_nonspin:st_caes_nonspin+num_res2_caes-1);
    CAES_nonspin2=x_ed(obj_Bn+num_gen:obj_Bn+num_gen+num_res1_caes-1);% oCT 09 2011
    CAES_upreg=Upreg((st_caes_reg-obj_RU)+1:end);
    CAES_downreg=Downreg((st_caes_reg-obj_RU)+1:end);
    CAES_Comspin=x_ed(obj_comA:obj_comA+num_comA-1);
    CAES_Comupreg=x_ed(obj_comRU:obj_comRU+num_comRU-1);
    CAES_Comdownreg=x_ed(obj_comRD:obj_comRD+num_comRD-1);
    

%     CAES_energy_profit_2(ii)=365*(CAES_dischar(1:24)'*LMP_2(1:24)-CAES_charge(1:24)'*LMP_2(1:24));
%     CAES_energy_profit_21(ii)=365*(CAES_dischar(49:72)'*LMP_21(1:24)-CAES_charge(49:72)'*LMP_21(1:24));
%     CAES_spin_profit_2(ii)=365*spin_bid(1,1:24)*CAES_spin(1:24);
%     CAES_spin_profit_21(ii)=365*spin_bid(1,1:24)*CAES_spin(49:72);
%     CAES_nonspin_profit_2(ii)=365*nonspin_bid(1,1:24)*CAES_nonspin(1:24);
%     CAES_nonspin_profit_21(ii)=365*nonspin_bid(1,1:24)*CAES_nonspin(49:72);
%     CAES_ancillary_profit_2(ii)=CAES_spin_profit_2(ii)+CAES_nonspin_profit_2(ii);
%     CAES_ancillary_profit_21(ii)=CAES_spin_profit_21(ii)+CAES_nonspin_profit_21(ii);
    CAES_energy_profit_2(ii,1)=(CAES_dischar(1:dpi)'*LMP_2-CAES_charge(1:dpi)'*LMP_2);
    CAES_energy_profit_21(ii,1)=(CAES_dischar(dpi+1:end)'*LMP_21-CAES_charge(dpi+1:end)'*LMP_21);
    %CAES_spin_profit_2(ii)=spin_bid(1,1:dpi)*CAES_spin(1:dpi);
    CAES_spin_profit_2(ii,1)=CAES_spin(1:dpi)'*MCP_1;
%    CAES_spin_profit_21(ii)=spin_bid(1,1:dpi)*CAES_spin(dpi+1:end);
    CAES_spin_profit_21(ii,1)=CAES_spin(dpi+1:end)'*MCP_1;
%     CAES_nonspin_profit_2(ii)=365*nonspin_bid(1,1:24)*CAES_nonspin(1:24);
%     CAES_nonspin_profit_21(ii)=365*nonspin_bid(1,1:24)*CAES_nonspin(49:72);
%     CAES_nonspin_profit_2(ii)=CAES_nonspin(1:dpi)'*MCP_3;
%     CAES_nonspin_profit_21(ii)=CAES_nonspin(dpi+1:end)'*MCP_3;
    CAES_nonspin2_profit_2(ii,1)=CAES_nonspin2(1:dpi)'*MCP_3;
    CAES_nonspin2_profit_21(ii,1)=CAES_nonspin2(dpi+1:end)'*MCP_3;
    
    
    %CAES_upreg_profit_21(ii)=reg_bid(1,1:dpi)*CAES_upreg(dpi+1:end);
    CAES_upreg_profit_21(ii,1)=CAES_upreg(dpi+1:end)'*MCP_ru;
%    CAES_downreg_profit_21(ii)=reg_bid(1,1:dpi)*CAES_downreg(dpi+1:end);
    CAES_downreg_profit_21(ii,1)=CAES_downreg(dpi+1:end)'*MCP_rd;
    
    %CAES_Comspin_profit_21(ii)=spin_bid(1,1:dpi)*CAES_Comspin(dpi+1:end);
    CAES_Comspin_profit_21(ii,1)=CAES_Comspin(dpi+1:end)'*MCP_1;
    
 %   CAES_Comupreg_profit_21(ii)=reg_bid(1,1:dpi)*CAES_Comupreg(dpi+1:end);
   CAES_Comupreg_profit_21(ii,1)=CAES_Comupreg(dpi+1:end)'*MCP_ru;
   
% CAES_Comdownreg_profit_21(ii)=reg_bid(1,1:dpi)*CAES_Comdownreg(dpi+1:end);
 CAES_Comdownreg_profit_21(ii,1)=CAES_Comdownreg(dpi+1:end)'*MCP_rd;
 
%     CAES_ancillary_profit_2(ii)=CAES_spin_profit_2(ii)+CAES_nonspin_profit_2(ii);
    CAES_ancillary_profit_21(ii,1)=CAES_spin_profit_21(ii,1)+CAES_nonspin2_profit_21(ii,1)+CAES_upreg_profit_21(ii,1)+...
        CAES_downreg_profit_21(ii,1)+CAES_Comupreg_profit_21(ii,1)+CAES_Comdownreg_profit_21(ii,1)+CAES_Comspin_profit_21(ii,1);


    for kt=1:tt_tmp
        dual_tmp_wind((kt-1)*expand_num+1:kt*expand_num)=v_ed(nst_tran+(Wind_tran(kt,1)-1)*expand_num:nst_tran+Wind_tran(kt,1)*expand_num-1);
    end
    %     Profit_mean=Profit_mean-c_ED(st_wind:st_wind+num_wind+num_wind_new-1)'*x_ed(st_wind:st_wind+num_wind+num_wind_new-1)+dual_tmp_wind'*x_ed(st_wind:st_wind+num_wind+num_wind_new-1);
    %     ttttmp(ii)=-c_ED(st_wind:st_wind+num_wind+num_wind_new-1)'*x_ed(st_wind:st_wind+num_wind+num_wind_new-1)+dual_tmp_wind'*x_ed(st_wind:st_wind+num_wind+num_wind_new-1);
end
% Profit_mean=Profit_mean/n;
primal_mean_UC=xx_UC/n;
primal_mean_ED=xx_ED/n;
dual_mean_ED=v_k_ED/n;
optimalv_mean_ED=op_c_ED/n;
optimalv_mean_UC=op_c_UC/n;
%%** Just for the case study
% CAES_charge=primal_mean_ED(817:912);
% CAES_dischar=primal_mean_ED(1823:1918);
% CAES_strlevel=primal_mean_ED(913:1006);
% Wind_output=primal_mean_ED(1679:1822);
% LMP_2=dual_mean_ED(49:96);
% LMP_21=dual_mean_ED(961:1008);

CAES_charge=primal_mean_ED(st_com:st_com+num_com-1);
CAES_dischar=primal_mean_ED(st_tur:st_tur+num_tur-1);
% CAES_strlevel=primal_mean_ED(st_res:st_res+num_res-1);
CAES_strlevel=[primal_mean_ED(st_res:st_res+dpi-2);primal_mean_ED(obj_final_st);primal_mean_ED(st_res+dpi-1:st_res+num_res-1);primal_mean_ED(obj_final_st+1)];
Wind_output=primal_mean_ED(st_wind:st_wind+num_wind+num_wind_new-1);
LMP_2=dual_mean_ED(nst_tran+(Caes1_n - 1)*dpi:Caes1_n * dpi);
LMP_21=dual_mean_ED(nst_tran+(Caes2_n - 1)*dpi:Caes2_n * dpi);

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
plot(1:1:dpi,CAES_charge(dpi+1:end),'r');
hold;
[AX,H1,H2]=plotyy(1:1:dpi,-CAES_dischar(dpi+1:end),1:1:dpi,LMP_21);
title('CAES output');
xlabel('Hours')
ylabel('MW')
legend('Charging','discharging','LMP');


% plot(1:1:24,CAES_charge(49:72));
% hold;
% [AX,H1,H2]=plotyy(1:1:24,-CAES_dischar(49:72),1:1:24,LMP_21(1:24));

% plot(1:1:24,CAES_charge(1:24));
% hold;
% [AX,H1,H2]=plotyy(1:1:24,-CAES_dischar(1:24),1:1:24,LMP_2(1:24));

% CAES_energy_profit_2_avg=mean(CAES_energy_profit_2);
% CAES_energy_profit_21_avg=mean(CAES_energy_profit_21);
% CAES_spin_profit_2_avg=mean(CAES_spin_profit_2);
% CAES_spin_profit_21_avg=mean(CAES_spin_profit_21);
% CAES_nonspin_profit_2_avg=mean(CAES_nonspin_profit_2);
% CAES_nonspin_profit_21_avg=mean(CAES_nonspin_profit_21);
% CAES_ancillary_profit_2_avg=mean(CAES_ancillary_profit_2);
% CAES_ancillary_profit_21_avg=mean(CAES_ancillary_profit_21);

CAES_energy_profit_2_avg=mean(CAES_energy_profit_2);
CAES_energy_profit_21_avg=mean(CAES_energy_profit_21);
CAES_ancillary_profit_21_avg=mean(CAES_ancillary_profit_21);

Load_unmet = primal_mean_ED(obj_loadunmet:obj_loadunmet+num_loadunmet-1);
for lc = 1:w*b*a %hours
    L_unmet(lc,1)=0;
    for lci = 1:24 % nodes
        L_unmet(lc,1) = L_unmet(lc,1)+Load_unmet(((lci-1)*w*b*a)+lc,1); %node 1 hour 1
    end
end

for windc = 1:w*b*a %hours
    Wind_spill(windc,1)=0;
    for wci = 1:3 % wind plants
        Wind_spill(windc,1) = Wind_spill(windc,1)+Wspillage(((wci-1)*w*b*a)+windc,1); %total wind 1 hour 1
    end
end

for windc = 1:w*b*a %hours
    Wind_sup(windc,1)=0;
    for wci = 1:3 % wind plants
        Wind_sup(windc,1) = Wind_sup(windc,1)+Wind_output(((wci-1)*w*b*a)+windc,1); %total wind 1 hour 1
    end
end

CAESspinp = sum(CAES_spin)*100/sum(spin);
CAESupregp = (sum(CAES_upreg)+sum(CAES_Comupreg))*100/(sum(Upreg)+sum(CAES_Comupreg));
CAESdownregp = (sum(CAES_downreg)+sum(CAES_Comdownreg))*100/(sum(Downreg)+sum(CAES_Comdownreg));
CAESnonspinp = (sum(CAES_nonspin2))*100/(sum(nonspin2));

for gi = 1:w*b*a,
    Gen_Coal(gi,1) = gen_avg(w*b*a*(2-1)+gi,1)+gen_avg(w*b*a*(4-1)+gi,1)+gen_avg(w*b*a*(8-1)+gi,1)+gen_avg(w*b*a*(9-1)+gi,1)+gen_avg(w*b*a*(12-1)+gi,1)+gen_avg(w*b*a*(13-1)+gi,1)+gen_avg(w*b*a*(14-1)+gi,1);
    Gen_Nuc(gi,1) = gen_avg(w*b*a*(10-1)+gi,1)+gen_avg(w*b*a*(11-1)+gi,1);
    Gen_NG(gi,1) = gen_avg(w*b*a*(5-1)+gi,1)+gen_avg(w*b*a*(6-1)+gi,1)+gen_avg(w*b*a*(7-1)+gi,1);
    Gen_Oil(gi,1) = gen_avg(w*b*a*(1-1)+gi,1)+gen_avg(w*b*a*(3-1)+gi,1);
end

tot_gen = spin+Upreg+[gen_avg;CAES_dischar];
% save CAES3_200_800_200_0W optimalv_mean_ED optimalv_mean_UC CAES_ancillary_profit_21_avg CAES_energy_profit_21_avg CAES_dischar CAES_strlevel Wind_sup LMP_21 Gen_Coal Gen_Nuc Gen_NG Gen_Oil CAESspinp CAESupregp CAESdownregp Wspillagep

%save YearlyWP22CAES100_k7
% %% Mar 09 2012 - Storage Duals
% 
% obj_row_StorTur = obj_row_S_pter_ED-10*dpi*2;
% StorTurC = dual_mean_ED(obj_row_StorTur:obj_row_StorTur+10*dpi*2-1);
% StorTurCS = dual_mean_ED(obj_row_S_pter_ED:obj_row_S_pter_ED+dpi*2-1);
% StorComC = dual_mean_ED(obj_row_StorCom:obj_row_StorCom+8*dpi*2-1);
% StorResC = dual_mean_ED(obj_row_StorRes:obj_row_StorRes+4*dpi*2-1);
% StorDRlimC = dual_mean_ED(obj_row_DRlim:obj_row_DRlim+dpi*2-1);
% 
% %% Stor res
%  count =0;
% for i = 1:96,
% count=count+1;
% StorResCL(i,1)=StorResC(count,1);
% count=count+1;
% StorResCS(i,1)=StorResC(count,1);
% count=count+1;
% StorResCT(i,1)=StorResC(count,1);
% count=count+1;
% StorResCC(i,1)=StorResC(count,1);
% end
% StorResCE=dual_mean_ED(nst_stor:nst_stor+ 2*dpi-1);
% 
% %% Stor Comp
%  count =0;
% for i = 1:96,
% count=count+1;
% StorComCL(i,1)=StorComC(count,1);
% count=count+1;
% StorComCS(i,1)=StorComC(count,1);
% count=count+1;
% StorComCa(i,1)=StorComC(count,1);
% count=count+1;
% StorComCc(i,1)=StorComC(count,1);
% count=count+1;
% StorComCd(i,1)=StorComC(count,1);
% count=count+1;
% StorComCe(i,1)=StorComC(count,1);
% count=count+1;
% StorComCR(i,1)=StorComC(count,1);
% count=count+1;
% StorComCD(i,1)=StorComC(count,1);
% end
% 
% %% Stor Tur
% count =0;
% for i = 1:96,
% count=count+1;
% StorTurCE(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCL(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCM(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCa(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCb(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCR(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCD(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCc(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCd(i,1)=StorTurC(count,1);
% count=count+1;
% StorTurCe(i,1)=StorTurC(count,1);
% end

toc

%%**
% theta=primal_mean(st_dm+num_dm+num_bus:st_dm+num_dm+num_bus+num_bus-1);  % num_bus==node_tran
% lambda=dual_mean(nst_tran:nst_tran+node_tran-1);
% lambda_wind=zeros(num_wind_new,1);
% % lambda_wind_year=zeros(num_wind_new_year,1);
% for i=1:tt_tmp
%     lambda_wind((i-1)*expand_num+1:i*expand_num)=dual_mean(nst_tran+(Wind_tran(i,1)-1)*expand_num:nst_tran+Wind_tran(i,1)*expand_num-1);
% end
% clear;


