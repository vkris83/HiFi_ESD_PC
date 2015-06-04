%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% SCUC (uses loadhourly.mat, Reg_req.mat)
% 1 year run: simulates 2-day runs over a year sequentially

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program for generating MPS format from
%   nodes.txt and arcs.txt
% This is the optimized version and use the 11 properties arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [num_arcs,num_nodes,num_br,num_bus,num_new]= creatmps(N_lines)
% clear;
%**************************************************************************
%**************************************************************************
% C{1,1}='Elec';C{1,2}=[1 12  4 3];C{1,3}={'y' 'm' 'w' 'c'};
% a=1;
% b=24;
C{1,1}='Elec';C{1,2}=[a b];C{1,3}={'d' 'h'};
C{2,1}='Nuc';C{2,2}=[a b];C{2,3}={'d' 'h'};
C{3,1}='Coal';C{3,2}=[a b];C{3,3}={'d' 'h'};
C{4,1}='Oil';C{4,2}=[a b];C{4,3}={'d' 'h'};
C{5,1}='Wind';C{5,2}=[a b];C{5,3}={'d' 'h'};
C{6,1}='Stor';C{6,2}=[a b];C{6,3}={'d' 'h'};
C{6,1}='NG';C{6,2}=[a b];C{6,3}={'d' 'h'};
OR1=350;

T=C{1,2};
W=C{1,3};
Num_st=0;
Num_rhs=0;
Num_ang=0;
t_num=1;
dr=0;
lol=300;
num_x=0;        % dummy fuel resource nodes
num_gen=0;      % num. of generators excluding wind units    
num_fuel=0;     % num. of fuel transportation lines (excuding dummy lines)
num_dm=0;       % num. of demand nodes
num_nodes=0;    % num. of all nodes except dummy nodes
num_arcs=0;     % num. of all arcs
num_br=0;       % number of existing electric branches
% num_fconst=0; % number of flow constraints
num_bus=0;      % num. of all angles/buses
num_new=0;      % num. of all potential lines
num_new_DC=0;
num_wind=0;     % num. of existing wind generators
num_wind_new=0; % num. of new wind generators
num_UC=0;
num_xcew=0;     % Number of wind units

st_arcs=1;
st_x=0;          % num_arcs - num_dm - num_br - num_gen - num_x +1;
st_gen=0;        % st_x + num_x;
st_br=0;         % st_gen + num_gen;
st_dm=0;         % st_br + num_br;
st_new=0;
st_new_DC=0;
st_wind=0;
st_wind_new=0;
st_xcew=0;

Pctg_nowind=0.07; % percentage of non-wind generation that should be secured by reserve
Pctg_wind=0.10;   % percentage of wind generation that should be secured by reserve
Max_gen=350;      % The capacity of the largest gen. unit 
obj_UC_num=0;
obj_pointer=0;
obj_enebid=0; 
obj_UC=0;         % The beginning of the UC variable
obj_X=0;
obj_Y=0;
obj_A=0;          % The beginning of the first reserve1(spin) variable
obj_B=0;          % The beginning of the first reserve2(non-spin) variable 
obj_SL=0;
flag_obj_UC=0;         % The beginning of the UC variable
flag_obj_X=0;
flag_obj_Y=0;
flag_obj_A=0;          % The beginning of the first reserve1(spin) variable
flag_obj_B=0;          % The beginning of the first reserve2(non-spin) variable
flag_obj_SL=0;
obj_variable={};
obj_row={};
obj_count=0; 
row_count=0;

node_wind=0;
node_prod=0;
node_stor=0;
node_tran=0;
node_gen=0;
nst_prod=0;
nst_stor=0;
nst_tran=0;
nst_gen=0;
nst_wind=0;
nst_wind_new=0;

flagres = 0; %reservoir
num_res =0;
st_res = 0;
flagtur = 0;%turbine
num_tur =0;
st_tur = 0;
flagcom = 0;%compressor
num_com =0;
st_com = 0;

expand_num=C{1,2}(1)*C{1,2}(2);

num_itg=0;      % number of integer variable
%**************************************************************************
% New  variables are anfle= AG+nodename(buses), New branches for DC OPF = BR+branchname(transmission)
%**************************************************************************

%   Read nodes.txt
fid=fopen('nodes.txt');
nodes=textscan(fid,'%s %s %s %s',-1);
fclose('all');
for k1=1:length(nodes)
    for k2=1:length(nodes{1})
        nodestmp{k2,k1}=nodes{k1}(k2);
    end;
end;
clear nodes;
%   Read arcs.txt
fid=fopen('arcs.txt');
arcs=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,-1);  % 32 colmns
fclose('all');
for k1=1:length(arcs)
    for k2=1:length(arcs{1})
        arcstmp{k2,k1}=arcs{k1}(k2);
    end;
end;
clear arcs;
%   Read arcsinitial.txt
fid=fopen('arcsinitial.txt');
arcsinitial=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,-1);  % 32 colmns
fclose('all');
for k1=1:length(arcsinitial)
    for k2=1:length(arcsinitial{1})
        arcsinitialtmp{k2,k1}=arcsinitial{k1}(k2);
    end;
end;
clear arcsinitial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Writes problem to a MPS file
% fid=fopen('C:\tomlab\problem_UC.mps', 'w');
fid=fopen('C:\tomlab\problem_UC.mps', 'w');

%   Write section ROWS
fprintf(fid,'NAME IES\nROWS\n N obj\n');
num_arcs=length(arcstmp(:,1)); % num. of all arcs
num_nodes=length(nodestmp(:,1)); % num. of all nodes
Name={};
fg01=0;
fg02=0;
fg03=0;
fg04=0;
fg05=0;
fg06=0;
for k1=1:length(nodestmp(:,1))
    if strncmp(cell2mat(nodestmp{k1,1}),'X',1)==0 % If it is not a dummy node
        fprintf(fid,' E %s\n',cell2mat(nodestmp{k1,1}));
        row_count=row_count+1;
        obj_row{row_count,1}={cell2mat(nodestmp{k1,1})};
        obj_row{row_count,2}={num2str(row_count)};       
    end
    % if (strncmp(cell2mat(nodestmp{k1,3}),'Gen',3)==1)&&(strncmp(cell2mat(nodestmp{k1,2}),'Wind',4)~=1)
    % if it is the first generator arc but not a wind arc
    if strncmp(cell2mat(nodestmp{k1,1}),'EL',2)==1 % If it is elec node
        num_bus=num_bus+1;
        node_tran=node_tran+1;
        if fg01==0
            nst_tran=k1;
            fg01=1;
        end
    end
    if strncmp(cell2mat(nodestmp{k1,3}),'Prod',4)==1 % If it is prod node
        node_prod=node_prod+1;
        if fg02==0
            nst_prod=k1;
            fg02=1;
        end
    end
    if strncmp(cell2mat(nodestmp{k1,1}),'ST',2)==1 % If it is stor node
        node_stor=node_stor+1;
        if fg03==0
            nst_stor=k1; % Shadow P1
            fg03=1;
        end
    end
    if strncmp(cell2mat(nodestmp{k1,3}),'Gen',3)==1 % If it is gen node
        node_gen=node_gen+1;
        if fg04==0
            nst_gen=k1;
            fg04=1;
        end
    end
    if strncmp(cell2mat(nodestmp{k1,2}),'Wind',4)==1 % If it is wind node
        node_wind=node_wind+1;
        if fg06==0
            nst_wind=k1;
            fg06=1;
        end
    end
end
flag01=0;
flag02=0;
flag03=0;
flag04=0;
flag05=0;
flag06=0;
flag07=0;
flag08=0;
flag09=0;
% NEW code, add N equations
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
            &&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==0)  % If it's a transmission line, add N equations   OK
        fprintf(fid,' E %s%s\n','BR',cell2mat(arcstmp{k1,1}));
        num_br=num_br+1;
        Name{num_br,1}=strcat('BR',cell2mat(arcstmp{k1,1}));
        row_count=row_count+1;
        obj_row{row_count,1}={cell2mat(arcstmp{k1,1})};
        obj_row{row_count,2}={num2str(row_count)};    
    end;
    
    if (strncmp(cell2mat(arcstmp{k1,1}),'X',1)==1)
        num_x=num_x+1;
        if (flag01==0)
            st_x=k1;
            flag01=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))&&((strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)~=1))
        num_gen=num_gen+1;
        if (flag02==0)
            st_gen=k1;
            flag02=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)
        
        if (flag03==0)
            st_br=k1;
            flag03=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'DD',2)==1)
        num_dm=num_dm+1;
        if (flag04==0)
            st_dm=k1;
            flag04=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
            &&(str2num(cell2mat(arcstmp{k1,9}))==0)&&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==0)
        num_new=num_new+1;
        if (flag05==0)
            st_new=k1;
            flag05=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
            &&(str2num(cell2mat(arcstmp{k1,9}))==0)&&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==1)
        num_new_DC=num_new_DC+1;
        if (flag08==0)
            st_new_DC=k1;
            flag08=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1)&&(str2double(cell2mat(arcstmp{k1,10}))<10^-5)
        num_wind=num_wind+1;
        if flag06==0
        st_wind=k1;
        flag06=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1)&&(str2double(cell2mat(arcstmp{k1,10}))>10^-5)
        num_wind_new=num_wind_new+1;
        if flag09==0
        st_wind_new=k1;
        flag09=1;
        end
    end
    if (strncmp(cell2mat(arcstmp{k1,2}),'XC',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EW',2)==1)
        num_xcew=num_xcew+1;
        if (flag07==0)
            st_xcew=k1;
            flag07=1;
        end
    end
    
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        num_res=num_res+1;
        if (flagres==0)
            st_res=k1;
            flagres=1;
        end
   end
   
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
        num_tur=num_tur+1;
        if (flagtur==0)
            st_tur=k1;
            flagtur=1;
        end
   end
    
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       num_com=num_com+1;
       if (flagcom==0)
           st_com=k1;
           flagcom=1;
       end
   end   
end;
if st_wind==0
    st_wind=st_wind_new;
end
nst_wind_new=nst_wind+num_wind;
num_br=num_br-num_new;  % In this way, num_br is the num. of all existing lines
invcost=zeros(num_new,1); % the investment cost vector
line_num=zeros(num_new,1); % number of lines for each kind of potential lines
line_sus=zeros(num_new,1); % susceptance of each potential line
bound=zeros(num_new,1);   % bound of each potential line
for i_tmp=st_new:st_new+num_new-1
    j_tmp=i_tmp-st_new+1;
    line_sus(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,11}));
    line_num(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,9}));
    invcost(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,10}));
    bound(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,8}));
end


for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc segment
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))
        fprintf(fid,' E %s%s\n',cell2mat(arcstmp{k1,1}),'E');  % this row is defined before. Arcname=Eng.bid 1 + Eng. bid 2 + Eng. bid 3
                row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'E')};
        obj_row{row_count,2}={num2str(row_count)};    
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for reserve 1
                row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
        obj_row{row_count,2}={num2str(row_count)};    
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'M'); 
                row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'M')};
        obj_row{row_count,2}={num2str(row_count)}; 
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'a');  % Max. cap for reserve 1
                row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'a')};
        obj_row{row_count,2}={num2str(row_count)};       
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'b');  % Max. cap for reserve 2
                        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'b')};
        obj_row{row_count,2}={num2str(row_count)}; 
        fprintf(fid,' E %s%s\n',cell2mat(arcstmp{k1,1}),'U');  % This row defines the startup and shutdown costs U(t)-U(t+1)=X-Y
                        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'U')};
        obj_row{row_count,2}={num2str(row_count)}; 
        %% Oct 03 2011 - Non-Spinning
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'N');
                        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'N')};
        obj_row{row_count,2}={num2str(row_count)}; 
        %% Oct 20 2011 - Non-Spinning
        fprintf(fid,' E %s%s\n',cell2mat(arcstmp{k1,1}),'Un');  % This row defines the startup and shutdown costs U(t)-U(t+1)=X-Y
                        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'Un')};
        obj_row{row_count,2}={num2str(row_count)}; 
        %   the above line is differnet from fprintf(fid,'E%s%s\n','E',cell2mat(nodestmp{k1,1}));,which means the KCLfor each node
        

        %% Generator Ramping - Mar 11 2011 - Ramp up & Ramp down
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'R');  %Ramp up
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'R')};
        obj_row{row_count,2}={num2str(row_count)}; 
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'D');  %Ramp down
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'D')};
        obj_row{row_count,2}={num2str(row_count)}; 
        
        %% Generator Regulation - Mar 27 2011 - Reg up & down
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'c');  %Reg up
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'c')};
        obj_row{row_count,2}={num2str(row_count)}; 
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'d');  %Reg down
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'d')};
        obj_row{row_count,2}={num2str(row_count)};
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'e');  %SR+RU+NSR<10*rr
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'e')};
        obj_row{row_count,2}={num2str(row_count)}; % Feb 18 2012 new
        
        %% Nov 4 2011 - min. up time
        fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'MU');  %Up time 1
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'MU')};
        obj_row{row_count,2}={num2str(row_count)};
        
               %% Nov 4 2011 - min. down time
        fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'MD');  %Down time 1
        row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'MD')};
        obj_row{row_count,2}={num2str(row_count)};
        
%         %%
%         %%April 2 2011 - Generator Up and down time
%         % up time S(i,t) >= U(i,t) - U(i,t-1)
%         fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'SG');  %Up time 1
%         row_count=row_count+1;
%         obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'SG')};
%         obj_row{row_count,2}={num2str(row_count)};
%         % Sigma(S) <= u(i,t)
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'SL');  %Up time 2
%         row_count=row_count+1;
%         obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'SL')};
%         obj_row{row_count,2}={num2str(row_count)};
% 
%         % down time H(i,t) >= U(i,t-1) - U(i,t)
%         fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'HG');  %Down time 1
%         row_count=row_count+1;
%         obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'HG')};
%         obj_row{row_count,2}={num2str(row_count)};
%         % Sigma(H) <= 1 - u(i,t)
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'HL');  %Down time 2
%         row_count=row_count+1;
%         obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'HL')};
%         obj_row{row_count,2}={num2str(row_count)};
% 
%         %%        
        
        
        
        num_UC=num_UC+1;
    end
end

obj_row_S_pter=row_count+1;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc segment
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))
        fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Max. cap for reserve 2
                row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
        obj_row{row_count,2}={num2str(row_count)}; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part add new constraints for reserve at each hour
bb={};
for k4=1:a*b
    bb{k4,1}={[cell2mat(W(1)) num2str(fix((k4-1)/b)+1,'%02.0f')]};
end 
for k4=1:a
    for k5=1:b
        bb{(k4-1)*b+k5,1}={[cell2mat(bb{(k4-1)*b+k5,1}) cell2mat(W(2)) num2str(k5,'%02.0f')]};
    end
end

for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'1');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'1')};
        obj_row{row_count,2}={num2str(row_count)}; 
end

for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'2');  % lower bound  >=     spin reserve >= 0.5*percentage of gen.output
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'2')};
        obj_row{row_count,2}={num2str(row_count)}; 
end

for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'3');  % lower bound  >=     non-spin reserve >= 0.5*largest gen. unit
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'3')};
        obj_row{row_count,2}={num2str(row_count)}; 
end
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'4');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'4')};
        obj_row{row_count,2}={num2str(row_count)}; 
end

%% Mar 27 2011 - Regulation up and down
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
        fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'RU');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'RU')};
        obj_row{row_count,2}={num2str(row_count)}; 
end

for k1=1:length(bb(:,1))  % this part add new constraints for reserve             
        fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'RD');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(bb{k1,1}),'RD')};
        obj_row{row_count,2}={num2str(row_count)}; 
end

%% Feb 20 2012 - Restrain down regulation from storage (else extra accumulation)
for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
               fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,2}),'5');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count=row_count+1;
        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,2}),'5')};
        obj_row{row_count,2}={num2str(row_count)};
        
        tmp_name=cell2mat(arcstmp{k1,3});
       len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1                    
%% Feb 20 2012 - down reg restrain 
           fprintf(fid,' L %s\n',strcat(tmp_name,'5'));  
           row_count=row_count+1;
           obj_row{row_count,1}={strcat(tmp_name,'5')};
           obj_row{row_count,2}={num2str(row_count)};
       end
   end
end       

%% Mar 30 2011 - Add compressor reserve and regulation capability

for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
%        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Compressor with regu/res
%        row_count=row_count+1;
%        obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
%        obj_row{row_count,2}={num2str(row_count)};    
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Min. cap for Compressor with regu/res
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
       obj_row{row_count,2}={num2str(row_count)};        
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'a');  % Max. cap for reserve 1
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'a')};
       obj_row{row_count,2}={num2str(row_count)};       
      
       
        %% Compressor Regulation - Mar 30 2011 - Reg up & down
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'c');  %Reg up
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'c')};
       obj_row{row_count,2}={num2str(row_count)}; 
        
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'d');  %Reg down
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'d')};
       obj_row{row_count,2}={num2str(row_count)}; 
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'e');  %RU+SR<10*rr
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'e')};
       obj_row{row_count,2}={num2str(row_count)};  % Feb 18 2012 new
       
        %% CAES Compressor Ramping - Wed, Mar 30 2011 - Ramp up & Ramp down
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'R');  %Ramp up
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'R')};
       obj_row{row_count,2}={num2str(row_count)}; 
        
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'D');  %Ramp down
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'D')};
       obj_row{row_count,2}={num2str(row_count)};     

   end 
end

for k1=1:length(arcstmp(:,1))  % Oct 12 2012 - COmp limit
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Compressor with regu/res
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
       obj_row{row_count,2}={num2str(row_count)};    
   end
end

%%
%% Mar 31 2011 - Add Storage reservoir limits with Creg- and Treg+, Tres+ - ROW

for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Reservoir
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
       obj_row{row_count,2}={num2str(row_count)};    
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Min. cap for Reservoir
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
       obj_row{row_count,2}={num2str(row_count)};     
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'T');  % Mar 18 2012 - Storage Res - min for Turb.DR
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'T')};
       obj_row{row_count,2}={num2str(row_count)};     
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'C');  % Mar 18 2012 - Storage Res - max for Comp RU, SR
       row_count=row_count+1;
       obj_row{row_count,1}={strcat(cell2mat(arcstmp{k1,1}),'C')};
       obj_row{row_count,2}={num2str(row_count)}; 
       
       tmp_name=cell2mat(arcstmp{k1,3});
       len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
%       if(strncmp(tmp_name(5:end),'d02h24',6)==1) %% ROW for final
           fprintf(fid,' L %s\n',strcat(tmp_name(1:4),tmp_name,'L'));  % Max. cap for Reservoir
           row_count=row_count+1;
           obj_row{row_count,1}={strcat(tmp_name(1:4),tmp_name,'L')};
           obj_row{row_count,2}={num2str(row_count)};    
           
           fprintf(fid,' G %s\n',strcat(tmp_name(1:4),tmp_name,'S'));  % Min. cap for Reservoir
           row_count=row_count+1;
           obj_row{row_count,1}={strcat(tmp_name(1:4),tmp_name,'S')};
           obj_row{row_count,2}={num2str(row_count)};         
           
           fprintf(fid,' G %s\n',strcat(tmp_name(1:4),tmp_name,'T'));  % Mar 18 2012 - Storage Res - min for Turb.DR
           row_count=row_count+1;
           obj_row{row_count,1}={strcat(tmp_name(1:4),tmp_name,'T')};
           obj_row{row_count,2}={num2str(row_count)};    
           
           fprintf(fid,' L %s\n',strcat(tmp_name(1:4),tmp_name,'C'));  % Mar 18 2012 - Storage Res - max for Comp RU, SR
           row_count=row_count+1;
           obj_row{row_count,1}={strcat(tmp_name(1:4),tmp_name,'C')};
           obj_row{row_count,2}={num2str(row_count)};             

       end
       
   end 
end
%%

%  Write section COLUMNS
fprintf(fid,'COLUMNS\n');

loop = a*b;
loop_c=0;

for k1=1:length(arcstmp(:,1))
    %%%%%%%%%%%%%% Objective function
    %********************************************************************************************
    % important: When mps file is changed to matrix form, the sequence of variable vector x is decided by the sequence that
    % these variables appear in the objective function.
    % In the objective function part, all transmission lines, existing and petential, are included.
    % Another important rule: The sequence of the equations is decided by
    % the sequence of row name appearing in "section ROWS" of mps file.
    %********************************************************************************************
    %     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)~=1)||(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1)
    fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
    fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,5}))+ str2double(cell2mat(arcstmp{k1,30}))*carbon_tax);
    fprintf(fid,'\n');
    obj_count=obj_count+1;
    obj_variable{obj_count,1}=cell2mat(arcstmp{k1,1});
    obj_variable{obj_count,2}={num2str(obj_count)};
    obj_pointer=obj_pointer+1;
    %%%%%%%%%%%%%% End of Objective function part 1
    k2=1;
    Num=0;
    
    x=cell2mat(arcstmp{k1,2});
    y=cell2mat(arcstmp{k1,3});
    
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
            &&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==0)
        % If it IS a transmission line
        % DC power flow part
        % Not a new line
        % Write the first part of DC power flow: eij
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name= BR+branch name
        fprintf(fid,' %16s                         1',strcat('BR',cell2mat(arcstmp{k1,1})));
        fprintf(fid,'\n');
        
        % Need to update the number of lines before here if I want to update them
        
        % Write the first part of DC power flow: -bij(anglei-anglej)
        % Save it to a matrix and print all of them afterward
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % need to update bij here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bij=str2num(cell2mat(arcstmp{k1,9}))*str2num(cell2mat(arcstmp{k1,11}));
        Num_ang=Num_ang+1;
        Angle_tmp(Num_ang,1)={strcat('AG',cell2mat(arcstmp{k1,2}))};
        Angle_tmp(Num_ang,2)={strcat('BR',cell2mat(arcstmp{k1,1}))};
        Angle_tmp(Num_ang,3)={num2str(-bij)};
        Num_ang=Num_ang+1;
        Angle_tmp(Num_ang,1)={strcat('AG',cell2mat(arcstmp{k1,3}))};
        Angle_tmp(Num_ang,2)={strcat('BR',cell2mat(arcstmp{k1,1}))};
        Angle_tmp(Num_ang,3)={num2str(bij)};
        
        % The first part of KCL part of DC power flow is written in the transportation part.
        % Write the second part of KCL part of DC power flow.
        % Here only existing lines considered
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
        fprintf(fid,' %16s                        -1',cell2mat(arcstmp{k1,2}));
        fprintf(fid,' %16s                         1',cell2mat(arcstmp{k1,3}));
        fprintf(fid,'\n');
        
    else
        if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'DD',2)==1)
            % Not a transmission line and  a demand line
            % for demand line, there will be no KCL equation for demand node
            fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
            fprintf(fid,' %16s                        -1',cell2mat(arcstmp{k1,2}));
            fprintf(fid,'\n');
        else
            if strncmp(cell2mat(arcstmp{k1,2}),'X',1)==0
                %   If it is not coming from a dummy node(Fuel production node)
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
                if strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1 % storage efficiency MCS
                    fprintf(fid,' %25s',(ke+str2double(cell2mat(arcstmp{k1,6}))));
                else
                    fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
                end
                fprintf(fid,' %16s                        -1',cell2mat(arcstmp{k1,2}));
                fprintf(fid,'\n');
            else  % From a dummy node
                if strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==0
                    fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
                    fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
                    fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
                    fprintf(fid,'\n');
                end;
            end;
        end
    end;
    
%% Mar 30 2011 - Compressor max and min equations with Res/Reg
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        
        % Compr + reg down < max.capacity
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'L'));
        fprintf(fid,'\n');       
        
        % Compr - reserve 1 - reg up > min.capacity
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,'\n'); 
        
        %% Ramp up % down - Mar 30 2011
        loop_c=loop_c+1;
        if (loop_c~=1 && loop_c~=loop)
            % Ramp up --> comp(t) - comp(t-1) <= rup;
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,' %16s                        -1',strcat(cell2mat(arcstmp{k1+1,1}),'R'));
            fprintf(fid,'\n');
            
            % Ramp down --> comp(t-1) - comp(t) <= rdown;
            
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
            fprintf(fid,' %25s',num2str(-1));
            fprintf(fid,' %16s                        1',strcat(cell2mat(arcstmp{k1+1,1}),'D'));
            fprintf(fid,'\n');
        else
            if(loop_c==1)
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R - up
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %comment
                fprintf(fid,' %25s',num2str(1));%comment
                fprintf(fid,' %16s                        -1',strcat(cell2mat(arcstmp{k1+1,1}),'R'));
                fprintf(fid,'\n');
                
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R- down
                fprintf(fid,' %16s                        1',strcat(cell2mat(arcstmp{k1+1,1}),'D'));
                fprintf(fid,'\n');
            else
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R - up
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
                fprintf(fid,' %25s',num2str(1));
                fprintf(fid,'\n');
                
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R- down
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
                fprintf(fid,' %25s',num2str(-1));
                fprintf(fid,'\n');

                loop_c=0;
            end
        end
               
    end
%%       

%% Mar 31 2011 - Reservoir max and min equations with Res/Reg
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        
        % Reservoir + reg down < max.capacity
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'L'));
        fprintf(fid,'\n');       
        
        % Reservoir - reg up - Sp. res > min.capacity
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,'\n'); 
        
              % Mar 18 2012 - Res(t-1) - tur DR> min
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'T'));
        fprintf(fid,'\n'); 
        
                      % Mar 18 2012 - Res(t-1) + com ur + com SR< max
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'C'));
        fprintf(fid,'\n'); 
    
    end
%%    
         
    % The previous part is similar as the Economic Dispatch part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The following part deals with the reserve equations that arcsname(energy bidding) participates
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        
        % Eng.bid1 + Eng.bid2 + Eng.bid3 - arcname=0
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
        fprintf(fid,' %25f',-1);
        fprintf(fid,'\n');
        
        % Energy bid + reserve 1  -max.capacity*U <0
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'L')); % Energy bid + reserve 1 -max.capacity*U <0
        fprintf(fid,'\n');
        % Energy bid + reserve 1  -min.capacity*U >0
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'S')); % Energy bid + reserve 1 -max.capacity*U <0
        fprintf(fid,'\n');                
        % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if  if reserve.indicator1 ==1)
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'M')); % Energy bid + reserve 1 -max.capacity*U <0
        fprintf(fid,'\n');
        
        % this part deal with the percentage reserve constraint (the second)
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25.4f',-0.5*Pctg_nowind);
        fprintf(fid,'\n');
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
        fprintf(fid,' %25.4f',-Pctg_nowind);% Feb 18 2012 new
        fprintf(fid,'\n');
        
        %% Ramp up % down - Mar 11 2011
        loop_c=loop_c+1;
        if (loop_c~=1 && loop_c~=loop)
            % Ramp up --> gen(t) - gen(t-1) <= rup;
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,' %16s                        -1',strcat(cell2mat(arcstmp{k1+1,1}),'R'));
            fprintf(fid,'\n');
            
            % Ramp down --> gen(t-1) - gen(t) <= rdown;
            
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
            fprintf(fid,' %25s',num2str(-1));
            fprintf(fid,' %16s                        1',strcat(cell2mat(arcstmp{k1+1,1}),'D'));
            fprintf(fid,'\n');
        else
            if(loop_c==1)
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R - up
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %comment
                fprintf(fid,' %25s',num2str(1));%comment
                fprintf(fid,' %16s                        -1',strcat(cell2mat(arcstmp{k1+1,1}),'R'));
                fprintf(fid,'\n');
                
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R- down
                fprintf(fid,' %16s                        1',strcat(cell2mat(arcstmp{k1+1,1}),'D'));
                fprintf(fid,'\n');
            else
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R - up
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
                fprintf(fid,' %25s',num2str(1));
                fprintf(fid,'\n');
                
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}))); % R- down
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
                fprintf(fid,' %25s',num2str(-1));
                fprintf(fid,'\n');

                loop_c=0;
            end
        end        
    end
    
    
    % this part deal with the percentage reserve constraint (the second)
    if strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1          % the percentage of wind generator output that needs to be backed up by spin reserve
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25.4f',-0.5*Pctg_wind);
        fprintf(fid,'\n');
    end
    if strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1        % the percentage of wind generator output that needs to be backed up by non-spin reserve
        fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));        % print the row name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
%        fprintf(fid,' %25.4f',-0*Pctg_wind);% Mar 25 2011
        fprintf(fid,' %25.4f',-Pctg_wind);% Feb 18 2012 new
        fprintf(fid,'\n');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following part define the energy bidding variables (3 for each gen.)
flag_enebid=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,20})));
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('1',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_enebid==0
            flag_enebid=1;
            obj_enebid=obj_count;
        end
%         % Energy bid + reserve 1 + reserve 2 -max.capacity*U <0
%         fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        % Energy bid + reserve 1 + reserve 2 -min.capacity*U >0
%         fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');

        fprintf(fid,' %16s',strcat('2',cell2mat(arcstmp{k1,1}))); % Second energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,23})));
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('2',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
%         % Energy bid + reserve 1 + reserve 2 -max.capacity*U <0
%         fprintf(fid,' %16s',strcat('2',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         % Energy bid + reserve 1 + reserve 2 -min.capacity*U >0
%         fprintf(fid,' %16s',strcat('2',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        fprintf(fid,' %16s',strcat('2',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('3',cell2mat(arcstmp{k1,1}))); % Third energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,26})));
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('3',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
%         % Energy bid + reserve 1 + reserve 2 -max.capacity*U <0
%         fprintf(fid,' %16s',strcat('3',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         % Energy bid + reserve 1 + reserve 2 -min.capacity*U >0
%         fprintf(fid,' %16s',strcat('3',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        fprintf(fid,' %16s',strcat('3',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following part define the reserve 1(spin reserve) and reserve 2(non-spin reserve) varibales
flag_reserve1=0;
flag_reserve2=0;
loop = a*b;
loop_c=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function      
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        %% Apr 05 2011 - hourly spin bids
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
        %%
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('A',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_A==0
            obj_A=obj_count;
            flag_obj_A=1;
        end
%         % Energy bid + reserve 1  -max.capacity*U <0 %% Oct 07 2011 Modified
%         fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        
        % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if  if reserve.indicator1 ==1)
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));      % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'M')); % Energy bid + reserve 1 -max.capacity*U <0
        fprintf(fid,'\n');
        
        % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'1'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        % Sum of Reserve biddings > OR2, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
                % Sum of Reserve biddings > OR2, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3')); %Feb 18 2012 new
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
               % Sum of Reserve biddings > OR2, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4')); %Feb 18 2012 new
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');

        
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'e'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new
        
%         %% Ramp up - Mar 11 2011 in Spin reserve
%         % Ramp up --> gen(t) - gen(t-1) + Res1(t) + Res2(t) <= rup;
%         
%         fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %includes hour=1
%         fprintf(fid,' %25s',num2str(1));
%         fprintf(fid,'\n');
      
%% Feb 15 2012
% %% Mar 31 2011 Turbine Sp. res inside reservoir        
%         if (strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1)
%             % Reservoir - Reg up - Sp. Res > min.
%             resr_name = cell2mat(arcstmp{k1,2});
%             fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); 
%             fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'S'));
%             fprintf(fid,' %25f',-1);
%             fprintf(fid,'\n');
%         end
        
        %% Feb 07 2012- Account for ancillary in reservoir dynamics
        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,2}));
                fprintf(fid,' %25f',-1*kf);
                fprintf(fid,'\n');
%           end 
        end
        
        if(loop_c == loop)
            loop_c=0;
        end
            
%%            
%         if (loop_c~=1 && loop_c~=loop)
%             fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
%             fprintf(fid,' %25s',num2str(1));
%             fprintf(fid,'\n');            
%         else
%             if(loop_c~=1)
%                 fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); % R - up
%                 fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
%                 fprintf(fid,' %25s',num2str(1));
%                 fprintf(fid,'\n');               
%                 loop_c=0;
% %             else % comment
% %                 fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); % R - up
% %                 fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
% %                 fprintf(fid,' %25s',num2str(1));
% %                 fprintf(fid,'\n');
%             end
%         end
    end
end

% loop=a*b;
% loop_c=0;
% 
% for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         tmp_name=cell2mat(arcstmp{k1,2});
%         % Claim this variable in the objective function
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));
%                 %% Apr 05 2011 - Non-Spin bid hourly
%         loop_c=loop_c+1;
%         fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
%         %%
% 
%         fprintf(fid,'\n');
%         obj_count=obj_count+1;
%         obj_variable{obj_count,1}=strcat('B',cell2mat(arcstmp{k1,1}));
%         obj_variable{obj_count,2}={num2str(obj_count)};
%            if flag_obj_B==0
%                obj_B=obj_count;
%                flag_obj_B=1;
%            end
%            
% %         % Energy bid + Nsp -min.capacity*U >0
% %         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
% %         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
% %         fprintf(fid,' %25f',1);
% %         fprintf(fid,'\n');
%         
%         % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if  if reserve.indicator1 ==1)
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));      % print the row name
%         fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'M')); % Energy bid + reserve 1 -max.capacity*U <0
%         fprintf(fid,'\n');
%         
%         % Sum of Reserve biddings > OR1, reserve 2 bidding
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
%         time_tmp=cell2mat(arcstmp{k1,2});
%         fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         
%         % Sum of Reserve biddings > OR2, reserve 2 bidding
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
%         fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         
% %% mar 11 2011 - Ramp up for NSR
% %         %loop_c=loop_c+1;
% %         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));
% %         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %includes hour=1
% %         fprintf(fid,' %25s',num2str(1));
% %         fprintf(fid,'\n');
% 
%         %% Feb 07 2012- Account for ancillary in reservoir dynamics
%         if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
% %            if loop_c<loop
%                 fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1}))); 
%                 fprintf(fid,' %16s',cell2mat(arcstmp{k1,2}));
%                 fprintf(fid,' %25f',-1);
%                 fprintf(fid,'\n');
% %           end 
%         end
% 
%         if(loop_c == loop)
%             loop_c=0;
%         end
% 
%         %% Feb 15 2012
% %         %% Oct 03 2011 Turbine Non Sp.res inside reservoir
% %         if (strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1)
% %             % Reservoir - Reg up - Sp. Res - Non Sp.Res > min.
% %             resr_name = cell2mat(arcstmp{k1,2});
% %             fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1}))); 
% %             fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'S'));
% %             fprintf(fid,' %25f',-1);
% %             fprintf(fid,'\n');
% %         end
%         
% %         if (loop_c~=1 && loop_c~=loop)
% %             fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));
% %             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
% %             fprintf(fid,' %25s',num2str(1));
% %             fprintf(fid,'\n');            
% %         else
% %             if(loop_c~=1)
% %                 fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1}))); % R - up
% %                 fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
% %                 fprintf(fid,' %25s',num2str(1));
% %                 fprintf(fid,'\n');               
% %                 loop_c=0;
% % %             else  %comment
% % %                 fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1}))); % R - up
% % %                 fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R'));
% % %                 fprintf(fid,' %25s',num2str(1));
% % %                 fprintf(fid,'\n');
% %             end
% %         end
%     end
% end


%% Mar 27 2011 Regulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following part define the Regulation variables
flag_obj_RU = 0;
obj_RU=0;
loop = a*b;
loop_c=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));
        %% Apr 05 2011 - Reg bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,31})));%reg_bid(loop_c));
        %%
%        fprintf(fid,'              obj %25.4f',30);
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('RU',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_RU ==0
            obj_RU =obj_count;
            flag_obj_RU =1;
        end
%         % Energy bid + reserve 1  -max.capacity*U <0
%         fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
        
        % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if  if reserve.indicator1 ==1)
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));      % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'M')); % Energy bid + reserve 1 -max.capacity*U <0
        fprintf(fid,'\n');
        
        % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'RU'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');        
        
                % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'1'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');      %Feb 18 2012 new  
        
                % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');        %Feb 18 2012 new
        
                        % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');      %Feb 18 2012 new  
        
                % Sum of Reserve biddings > OR1, reserve 1 bidding
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');        %Feb 18 2012 new

        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'c'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % ru+sr+nsr<10*rr
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'e'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new        
        
        %% Ramp up - Mar 27 2011 in Regulation
        % Ramp up --> gen(t) - gen(t-1) + Res1(t) + Res2(t) <= rup;

%         fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %includes hour=1
%         fprintf(fid,' %25s',num2str(1));
%         fprintf(fid,'\n');

        %% Feb 07 2012- Account for ancillary in reservoir dynamics
        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,2}));
                fprintf(fid,' %25f',-1*kf);
                fprintf(fid,'\n');
%            end
        end
        
        if(loop_c == loop)
            loop_c=0;
        end
        
        %% Feb 15 2012
%         %% Mar 31 2011 Turbine Up regulation inside reservoir        
%         if (strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1)
%             % Reservoir - Reg up - Sp. Res > min.
%             resr_name = cell2mat(arcstmp{k1,2});
%             fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1}))); 
%             fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'S'));
%             fprintf(fid,' %25f',-1);
%             fprintf(fid,'\n');
%         end
        %%
    end
end

%% Mar 27 2011 Down regulation
flag_obj_RD = 0;
obj_RD=0;

loop=a*b;
loop_c=0;

for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));
                %% Apr 05 2011 - Reg bid hourly
        loop_c=loop_c+1;
        %% Feb 22 2012 - Down regulation from turbine (only if greater than bid)
%         if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%             fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,31}))+str2double(cell2mat(arcstmp{k1,5})));%reg_bid(loop_c));
%         else
            fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,31})));%reg_bid(loop_c));
%         end
        %%
%        fprintf(fid,'              obj %25.4f',30);
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('RD',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_RD==0
            obj_RD=obj_count;
            flag_obj_RD=1;
        end
        % Energy bid -min.capacity*U >0
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,' %25f',-1);
        fprintf(fid,'\n');
        
        % Sum of Reserve biddings > OR1, reserve 2 bidding
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'RD'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'d'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
%% mar 27 2011 - Regulation down 
 
%         if(loop_c~=1)
%             fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D')); %includes hour=1
%             fprintf(fid,' %25s',num2str(1));
%             fprintf(fid,'\n');
%         end
        
        %% Feb 07 2012- Account for ancillary in reservoir dynamics
        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
 %           if loop_c<loop
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,2}));
                fprintf(fid,' %25f',1*kf);
                fprintf(fid,'\n');
                %% Feb 20 2012- restrain down reg in stor
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,2}),'5'));
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');                
  %          end
        end
        
                %% Mar 18 2012 Turbine Down regulation inside reservoir (t-1)       
        if (strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1)
            % Reservoir(t-1) - Reg down(t) > min.
            %if(loop_c>1)
            resr_name = cell2mat(arcstmp{k1,2});
            fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1}))); 
            fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'T'));
            fprintf(fid,' %25f',-1);
            fprintf(fid,'\n');
            %end
        end

        if(loop_c == loop)
            loop_c=0;
        end

    end
end

%%

flag_obj_comA = 0;
obj_comA = 0;
num_comA = 0;

loop = a*b;
loop_c = 0;
%% Mar 30 2011 - Compressor min equation with Spinning res
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        %% Apr 05 2011 - Hourly spin bids
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
%%
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('A',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        num_comA = num_comA+1;
        if flag_obj_comA==0
            obj_comA=obj_count;
            flag_obj_comA=1;
        end
        
        % Compr - up Reg - Res > min.
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,' %25f',-1);
        fprintf(fid,'\n');
     
        % Sum of Compr Reserve biddings > OR1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'1'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        % Sum of Compr Reserve biddings > OR2, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
                % Sum of Compr Reserve biddings > OR1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');% Feb 18 2012 new
        
        % Sum of Compr Reserve biddings > OR2, reserve 1 bidding
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new
        
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'e'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new        
        
%% Mar 30 2011 - Compressor Reserve - ramp down rate 
 
%         if(loop_c~=1)
%             fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
%             fprintf(fid,' %25s',num2str(1));
%             fprintf(fid,'\n');
%         end

        %% Feb 07 2012- Account for ancillary in reservoir dynamics
%        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
                fprintf(fid,' %25f',(-kf*(ke+str2double(cell2mat(arcstmp{k1,6})))));
                fprintf(fid,'\n');
%           end 
%        end
        
                %% Mar 18 2012 Turbine Down regulation inside reservoir (t-1)       
             % Reservoir(t-1) + Com Reg up(t)+ com SR (t) <max
            %if(loop_c>1)
                resr_name = cell2mat(arcstmp{k1,3});
                fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'C'));
                fprintf(fid,' %25s',(ke+str2double(cell2mat(arcstmp{k1,6}))));
                fprintf(fid,'\n');
            %end        
        
        if(loop_c == loop)
            loop_c=0;
        end
    end
end
%%             

flag_obj_comRU = 0;
obj_comRU = 0;
num_comRU = 0;

loop = a*b;
loop_c=0;

%% Mar 30 2011 - Compressor min equation with Up regu
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));
        %% Apr 05 2011 - Reg bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,31})));%reg_bid(loop_c));
        %%
 %       fprintf(fid,'              obj %25.4f',30);
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('RU',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        num_comRU = num_comRU+1;
        if flag_obj_comRU==0
            obj_comRU=obj_count;
            flag_obj_comRU=1;
        end
        
        % Compr - up Reg - Res > min.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,' %25f',-1);
        fprintf(fid,'\n');
     
        % Sum of Up Regulation biddings > RU1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'RU'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
                % Sum of Up Regulation biddings > RU1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'1'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');% Feb 18 2012 new
        
                % Sum of Up Regulation biddings > RU1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'2'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new
        
                        % Sum of Up Regulation biddings > RU1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');% Feb 18 2012 new
        
                % Sum of Up Regulation biddings > RU1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new
        
      
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'c'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'e'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new
        
%% Mar 30 2011 - Compressor Up Reg - ramp down rate 
       
%         if(loop_c~=1)
%             fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'D'));
%             fprintf(fid,' %25s',num2str(1));
%             fprintf(fid,'\n');
%         end
        

        %% Feb 07 2012- Account for ancillary in reservoir dynamics
       % if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
                fprintf(fid,' %25f',(-kf * (ke+str2double(cell2mat(arcstmp{k1,6})))));
                fprintf(fid,'\n');
%            end
       % end
        
                       %% Mar 18 2012 Turbine Down regulation inside reservoir (t-1)       
             % Reservoir(t-1) + Com Reg up(t)+ com SR (t) <max
            %if(loop_c>1)
                resr_name = cell2mat(arcstmp{k1,3});
                fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'C'));
                fprintf(fid,' %25s',(ke+str2double(cell2mat(arcstmp{k1,6}))));
                fprintf(fid,'\n');
            %end
            
        if(loop_c == loop)
            loop_c=0;
        end
        
    end
end
%%             

flag_obj_comRD = 0;
obj_comRD = 0;
num_comRD = 0;

loop = a*b;
loop_c=0;

%% Mar 30 2011 - Compressor max equation with Down regu
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));
                %% Apr 05 2011 - Reg bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,31})));%/str2double(cell2mat(arcstmp{k1,6})));%reg_bid(loop_c));% Feb 18 2012 new
        %%
%        fprintf(fid,'              obj %25.4f',30);
        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('RD',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        num_comRD = num_comRD+1;
        if flag_obj_comRD==0
            obj_comRD=obj_count;
            flag_obj_comRD=1;
        end
        
        % Compr + down Reg < max.
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
     
        % Sum of Down Regulation biddings > OR1, reserve 1 bidding from compr.
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % print the reg 2 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'RD'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
               
        fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'d'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
       
        %% Mar 30 2011 - Compressor Down Reg - ramp up rate 
        % Ramp up --> gen(t) - gen(t-1) + Res1(t) + Res2(t) <= rup;
 
%         fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'R')); %includes hour=1
%         fprintf(fid,' %25s',num2str(1));
%         fprintf(fid,'\n');

  %% Feb 07 2012- Account for ancillary in reservoir dynamics
       % if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
                fprintf(fid,' %25f',kf * (ke+str2double(cell2mat(arcstmp{k1,6}))));
                fprintf(fid,'\n');                
              
                %%Feb 20 2012 Restrain down reg in stor
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,3}),'5'));
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');                            
                
%            end
       % end
        if(loop_c == loop)
            loop_c=0;
        end
        %% Feb 15 2012
%         %% Mar 31 2011 Compressor Down regulation in reservoir
%         resr_name=cell2mat(arcstmp{k1,3});
%         fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
%         fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'L'));
%         fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,6})));
%         fprintf(fid,'\n');
    end
end
%%             

%% Mar 31 2011 Storage reservoir final - STxx STxx d<final> h<final> 
obj_final_st = 0;
flag_final_st = 0;
final_c = 0;
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_name=cell2mat(arcstmp{k1,3});
        len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
%        if(strncmp(tmp_name(5:end),'d02h24',6)==1)
            final_c = final_c+1;
            % Claim this variable in the objective function
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,5})));
            fprintf(fid,'\n');
            obj_count=obj_count+1;
            obj_variable{obj_count,1}=strcat(tmp_name(1:4),tmp_name);
            obj_variable{obj_count,2}={num2str(obj_count)};
            if flag_final_st==0
                obj_final_st=obj_count;
                flag_final_st=1;
            end
            
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,' %16s',cell2mat(arcstmp{k1,3}));
            fprintf(fid,' %25s',num2str(-1));
            fprintf(fid,'\n');
            
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name,'L'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,'\n');
            
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name,'S'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,'\n');
            
                         % Mar 18 2012 - Res(t-1) - tur DR> min
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name,'T'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,'\n');
            
                      % Mar 18 2012 - Res(t-1) + com ur + com SR< max
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name,'C'));
            fprintf(fid,' %25s',num2str(1));
            fprintf(fid,'\n');
            
       
            %% start here - add this variable to ST15d02h24 ROW , BND, and
            %% ST15ST15d02h24L and S            
        end
    end
end

%%

loop=a*b;
loop_c=0;
flag_obj_Bn=0;
obj_Bn=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));
                %% Apr 05 2011 - Non-Spin bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,28})));%nonspin_bid(loop_c));
        %%

        fprintf(fid,'\n');
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('Bn',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_Bn==0
            obj_Bn=obj_count;
            flag_obj_Bn=1;
        end
        
        % Feb 18 2012 new
        % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if  if reserve.indicator1 ==1)
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));      % print the row name
        fprintf(fid,' %16s                         1',strcat(cell2mat(arcstmp{k1,1}),'M')); % Energy bid + reserve 1 -max.capacity <0
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'b'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'e'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n'); % Feb 18 2012 new        
        
        % Sum of Reserve biddings > OR1, reserve 2 bidding
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        time_tmp=cell2mat(arcstmp{k1,2});
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'3'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
        % Sum of Reserve biddings > OR2, reserve 2 bidding
        fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(time_tmp(length(time_tmp)-6+1:length(time_tmp)),'4'));
        fprintf(fid,' %25f',1);
        fprintf(fid,'\n');
        
  %% Feb 07 2012- Account for ancillary in reservoir dynamics
        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
%            if loop_c<loop
                fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',cell2mat(arcstmp{k1,2}));
                fprintf(fid,' %25f',-1*kf);
                fprintf(fid,'\n');
%            end
        end
       
       
        if(loop_c == loop)
            loop_c=0;
        end
%% Feb 15 2012
%         %% Oct 03 2011 Turbine Non Sp.res inside reservoir
%         if (strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1)
%             % Reservoir - Reg up - Sp. Res - Non Sp.Res_ON - Non Sp.Res_OFF > min.
%             resr_name = cell2mat(arcstmp{k1,2});
%             fprintf(fid,' %16s',strcat('Bn',cell2mat(arcstmp{k1,1}))); 
%             fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'S'));
%             fprintf(fid,' %25f',-1);
%             fprintf(fid,'\n');
%         end      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following part define the unit commitment variables (for each gen. except for GenWind)
flag_UC=0;
loop=a*b;
loop_c=0;
count=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        % Claim the UC variable in the objective function
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        len_arc=length(cell2mat(arcstmp{k1,1}));
        arc_tmp=cell2mat(arcstmp{k1,1});
        fprintf(fid,'              obj %25.4f',0);
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('U',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        fprintf(fid,'\n');
        if flag_UC==0 
            obj_UC=obj_count;
            flag_UC=1;
        end
        
        % This part deal with the equation U(t)-U(t-1)=X-Y startup and
        % shutdown cost %% Mar 31 2011 correction
        % 1. U(t)-U(t+1)=X-Y, 2. d01h24 and d02h24 wrong, only X-Y
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));
        fprintf(fid,'%25.4f',1);
        fprintf(fid,'\n');
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)~=1
            fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1+1,1}),'U'));
            fprintf(fid,'%25.4f',-1);
            fprintf(fid,'\n');
        end
                 
        % Energy bid - max.capacity <0 %% Oct 07 2011 - Modified
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,8})));
        fprintf(fid,'\n');
        
        % Energy bid  - min.capacity >0
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,7})));
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,14}))/6); %% Nov 05
        fprintf(fid,'\n');

        %% Mar 27 Regulation            
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'c')); %Up regu
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,14}))/12);
        fprintf(fid,'\n');
            
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'d')); % Down regu
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,15}))/12);
        fprintf(fid,'\n');
        
        
        %% Oct 03 2011 Non-Spinning
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'N'));
        fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,13})));
        fprintf(fid,'\n');

%% Nov 04 2011 - Min Up/down time
        loop_c= loop_c + 1;
        count = count+1;
        
        %% Nov 04 2011 - Min Up time
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MU'));
        fprintf(fid,'%25.4f',1);
        fprintf(fid,'\n');
                  
%        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',1),'h',sprintf('%02.0f',1)),6)~=1
        for mu_c=1:UT(count,1)-1
            if (loop_c-mu_c)>0
                fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1-mu_c,1}),'MU'));
                fprintf(fid,'%25.4f',1);
                fprintf(fid,'\n');
            end
        end
     %% Nov 04 2011 - Min down time 
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MD'));
        fprintf(fid,'%25.4f',-1);
        fprintf(fid,'\n');
                  
%        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',1),'h',sprintf('%02.0f',1)),6)~=1
        for md_c=1:DT(count,1)-1
            if (loop_c-md_c)>0
                fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1-md_c,1}),'MD'));
                fprintf(fid,'%25.4f',-1);
                fprintf(fid,'\n');
            end
        end
            
        if loop_c==loop
            loop_c=0;
        end
        
    %% Oct 12 2012 - Disjoint operation
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
        % Comp bid  + reg- < max.capacity (1-U-Uo)
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(arc_tmp(5:8),arc_tmp(1:4),arc_tmp(9:len_arc),'L'));
        fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,8})));
        fprintf(fid,'\n');
        
        % Comp bid  - sr - reg+ >min.capacity (1-U-Uo)
        fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(arc_tmp(5:8),arc_tmp(1:4),arc_tmp(9:len_arc),'S'));
        fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,7})));
        fprintf(fid,'\n');
    end
      %% Oct 12 2012 - Disjoint operation  
        
% %%        % Min. Up time - April 2 2011
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'SG'));
%         fprintf(fid,'%25.4f',-1);
%         fprintf(fid,'\n');
%         if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)~=1
%             fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1+1,1}),'SG'));
%             fprintf(fid,'%25.4f',1);
%             fprintf(fid,'\n');
%         end
%         
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'SL'));
%         fprintf(fid,'%25.4f',-1);
%         fprintf(fid,'\n');
%         
% %%        % Min. Down time - April 2 2011
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'HG'));
%         fprintf(fid,'%25.4f',1);
%         fprintf(fid,'\n');
%         if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)~=1
%             fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1+1,1}),'HG'));
%             fprintf(fid,'%25.4f',-1);
%             fprintf(fid,'\n');
%         end
%         
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'HL')); %% Needs RHS =1
%         fprintf(fid,'%25.4f',1);
%         fprintf(fid,'\n');
        
    end
  
end
obj_UC_num=obj_count-obj_UC+1;
loop=a*b;
loop_c=0;
count=0;
flag_UCNsp=0; %% Oct 03 - Non-Spinning
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        % Claim the UC variable in the objective function
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        len_arc=length(cell2mat(arcstmp{k1,1})); 
        arc_tmp=cell2mat(arcstmp{k1,1});
        fprintf(fid,'              obj %25.4f',0);
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('N',cell2mat(arcstmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        fprintf(fid,'\n');
        if flag_UCNsp==0 
            obj_UCNsp=obj_count;
            flag_UCNsp=1;
        end
        
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'b'));
        fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,14}))/6);
        fprintf(fid,'\n');

        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'N'));
        fprintf(fid,'%25.4f',1);%% Nov 05
        fprintf(fid,'\n');
        
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'Un'));
        fprintf(fid,'%25.4f',1);
        fprintf(fid,'\n');
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)~=1
            fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
            fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1+1,1}),'Un'));
            fprintf(fid,'%25.4f',-1);
            fprintf(fid,'\n');
        end
        
        %% Nov 04 2011 - Min Up time
        loop_c= loop_c + 1;
        count = count+1;
        
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MU'));
        fprintf(fid,'%25.4f',1);
        fprintf(fid,'\n');
%        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',1),'h',sprintf('%02.0f',1)),6)~=1
        for mu_c=1:UT(count,1)-1
            if (loop_c-mu_c)>0
                fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1-mu_c,1}),'MU'));
                fprintf(fid,'%25.4f',1);
                fprintf(fid,'\n');
            end
        end
        
        %% Min. down time
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MD'));
        fprintf(fid,'%25.4f',-1);
        fprintf(fid,'\n');
%        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',1),'h',sprintf('%02.0f',1)),6)~=1
        for md_c=1:DT(count,1)-1
            if (loop_c-md_c)>0
                fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1-md_c,1}),'MD'));
                fprintf(fid,'%25.4f',-1);
                fprintf(fid,'\n');
            end
        end
            
        if loop_c==loop
            loop_c=0;
        end
        
            %% Oct 12 2012 - Disjoint operation
            if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
        % Comp bid  + reg- < max.capacity (1-U-Uo)
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(arc_tmp(5:8),arc_tmp(1:4),arc_tmp(9:len_arc),'L'));
        fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,8})));
        fprintf(fid,'\n');
        
        % Comp bid  - sr - reg+ >min.capacity (1-U-Uo)
        fprintf(fid,' %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %16s',strcat(arc_tmp(5:8),arc_tmp(1:4),arc_tmp(9:len_arc),'S'));
        fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,7})));
        fprintf(fid,'\n');
            end
      %% Oct 12 2012 - Disjoint operation 
        
    end
end
obj_UCNsp_num=obj_count-obj_UCNsp+1;
loop=a*b;
loop_c=0;
count=0;
for k1=1:length(arcstmp(:,1))   % This part define the startup and shutdown cost
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        % First define the variable in the objective function 
        fprintf(fid,' %16s',strcat('X',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,16})));
        fprintf(fid,'\n');   
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('X',cell2mat(arcstmp{k1,1})); % Startup cost
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_X==0
            obj_X=obj_count;
            flag_obj_X=1;
        end
        fprintf(fid,' %16s',strcat('X',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));
        fprintf(fid,' %25.4f',-1);
        fprintf(fid,'\n');    
        
        %% Nov 4 2011 Min up time
        loop_c=loop_c+1;
        count=count+1;
        fprintf(fid,' %16s',strcat('X',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MU'));
        fprintf(fid,' %25.4f',-min(UT(count,1),loop+1-loop_c));
        fprintf(fid,'\n'); 
        
        if loop_c==loop
            loop_c=0;
        end
        
    end
end

loop = a*b;
loop_c = 0;
count=0;
for k1=1:length(arcstmp(:,1))  % This part define the startup and shutdown cost
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        fprintf(fid,' %16s',strcat('Y',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,17})));
        fprintf(fid,'\n'); 
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('Y',cell2mat(arcstmp{k1,1})); % Shutdown cost
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_Y==0
            obj_Y=obj_count;
            flag_obj_Y=1;
        end
        fprintf(fid,' %16s',strcat('Y',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));
        fprintf(fid,' %25.4f',1);
        fprintf(fid,'\n'); 
        
        
        %% Nov 4 2011 Min down time
        loop_c=loop_c+1;
        count=count+1;
        fprintf(fid,' %16s',strcat('Y',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MD'));
        fprintf(fid,' %25.4f',-min(DT(count,1),loop+1-loop_c));
        fprintf(fid,'\n'); 
        
        if loop_c==loop
            loop_c=0;
        end
    end
end

%% Oct 20 2011 (NSp start up and shut down)
flag_obj_Xn=0;
obj_Xn=0;
flag_obj_Yn=0;
obj_Yn=0;
loop=a*b;
loop_c=0;
count=0;
for k1=1:length(arcstmp(:,1))   % This part define the startup and shutdown cost
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        % First define the variable in the objective function 
        fprintf(fid,' %16s',strcat('Xn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,16})));
        fprintf(fid,'\n');   
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('Xn',cell2mat(arcstmp{k1,1})); % Startup cost
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_Xn==0
            obj_Xn=obj_count;
            flag_obj_Xn=1;
        end
        fprintf(fid,' %16s',strcat('Xn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'Un'));
        fprintf(fid,' %25.4f',-1);
        fprintf(fid,'\n');      
        
                %% Nov 4 2011 Min up time
        loop_c=loop_c+1;
        count=count+1;
        fprintf(fid,' %16s',strcat('Xn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MU'));
        fprintf(fid,' %25.4f',-min(UT(count,1),loop+1-loop_c));
        fprintf(fid,'\n'); 
        
        if loop_c==loop
            loop_c=0;
        end
    end
end

loop=a*b;
loop_c=0;
count=0;
for k1=1:length(arcstmp(:,1))  % This part define the startup and shutdown cost
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        fprintf(fid,' %16s',strcat('Yn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,17})));
        fprintf(fid,'\n'); 
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('Yn',cell2mat(arcstmp{k1,1})); % Shutdown cost
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_Yn==0
            obj_Yn=obj_count;
            flag_obj_Yn=1;
        end
        fprintf(fid,' %16s',strcat('Yn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'Un'));
        fprintf(fid,' %25.4f',1);
        fprintf(fid,'\n'); 
        
        %% Nov 4 2011 Min down time
        loop_c=loop_c+1;
        count=count+1;
        fprintf(fid,' %16s',strcat('Yn',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'MD'));
        fprintf(fid,' %25.4f',-min(DT(count,1),loop+1-loop_c));
        fprintf(fid,'\n'); 
        
        if loop_c==loop
            loop_c=0;
        end
    end
end


% %% Apr 2 2011 S variable definitions - Generator Up time
% 
% flag_S_up = 0;
% num_S_up = 0;
% obj_S_up = 0;
% UT = 1;
% % load UpDn
% for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         % Claim the S variable in the objective function
%         fprintf(fid,' %16s',strcat('S',cell2mat(arcstmp{k1,1})));
%         len_arc=length(cell2mat(arcstmp{k1,1})); 
%         arc_tmp=cell2mat(arcstmp{k1,1});
%         fprintf(fid,'              obj %25.4f',0);
%         obj_count=obj_count+1;
%         obj_variable{obj_count,1}=strcat('S',cell2mat(arcstmp{k1,1}));
%         obj_variable{obj_count,2}={num2str(obj_count)};
%         fprintf(fid,'\n');
%         num_S_up =num_S_up +1;
%         if flag_S_up==0 
%             obj_S_up=obj_count;
%             flag_S_up=1;
%         end
%         
%         % S in SG - up time
%         fprintf(fid,' %16s',strcat('S',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'SG'));
%         fprintf(fid,'%25.4f',1);
%         fprintf(fid,'\n');
% 
%         % S in SL - sigma (S)- up time
%         for c_ut = 1:UT%(num_S_up,1)
%             arc_tmpi = cell2mat(arcstmp{k1+c_ut-1,1});
%             len_arci = length(arc_tmpi);
%             fprintf(fid,' %16s',strcat('S',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(arc_tmpi,'SL'));
%             fprintf(fid,'%25.4f',1);
%             fprintf(fid,'\n'); 
%             if strncmp(arc_tmpi(len_arci-5:len_arci),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
%                 break;
%             end
%         end               
%     end
% end
% 
% %% Apr 2 2011 H variable definitions - Generator Down time
% 
% flag_H_down = 0;
% num_H_down = 0;
% obj_H_down = 0;
% DT = 1;
% for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         % Claim the S variable in the objective function
%         fprintf(fid,' %16s',strcat('H',cell2mat(arcstmp{k1,1})));
%         len_arc=length(cell2mat(arcstmp{k1,1})); 
%         arc_tmp=cell2mat(arcstmp{k1,1});
%         fprintf(fid,'              obj %25.4f',0);
%         obj_count=obj_count+1;
%         obj_variable{obj_count,1}=strcat('H',cell2mat(arcstmp{k1,1}));
%         obj_variable{obj_count,2}={num2str(obj_count)};
%         fprintf(fid,'\n');
%         num_H_down =num_H_down +1;
%         if flag_H_down==0 
%             obj_H_down=obj_count;
%             flag_H_down=1;
%         end
%         
%         % H in HG - Down time
%         fprintf(fid,' %16s',strcat('H',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'HG'));
%         fprintf(fid,'%25.4f',1);
%         fprintf(fid,'\n');
% 
%         % H in HL - sigma (H)- down time
%         for c_dt = 1:DT%(num_H_down,1)
%             arc_tmpi = cell2mat(arcstmp{k1+c_dt-1,1});
%             len_arci = length(arc_tmpi);
%             fprintf(fid,' %16s',strcat('H',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(arc_tmpi,'HL'));
%             fprintf(fid,'%25.4f',1);
%             fprintf(fid,'\n'); 
%             if strncmp(arc_tmpi(len_arci-5:len_arci),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
%                 break;
%             end
%         end               
%     end
% end
% %%

for k1=1:length(nodestmp(:,1))
    %%%%%%%%%%%%%% Objective function (Loss of load penalty part)
    %*******************************************************************************************
    % important: the way mps file is changed to matrix form is decided by the sequence that
    % these variables appear in the objective function
    %*******************************************************************    *************************
    if (strncmp(cell2mat(nodestmp{k1,1}),'EL',2)==1) % If it is a electric node
        fprintf(fid,' %16s',strcat('SL',cell2mat(nodestmp{k1,1})));
        tmp_name=cell2mat(nodestmp{k1,1});
        fprintf(fid,'              obj %25s\n',num2str(lol));    % discount rate NOT considered
        obj_count=obj_count+1;
        obj_variable{obj_count,1}=strcat('SL',cell2mat(nodestmp{k1,1}));
        obj_variable{obj_count,2}={num2str(obj_count)};
        if flag_obj_SL==0
            obj_SL=obj_count;
            flag_obj_SL=1;
        end
        %********************************************************************************************
        % Only in the electric part
        %********************************************************************************************
        fprintf(fid,' %16s',strcat('SL',cell2mat(nodestmp{k1,1})));
        fprintf(fid,' %16s                         1\n',cell2mat(nodestmp{k1,1}));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of (Loss of load penalty part)
    end
end


% An external function named sortcell is used
Sorted_angle=sortcell(Angle_tmp);
for kk=1:length(Sorted_angle)
    fprintf(fid,' %16s',cell2mat(Sorted_angle(kk,1)));
    fprintf(fid,' %16s %25s',cell2mat(Sorted_angle(kk,2)),cell2mat(Sorted_angle(kk,3)));
    fprintf(fid,'\n');
end


%   Write section RHS
fprintf(fid,'RHS \n');
k2=0;
for k1=1:length(nodestmp(:,1))
    if str2double(cell2mat(nodestmp{k1,4}))~=0 % the value is not 0; If there is a demand associated to this node
        if k2==0
            fprintf(fid,'\n rhs %16s',cell2mat(nodestmp{k1,1}));
            fprintf(fid,' %25s',cell2mat(nodestmp{k1,4}));
            k2=1;
        else
            fprintf(fid,' %16s',cell2mat(nodestmp{k1,1}));
            fprintf(fid,' %25s',cell2mat(nodestmp{k1,4}));
            k2=0;
        end
    end;
    if strncmp(cell2mat(nodestmp{k1,1}),'ST',2)==1  % Set the initial  storage level
        Num_rhs=Num_rhs+1;
        if Num_rhs==1
            fprintf(fid,'\n rhs %16s',cell2mat(nodestmp{k1,1}));
            fprintf(fid,' %25s',cell2mat(nodestmp{k1,4}));
        end
    end
end;

%% Apr 05 2011 - Regulation amt.
% load Reg_req
% load Regh % Apr 10
% NLVhr = NLVhr(1:dpi,1);%% new
NLVhr = Yregh(1+(ii-1)*dpi:dpi*ii,1);%% new

loop=a*b;
loop_c=0;

for k1=1:length(bb(:,1))  % this part add new constraints for reserve
%% Mar 27 Regulation up and down requirements
%% Apr 05 2011 Reg amt.
    loop_c=loop_c+1;
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'1');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',0.5*Max_gen+3*NLVhr(loop_c,1));
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'3');  % lower bound  >=     non-spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',Max_gen+3*NLVhr(loop_c,1));% Oct 04 2011
    
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'RU');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));%30); % 1% of peak load
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'RD');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));
    % Feb 18 2012 new
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'2');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));%30); % 1% of peak load
    fprintf(fid,'\n rhs %15s%s',cell2mat(bb{k1,1}),'4');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
    fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));
    if loop_c==loop
        loop_c=0;
    end
end

loop_c=0;
loop=a*b;
count=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%     tmp_E=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9}));
%     fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
%     fprintf(fid,' %25.4f',tmp_E);
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'M');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));
    
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'N');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,13}));
        
        %% Nov 4 2011 Min down time
        loop_c=loop_c+1;
        count=count+1;
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'MD');
        fprintf(fid,' %25.4f',-min(DT(count,1),loop+1-loop_c)); 
        if loop_c==loop
            loop_c=0;
        end
    end
end

%% Mar 30 2011 - Compressor input + reg + res equations RHS

for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'L');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'S');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,7}));
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'a');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));

        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'c');
        fprintf(fid,' %25f',str2double(cell2mat(arcstmp{k1,8})));
                
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'d');
        fprintf(fid,' %25f',str2double(cell2mat(arcstmp{k1,8})));
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'e');
        fprintf(fid,' %25f',str2double(cell2mat(arcstmp{k1,14}))/6); % Feb 18 2012 new        
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'R'); % Mar 30
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,14}));
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'D');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,15}));
        
     end
end

% Feb 20 2012 new - down reg restrain                     
loop=a*b;
loop_c=0;
%% Mar 31 2011 - Reservoir with reg + res equations - RHS
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'L');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));
        
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'S');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,7}));
        
        % Mar 18 2012 - Res (t-1)-T> min
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'T');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,7}));
        
                % Mar 18 2012 - Res (t-1)+C> max
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'C');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));

        len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});                
        tmp_name=cell2mat(arcstmp{k1,3});
        
        loop_c=loop_c+1;
        % Feb 20 2012 new - down reg restrain             
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,2}),'5');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
        fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));  
        
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
 %       if(strncmp(tmp_name(5:end),'d02h24',6)==1) %% RHS for final
	        fprintf(fid,'\n rhs %15s%s',strcat(tmp_name(1:4),tmp_name,'L'));
            fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));
        
	        fprintf(fid,'\n rhs %15s%s',strcat(tmp_name(1:4),tmp_name,'S'));
            fprintf(fid,' %25s',cell2mat(arcstmp{k1,7}));
                        
            % Mar 18 2012 - Res (t-1)-T> min
	        fprintf(fid,'\n rhs %15s%s',strcat(tmp_name(1:4),tmp_name,'T'));
            fprintf(fid,' %25s',cell2mat(arcstmp{k1,7}));

            % Mar 18 2012 - Res (t-1)+C> max
	        fprintf(fid,'\n rhs %15s%s',strcat(tmp_name(1:4),tmp_name,'C'));
            fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));
            
                    loop_c=loop_c+1;          
            fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,3}),'5');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
            fprintf(fid,' %25.4f',3*NLVhr(loop_c,1));  
        end
        
        if loop_c==loop
            loop_c=0;
        end
    end
end

%% RHS for Ramp rates - Mar 11 2011
for k1=1:length(arcstmp(:,1))  
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
    fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'R');
    %fprintf(fid,' %25s',num2str(30));%cell2mat(arcstmp{k1,14}));
    fprintf(fid,' %25s',cell2mat(arcstmp{k1,14}));
    fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'D');
    fprintf(fid,' %25s',cell2mat(arcstmp{k1,15}));
    % Feb 18 2012 new
    fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'e');
    fprintf(fid,' %25f',str2double(cell2mat(arcstmp{k1,14}))/6);        
    end
end


% %% April 2 2011 - Down time Sigma(H)
% for k1=1:length(arcstmp(:,1))  
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'HL');
%         fprintf(fid,' %25s',num2str(1));
%     end
% end

        
fprintf(fid,'\n');

%   Write section BOUNDS

fprintf(fid,'\nBOUNDS\n');
%% Apr 05 2011 - add new load
%  load loadhourly
% loadhr = loadhr(1:dpi,1);% new
% load Yearly
loadhr = Yloadh(1+(ii-1)*dpi:dpi*ii,1)*3000/(max(Yloadh));% new
loop_c=0;
loop = a*b;
numl = 0; %% load number - 24 loads
for k1=1:length(arcstmp(:,1))
    if strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==0
        %     if str2double(cell2mat(arcstmp{k1,9}))~=0   % for existing lines
        if str2double(cell2mat(arcstmp{k1,7}))~=0   % The default setting is 0 to infinity, so no need to set up.
            fprintf(fid,' LO bnd %16s',cell2mat(arcstmp{k1,1}));
            if strncmp(cell2mat(arcstmp{k1,3}),'DD',2)==0 %% Apr 05 2011 - load
                tmp=str2num(cell2mat(arcstmp{k1,7}))*str2num(cell2mat(arcstmp{k1,9}));
                if  strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=0
                    fprintf(fid,' %25f\n',0);
                else
                    fprintf(fid,' %25s\n',num2str(tmp));
                end
            else
                %% Apr 05 2011 - add new load
                if loop_c==0
                    numl = numl+1; %% load number
                end
                loop_c=loop_c+1;
                fprintf(fid,' %25f\n',0.6*loadhr(loop_c,1)*stressf(numl,1));
                fprintf(fid,' UP bnd %16s',cell2mat(arcstmp{k1,1}));
                fprintf(fid,' %25f\n',0.6*loadhr(loop_c,1)*stressf(numl,1));
                if loop_c==loop
                    loop_c=0;%% hour number
                end
            %% here load
            end
        end;
        if strncmp(cell2mat(arcstmp{k1,8}),'Inf',3)==0  % If upper bound not infinity            
            if strncmp(cell2mat(arcstmp{k1,3}),'DD',2)==0 %% Apr 05 2011
                tmp=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9}));
                fprintf(fid,' UP bnd %16s',cell2mat(arcstmp{k1,1}));
                fprintf(fid,' %25s\n',num2str(tmp));
            end       
        end;

    elseif (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_L0=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9})); % Uppwer bound
        tmp_L1=str2num(cell2mat(arcstmp{k1,19}))*str2num(cell2mat(arcstmp{k1,9})); % Uppwer bound
        tmp_L2=str2num(cell2mat(arcstmp{k1,22}))*str2num(cell2mat(arcstmp{k1,9})); % Uppwer bound
        tmp_L3=str2num(cell2mat(arcstmp{k1,25}))*str2num(cell2mat(arcstmp{k1,9})); % Uppwer bound
        tmp_S1=str2num(cell2mat(arcstmp{k1,18}))*str2num(cell2mat(arcstmp{k1,9})); % Lower bound
        tmp_S2=str2num(cell2mat(arcstmp{k1,21}))*str2num(cell2mat(arcstmp{k1,9})); % Lower bound
        tmp_S3=str2num(cell2mat(arcstmp{k1,24}))*str2num(cell2mat(arcstmp{k1,9})); % Lower bound 
        fprintf(fid,' UP bnd %16s',strcat('1',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25s\n',num2str(tmp_L1));        
        fprintf(fid,' UP bnd %16s',strcat('2',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25s\n',num2str(tmp_L2));  
        fprintf(fid,' UP bnd %16s',strcat('3',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25s\n',num2str(tmp_L3));
        fprintf(fid,' UP bnd %16s',cell2mat(arcstmp{k1,1}));
        fprintf(fid,' %25s\n',num2str(tmp_L0));
        
        fprintf(fid,' UP bnd %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,12}))));
        fprintf(fid,' UP bnd %16s',strcat('RU',cell2mat(arcstmp{k1,1})));%Mar 27 Regu up
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,32})))); %% Mar 17 2012 
        fprintf(fid,' UP bnd %16s',strcat('RD',cell2mat(arcstmp{k1,1})));%Mar 27 Regu down
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,32}))));
%         fprintf(fid,' UP bnd %16s',strcat('B',cell2mat(arcstmp{k1,1})));%Oct 07 2011
%         fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,12}))));
        fprintf(fid,' UP bnd %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));%Oct 07 2011
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,13}))));
        
        fprintf(fid,' BV bnd %16s',strcat('U',cell2mat(arcstmp{k1,1})));
        fprintf(fid,'\n');
%% Oct 03 2011 Non-Spinning
        fprintf(fid,' BV bnd %16s',strcat('N',cell2mat(arcstmp{k1,1})));
        fprintf(fid,'\n');
        fprintf(fid,' BV bnd %16s',strcat('X',cell2mat(arcstmp{k1,1})));
        fprintf(fid,'\n'); 
        fprintf(fid,' BV bnd %16s',strcat('Y',cell2mat(arcstmp{k1,1})));
        fprintf(fid,'\n');  
    else
        if strncmp(cell2mat(arcstmp{k1,8}),'Inf',3)==0  % If upper bound not infinity
            tmp=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9}));
            fprintf(fid,' UP bnd %16s',cell2mat(arcstmp{k1,1}));
            fprintf(fid,' %25s\n',num2str(tmp));
        end;
    end
end;

for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_L0=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9})); % Uppwer bound
%% Apr 10 - Ramp limits on ancillary
        fprintf(fid,' UP bnd %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25f\n',str2double(cell2mat(arcstmp{k1,12}))*str2double(cell2mat(arcstmp{k1,14}))/6);
        fprintf(fid,' UP bnd %16s',strcat('RU',cell2mat(arcstmp{k1,1})));%Mar 27 Regu up
        fprintf(fid,' %25f\n',str2double(cell2mat(arcstmp{k1,32}))*str2double(cell2mat(arcstmp{k1,14}))/12);
        fprintf(fid,' UP bnd %16s',strcat('RD',cell2mat(arcstmp{k1,1})));%Mar 27 Regu down
        fprintf(fid,' %25f\n',str2double(cell2mat(arcstmp{k1,32}))*str2double(cell2mat(arcstmp{k1,15}))/12);
    end
end
        
%% Mar 31 2011 Reservoir final limit
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        tmp_name=cell2mat(arcstmp{k1,3});
        len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
  %      if(strncmp(tmp_name(5:end),'d02h24',6)==1)
              fprintf(fid,' UP bnd %16s',strcat(tmp_name(1:4),tmp_name));
	        fprintf(fid,' %25s\n',cell2mat(arcstmp{k1,8}));
	  end
    end
end
%%
        
for k2=1:length(nodestmp(:,1))
    if strncmp(cell2mat(nodestmp{k2,1}),'EL',2)==1
        fprintf(fid,' LO bnd %16s',strcat('AG',cell2mat(nodestmp{k2,1})));
        fprintf(fid,' %25s\n',num2str(-3.14));
        fprintf(fid,' UP bnd %16s',strcat('AG',cell2mat(nodestmp{k2,1})));
        fprintf(fid,' %25s\n',num2str(3.14));
    end
end
fprintf(fid,'ENDATA\n');
fclose('all');

% cd ('C:\tomlab')
% startup;
[F_UC1,c_UC1,A_UC1,b_L_UC1,b_U_UC1,x_L_UC1,x_U_UC1,IntVars_UC1] = cpx2mat('C:\tomlab\problem_UC.mps',0);



% Remember: X is determined, then everything evolves around X.
% X is comprised of 1. arcs (num_arcs) 2. loss of node penalty part (num_bus) 3. angle (num_bus)
% clear nodestmp;
% clear arcstmp;


% [x_uc,slack_uc,v_uc,nouse_uc,objv_uc,nouse2_uc,nouse3_uc,inform_uc] = cplex(c_UC,A_UC,x_L_UC,x_U_UC,b_L_UC,b_U_UC,[],[],[],[],IntVars_UC);
