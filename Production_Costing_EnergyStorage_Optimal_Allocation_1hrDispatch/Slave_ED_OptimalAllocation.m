%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% SCED (uses loadhourly.mat, Reg_req.mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Program for generating MPS format from
%   nodes.txt and arcs.txt
% This is the optimized version and use the 11 properties arcs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [num_arcs,num_nodes,num_br,num_bus,num_new]= creatmps(N_lines)
% clear

%**************************************************
% Differences between problem_ED and problem_UC
% 1. ED don't have UC variables
% 2. In ED, commented out the Angle_tmp part 
%**************************************************************************
%**************************************************************************
% % C{1,1}='Elec';C{1,2}=[1 12  4 3];C{1,3}={'y' 'm' 'w' 'c'};
% a=1;
% b=24;
% C{1,1}='Elec';C{1,2}=[a b];C{1,3}={'d' 'h'};
% C{2,1}='Nuc';C{2,2}=[a b];C{2,3}={'d' 'h'};
% C{3,1}='Coal';C{3,2}=[a b];C{3,3}={'d' 'h'};
% C{4,1}='Oil';C{4,2}=[a b];C{4,3}={'d' 'h'};
% C{5,1}='Wind';C{5,2}=[a b];C{5,3}={'d' 'h'};
% C{6,1}='Stor';C{6,2}=[a b];C{6,3}={'d' 'h'};
% C{6,1}='NG';C{6,2}=[a b];C{6,3}={'d' 'h'};
% OR1=350;
% 
% T=C{1,2};
% W=C{1,3};
% Num_st=0;
% Num_rhs=0;
% Num_ang=0;
% t_num=1;
% dr=0;
% lol=500;
% num_x=0;  % dummy fuel resource nodes
% num_gen=0;  % num. of generators
% num_fuel=0; % num. of fuel transportation lines (excuding dummy lines)
% num_dm=0;   % num. of demand nodes
% num_nodes=0; % num. of all nodes except dummy nodes
% num_arcs=0;% num. of all arcs
% num_br=0; % number of existing electric branches
% % num_fconst=0; % number of flow constraints
% num_bus=0; % num. of all angles/buses
% num_new=0; % num. of all potential lines
% num_new_DC=0;
% num_wind=0; % num. wind generators
% num_UC=0;
% 
% st_arcs=1;
% st_x=0;          %num_arcs - num_dm - num_br - num_gen - num_x +1;
% st_gen=0;        %st_x + num_x;
% st_br=0;         %st_gen + num_gen;
% st_dm=0;         %st_br + num_br;
% st_new=0;
% st_new_DC=0;
% st_wind=0;
% num_xcew=0;
% st_xcew=0;
% 
% Pctg_nowind=0.07; % percentage of non-wind generation that should be secured by reserve
% Pctg_wind=0.10;   % percentage of wind generation that should be secured by reserve
% Max_gen=350;      % The capacity of the largest gen. unit 
% obj_UC=0;         % The beginning of the UC variable
% obj_UC_num=0;
% obj_pointer=0;
% obj_enebid=0;
% obj_reserve1=0;   % The beginning of the first reserve1(spin) variable 
% obj_reserve2=0;   % The beginning of the first reserve2(non-spin) variable 
% 
% 
% node_wind=0;
% node_prod=0;
% node_stor=0;
% node_tran=0;
% node_gen=0;
% nst_prod=0;
% nst_stor=0;
% nst_tran=0;
% nst_gen=0;
% nst_wind=0;
% expand_num=C{1,2}(1)*C{1,2}(2);
% 
% num_itg=0;      % number of integer variable
obj_row_ED={};
obj_variable_ED={};
obj_count_ED=0; 
row_count_ED=0;

num_br=0;
Num_rhs=0;
%**************************************************************************
% New  variables are anfle= AG+nodename(buses), New branches for DC OPF = BR+branchname(transmission)
%**************************************************************************
% 
% %   Read nodes.txt
% fid=fopen('nodes.txt');
% nodes=textscan(fid,'%s %s %s %s',-1);
% fclose('all');
% for k1=1:length(nodes)
%     for k2=1:length(nodes{1})
%         nodestmp{k2,k1}=nodes{k1}(k2);
%     end;
% end;
% clear nodes;
% %   Read arcs.txt
% fid=fopen('arcs.txt');
% arcs=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,-1);  % 29 colmns
% fclose('all');
% for k1=1:length(arcs)
%     for k2=1:length(arcs{1})
%         arcstmp{k2,k1}=arcs{k1}(k2);
%     end;
% end;
%  clear arcs;
% %   Read arcsinitial.txt
% fid=fopen('arcsinitial.txt');
% arcsinitial=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s' ,-1);  % 29 colmns
% fclose('all');
% for k1=1:length(arcsinitial)
%     for k2=1:length(arcsinitial{1})
%         arcsinitialtmp{k2,k1}=arcsinitial{k1}(k2);
%     end;
% end;
% clear arcsinitial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Writes problem to a MPS file
fid=fopen('C:\tomlab\problem_ED_OptimalAllocation.mps', 'w');
% fid=fopen('C:\Users\vkrish\Documents\MATLAB\tomlab\problem_ED.mps', 'w');

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
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={cell2mat(nodestmp{k1,1})};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};   
    end;
%     % if (strncmp(cell2mat(nodestmp{k1,3}),'Gen',3)==1)&&(strncmp(cell2mat(nodestmp{k1,2}),'Wind',4)~=1)
%     % if it is the first generator arc but not a wind arc
%     if strncmp(cell2mat(nodestmp{k1,1}),'EL',2)==1 % If it is elec node
%         num_bus=num_bus+1;
%         node_tran=node_tran+1;
%         if fg01==0
%             nst_tran=k1;
%             fg01=1;
%         end
%     end;
%     if strncmp(cell2mat(nodestmp{k1,3}),'Prod',4)==1 % If it is prod node
%         node_prod=node_prod+1;
%         if fg02==0
%             nst_prod=k1;
%             fg02=1;
%         end
%     end;
%     if strncmp(cell2mat(nodestmp{k1,3}),'Stor',4)==1 % If it is stor node
%         node_stor=node_stor+1;
%         if fg03==0
%             nst_stor=k1;
%             fg03=1;
%         end
%     end;
%     if strncmp(cell2mat(nodestmp{k1,3}),'Gen',3)==1 % If it is gen node
%         node_gen=node_gen+1;
%         if fg04==0
%             nst_gen=k1;
%             fg04=1;
%         end
%     end;
%     if strncmp(cell2mat(nodestmp{k1,2}),'Wind',4)==1 % If it is wind node
%         node_wind=node_wind+1;
%         if fg06==0
%             nst_wind=k1;
%             fg06=1;
%         end
%     end;
    
end
flag01=0;
flag02=0;
flag03=0;
flag04=0;
flag05=0;
flag06=0;
flag07=0;
flag08=0;

% NEW code, add N equations
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
            &&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==0)  % If it's a transmission line, add N equations   OK
        fprintf(fid,' E %s%s\n','BR',cell2mat(arcstmp{k1,1}));
        num_br=num_br+1;
        Name{num_br,1}=strcat('BR',cell2mat(arcstmp{k1,1}));
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat('BR',cell2mat(arcstmp{k1,1}))};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
    end;
    
%     if (strncmp(cell2mat(arcstmp{k1,1}),'X',1)==1)
%         num_x=num_x+1;
%         if (flag01==0)
%             st_x=k1;
%             flag01=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))
%         num_gen=num_gen+1;
%         if (flag02==0)
%             st_gen=k1;
%             flag02=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)
%         
%         if (flag03==0)
%             st_br=k1;
%             flag03=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'DD',2)==1)
%         num_dm=num_dm+1;
%         if (flag04==0)
%             st_dm=k1;
%             flag04=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
%             &&(str2num(cell2mat(arcstmp{k1,9}))==0)&&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==0)
%         num_new=num_new+1;
%         if (flag05==0)
%             st_new=k1;
%             flag05=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,2}),'EL',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EL',2)==1)...
%             &&(str2num(cell2mat(arcstmp{k1,9}))==0)&&(strncmp(cell2mat(arcstmp{k1,11}),'DC',2)==1)
%         num_new_DC=num_new_DC+1;
%         if (flag08==0)
%             st_new_DC=k1;
%             flag08=1;
%         end
%     end
%     if (strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1)&&(flag06==0)
%         st_wind=k1;
%         flag06=1;
%     end
%     if strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1
%         num_wind=num_wind+1;
%     end
%     if (strncmp(cell2mat(arcstmp{k1,2}),'XC',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'EW',2)==1)
%         num_xcew=num_xcew+1;
%         if (flag07==0)
%             st_xcew=k1;
%             flag07=1;
%         end
%     end
%     
end
% num_br=num_br-num_new;  % In this way, num_br is the num. of all existing lines
% invcost=zeros(num_new,1); % the investment cost vector
% line_num=zeros(num_new,1); % number of lines for each kind of potential lines
% line_sus=zeros(num_new,1); % susceptance of each potential line
% bound=zeros(num_new,1);   % bound of each potential line
% for i_tmp=st_new:st_new+num_new-1
%     j_tmp=i_tmp-st_new+1;
%     line_sus(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,11}));
%     line_num(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,9}));
%     invcost(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,10}));
%     bound(j_tmp,1)=str2double(cell2mat(arcstmp{i_tmp,8}));
% end
obj_row_StorTur=0;

for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc segment
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))
%         fprintf(fid,' E %s%s\n',cell2mat(arcstmp{k1,1}),'E');  % this row is defined before. Arcname=Eng.bid 1 + Eng. bid 2 + Eng. bid 3
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for reserve 1
%         fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Max. cap for reserve 2
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'a');  % Max. cap for reserve 1
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'b');  % Max. cap for reserve 2
%         fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'M');  %

        % Feb 23 2012 Shadow P2 -- 10 constraints
          if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
              obj_row_StorTur=row_count_ED+1;
          end
          
          fprintf(fid,' E %s%s\n',cell2mat(arcstmp{k1,1}),'E');  % this row is defined before. Arcname=Eng.bid 1 + Eng. bid 2 + Eng. bid 3
                row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'E')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for reserve 1
                row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'M'); 
                                row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'M')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};         

        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'a');  % Max. cap for reserve 1
                row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'a')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};       
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'b');  % Max. cap for reserve 2
                        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'b')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 

        %   the above line is differnet from
        %   fprintf(fid,'E%s%s\n','E',cell2mat(nodestmp{k1,1}));,which means the KCLfor each node
        
        %% Generator Ramping - Mar 11 2011 - Ramp up & Ramp down
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'R');  %Ramp up
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'R')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'D');  %Ramp down
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'D')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
        
        
        %% Generator Regulation - Mar 27 2011 - Reg up & down
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'c');  %Reg up
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'c')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'d');  %Reg down
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'d')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
        
        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'e');  %SR+RU+NSR<10*rr
        row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'e')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; % Feb 18 2012 new        
    end
end

obj_row_S_pter_ED=row_count_ED+1; % Shadow P3 - Turbine min. cons.
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc segment
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&((strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1))
        fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Max. cap for reserve 2
                row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
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

obj_row_1=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
%     fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'1');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
%     fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'2');  % lower bound  >=     spin reserve >= 0.5*percentage of gen.output
%     fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'3');  % lower bound  >=     non-spin reserve >= 0.5*largest gen. unit
%     fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'4');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
     fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'1');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'1')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end

obj_row_2=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'2');  % lower bound  >=     spin reserve >= 0.5*percentage of gen.output
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'2')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end

obj_row_3=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'3');  % lower bound  >=     non-spin reserve >= 0.5*largest gen. unit
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'3')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end

obj_row_4=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
    fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'4');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'4')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end
        
        %% Mar 27 2011 - Regulation up and down
        obj_row_ru=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve
        fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'RU');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'RU')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end

obj_row_rd=row_count_ED+1;
for k1=1:length(bb(:,1))  % this part add new constraints for reserve        
        fprintf(fid,' G %s%s\n',cell2mat(bb{k1,1}),'RD');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
                    row_count_ED=row_count_ED+1;
        obj_row_ED{row_count_ED,1}={strcat(cell2mat(bb{k1,1}),'RD')};
        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
end

obj_row_DRlim=row_count_ED+1; % Shadow P4
%% Feb 20 2012 - Restrain down regulation from storage (else extra accumulation)
for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,2}),'5');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,2}),'5')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
        
              %% DR comp + UR tur <= rating Oct 10 2012
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,2}),'6');  % lower bound  >=     non-spin reserve >= 0.5*percentage of gen.output
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,2}),'6')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
       
        tmp_name=cell2mat(arcstmp{k1,3});
       len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1                    
%% Feb 20 2012 - down reg restrain 
           fprintf(fid,' L %s\n',strcat(tmp_name,'5'));  
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name,'5')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
           
           fprintf(fid,' L %s\n',strcat(tmp_name,'6'));  
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name,'6')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};
       end
   end
end

%% Mar 30 2011 - Add compressor reserve and regulation capability
obj_row_StorCom=row_count_ED+1; % Shadow P5-- 8
for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
%        fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Compressor with regu/res
%        row_count_ED=row_count_ED+1;
%        obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
%        obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Min. cap for Compressor with regu/res
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};        
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'a');  % Max. cap for reserve 1
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'a')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};       
      
        
        %% Compressor Regulation - Mar 30 2011 - Reg up & down
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'c');  %Reg up
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'c')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
        
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'d');  %Reg down
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'d')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'e');  %RU+SR<10*rr
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'e')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};  % Feb 18 2012 new       
       
       %% CAES Compressor Ramping - Wed, Mar 30 2011 - Ramp up & Ramp down
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'R');  %Ramp up
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'R')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
        
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'D');  %Ramp down
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'D')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
   end 
end

obj_row_StorComLDJ=row_count_ED+1; % Oct 12 2012 - COmp limit - disj
for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)~=1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Compressor with regu/res
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
   end
end

%%
%% Mar 31 2011 - Add Storage reservoir limits with Creg- and Treg+, Tres+ - ROW
obj_row_StorRes=row_count_ED+1; % Shadow P6 - 2
for k1=1:length(arcstmp(:,1))  
   if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'L');  % Max. cap for Reservoir
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'L')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'S');  % Min. cap for Reservoir
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'S')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};     
       
       fprintf(fid,' G %s%s\n',cell2mat(arcstmp{k1,1}),'T');  % Mar 18 2012 - Storage Res - min for Turb.DR
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'T')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};     
       
       fprintf(fid,' L %s%s\n',cell2mat(arcstmp{k1,1}),'C');  % Mar 18 2012 - Storage Res - max for Comp RU, SR
       row_count_ED=row_count_ED+1;
       obj_row_ED{row_count_ED,1}={strcat(cell2mat(arcstmp{k1,1}),'C')};
       obj_row_ED{row_count_ED,2}={num2str(row_count_ED)}; 
       
        len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        tmp_name=cell2mat(arcstmp{k1,3});
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1;
%       if(strncmp(tmp_name(5:end),'d02h24',6)==1) %% ROW for final
           fprintf(fid,' L %s\n',strcat(tmp_name(1:4),tmp_name,'L'));  % Max. cap for Reservoir
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name(1:4),tmp_name,'L')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
           
           fprintf(fid,' G %s\n',strcat(tmp_name(1:4),tmp_name,'S'));  % Min. cap for Reservoir
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name(1:4),tmp_name,'S')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};     
           
           fprintf(fid,' G %s\n',strcat(tmp_name(1:4),tmp_name,'T'));  % Mar 18 2012 - Storage Res - min for Turb.DR
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name(1:4),tmp_name,'T')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};    
           
           fprintf(fid,' L %s\n',strcat(tmp_name(1:4),tmp_name,'C'));  % Mar 18 2012 - Storage Res - max for Comp RU, SR
           row_count_ED=row_count_ED+1;
           obj_row_ED{row_count_ED,1}={strcat(tmp_name(1:4),tmp_name,'C')};
           obj_row_ED{row_count_ED,2}={num2str(row_count_ED)};                       
       end
       
   end 
end
%%

%   Write section COLUMNS
fprintf(fid,'COLUMNS\n');

loop = a*b; %mar 11
loop_c=0;

for k1=1:length(arcstmp(:,1))
    %%%%%%%%%%%%%% Objective function
    %********************************************************************************************
    % important: When mps file is changed to matrix form, the sequence of variable vector x is decided by the sequence that
    % these variables appear in the objective function.
    % In the objective function part, all transmission lines, existing and petential, are included.
    % Another important rule: The sequence of the equations is decided by
%     % the sequence of row name appearing in "section ROWS" of mps file.
    %********************************************************************************************
    %     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)~=1)||(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)==1)
    fprintf(fid,' %16s',cell2mat(arcstmp{k1,1}));
    fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,5}))+ str2double(cell2mat(arcstmp{k1,30}))*carbon_tax);
    fprintf(fid,'\n');
        obj_count_ED=obj_count_ED+1;
    obj_variable_ED{obj_count_ED,1}=cell2mat(arcstmp{k1,1});
    obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
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
%         bij=str2num(cell2mat(arcstmp{k1,9}))*str2num(cell2mat(arcstmp{k1,11}));
%         Num_ang=Num_ang+1;
%         Angle_tmp(Num_ang,1)={strcat('AG',cell2mat(arcstmp{k1,2}))};
%         Angle_tmp(Num_ang,2)={strcat('BR',cell2mat(arcstmp{k1,1}))};
%         Angle_tmp(Num_ang,3)={num2str(-bij)};
%         Num_ang=Num_ang+1;
%         Angle_tmp(Num_ang,1)={strcat('AG',cell2mat(arcstmp{k1,3}))};
%         Angle_tmp(Num_ang,2)={strcat('BR',cell2mat(arcstmp{k1,1}))};
%         Angle_tmp(Num_ang,3)={num2str(bij)};
        
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
                fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
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
        fprintf(fid,' %25.4f',-Pctg_wind);% Feb 18 2012 new
        fprintf(fid,'\n');
    end
end

% for k1=1:length(arcstmp(:,1))   % This part define the startup and shutdown cost
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         % First define the variable in the objective function 
%         fprintf(fid,' %16s',strcat('X',cell2mat(arcstmp{k1,1}))); % First energy bidding
%         fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,16})));
%         fprintf(fid,'\n');     
% %         fprintf(fid,' %16s',strcat('X',cell2mat(arcstmp{k1,1}))); % First energy bidding
% %         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));
% %         fprintf(fid,' %25.4f',-1);
% %         fprintf(fid,'\n');      
%     end
% end
% for k1=1:length(arcstmp(:,1))  % This part define the startup and shutdown cost
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         fprintf(fid,' %16s',strcat('Y',cell2mat(arcstmp{k1,1}))); % First energy bidding
%         fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,17})));
%         fprintf(fid,'\n'); 
%         fprintf(fid,' %16s',strcat('Y',cell2mat(arcstmp{k1,1}))); % First energy bidding
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));
%         fprintf(fid,' %25.4f',1);
%         fprintf(fid,'\n'); 
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following part define the energy bidding variables (3 for each gen.)
flag_enebid=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1}))); % First energy bidding
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,20})));
        fprintf(fid,'\n');
                obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('1',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
%         % Energy bid + reserve 1 + reserve 2 -max.capacity*U <0
%         fprintf(fid,' %16s',strcat('1',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,' %25f',1);
%         fprintf(fid,'\n');
%         % Energy bid + reserve 1 + reserve 2 -min.capacity*U >0
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
                        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('2',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
%         Energy bid + reserve 1 + reserve 2 -max.capacity*U <0
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
                      obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('3',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
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
num_res1=0;
st_spin =0;

st_caes_spin =0;
flag_reserve1_Caes=0;
num_res1_caes=0;
st_caes_nonspin =0;
flag_reserve2_Caes=0;
num_res2_caes=0;

loop = a*b; %mar 11
loop_c=0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
        tmp_name=cell2mat(arcstmp{k1,2});
        % Claim this variable in the objective function
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        %% Apr 05 2011 - Spin bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
        %%
        fprintf(fid,'\n');
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('A',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        
        num_res1=num_res1+1;
        if flag_reserve1==0 
            st_spin=obj_count_ED;
            flag_reserve1=1;            
        end
        
        if strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1
            num_res1_caes=num_res1_caes+1;
            if flag_reserve1_Caes==0 
                st_caes_spin=obj_count_ED;
                flag_reserve1_Caes=1;
            end
        end
        
%         % Energy bid + reserve 1  -max.capacity*U <0
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
        
        %% Ramp up - Mar 11 2011 in Spining reserve
        % Ramp up --> gen(t) - gen(t-1) + Res1(t) + Res2(t) <= rup;

%         fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); % Apr10
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
                fprintf(fid,' %25f',-1);
                fprintf(fid,'\n');
%           end 
        end
        
        if(loop_c == loop)
            loop_c=0;
        end
        
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

% loop = a*b;
% loop_c=0;
% for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         tmp_name=cell2mat(arcstmp{k1,2});
%         % Claim this variable in the objective function
%         fprintf(fid,' %16s',strcat('B',cell2mat(arcstmp{k1,1})));
%          %% Apr 05 2011 - Non-Spin bid hourly
%         loop_c=loop_c+1;
%         fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
%         %%
%  
%         fprintf(fid,'\n');
%                         obj_count_ED=obj_count_ED+1;
%         obj_variable_ED{obj_count_ED,1}=strcat('B',cell2mat(arcstmp{k1,1}));
%         obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
%         
%         if strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1
%             num_res2_caes=num_res2_caes+1;
%             if flag_reserve2_Caes==0 
%                 st_caes_nonspin=obj_count_ED;
%                 flag_reserve2_Caes=1;
%             end
%         end
% 
% %         % Energy bid -min.capacity*U >0
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
%         %% Ramp up - Mar 11 2011 in NSR
%         % Ramp up --> gen(t) - gen(t-1) + Res1(t) + Res2(t) <= rup;
% %        % loop_c=loop_c+1;
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
% %% Feb 15 2012
% %                 %% Oct 03 2011 Turbine Non Sp.res inside reservoir
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
% % %             else % comment
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
num_reg = 0;
flag_reg_Caes=0;
st_caes_reg=0;
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
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('RU',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        if flag_obj_RU ==0
            obj_RU =obj_count_ED;
            flag_obj_RU =1;
        end
        num_reg = num_reg+1;
        
        if strncmp(cell2mat(arcstmp{k1,4}),'GenS',4)==1
            if flag_reg_Caes==0 
                st_caes_reg=obj_count_ED;
                flag_reg_Caes=1;
            end
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

       
        %% Ramp up - Mar 11 2011 in Spin reserve
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
                fprintf(fid,' %25f',-1);
                fprintf(fid,'\n');
%            end
        end
        
        if(loop_c == loop)
            loop_c=0;
        end
        
         %% Oct 10 2012- DR comp + UR tur<rating
        if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)~=1)
                 fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,2}),'6'));
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');      
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
%         %%
        
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
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('RD',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        if flag_obj_RD==0
            obj_RD=obj_count_ED;
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
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');
  %          end
                  %% Feb 20 2012- restrain down reg in stor
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,2}),'5'));
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');      
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
                %% Apr 05 2011 - Spin bid hourly
        loop_c=loop_c+1;
        fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,27})));%spin_bid(loop_c));
        %%
 
        fprintf(fid,'\n');
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('A',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        num_comA = num_comA+1;
        if flag_obj_comA==0
            obj_comA=obj_count_ED;
            flag_obj_comA=1;
        end
        
        % Compr - up Reg - Res > min.
        fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1})));        % print the reserve 1 name
        fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
        fprintf(fid,' %25f',-1);
        fprintf(fid,'\n');
     
        % Sum of Reserve biddings > OR1, reserve 1 bidding from compr.
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
                fprintf(fid,' %25s',strcat('-',cell2mat(arcstmp{k1,6})));
                fprintf(fid,'\n');
%           end 
%        end
        
                %% Mar 18 2012 Turbine Down regulation inside reservoir (t-1)       
             % Reservoir(t-1) + Com Reg up(t)+ com SR (t) <max
            %if(loop_c>1)
                resr_name = cell2mat(arcstmp{k1,3});
                fprintf(fid,' %16s',strcat('A',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'C'));
                fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
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
loop_c = 0;
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
%        fprintf(fid,'              obj %25.4f',30);
        fprintf(fid,'\n');
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('RU',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        num_comRU = num_comRU+1;
        if flag_obj_comRU==0
            obj_comRU=obj_count_ED;
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
                fprintf(fid,' %25s',strcat('-',cell2mat(arcstmp{k1,6})));
                fprintf(fid,'\n');
%            end
       % end
       
                       %% Mar 18 2012 Turbine Down regulation inside reservoir (t-1)       
             % Reservoir(t-1) + Com Reg up(t)+ com SR (t) <max
            %if(loop_c>1)
                resr_name = cell2mat(arcstmp{k1,3});
                fprintf(fid,' %16s',strcat('RU',cell2mat(arcstmp{k1,1}))); 
                fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'C'));
                fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
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
loop_c = 0;
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
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('RD',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        num_comRD = num_comRD+1;
        if flag_obj_comRD==0
            obj_comRD=obj_count_ED;
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
                fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
                fprintf(fid,'\n');
                
               %%Feb 20 2012 Restrain down reg in stor
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,3}),'5'));
                fprintf(fid,' %25f',1);
                fprintf(fid,'\n');           
                
                %%Oct 10 2012 RD comp + RU tur<rating
                fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 2's bidding limit
                fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,3}),'6'));
                fprintf(fid,' %25s',cell2mat(arcstmp{k1,6}));
                fprintf(fid,'\n'); 
%            end
       % end
       
        if(loop_c == loop)
            loop_c=0;
        end
        
    %% Feb 15 2012            
        %% Mar 31 2011 Compressor Down regulation in reservoir
%         resr_name=cell2mat(arcstmp{k1,3});
%         fprintf(fid,' %16s',strcat('RD',cell2mat(arcstmp{k1,1})));        % Reserve 1's bidding limit
%         fprintf(fid,' %16s',strcat(resr_name(1:4),resr_name,'L'));
%         fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,6})));
%         fprintf(fid,'\n');
%         
    end
end
%%             

%% Mar 31 2011 Storage reservoir final - STxx STxx d<final> h<final> 
obj_final_st = 0;
flag_final_st = 0;
final_c = 0;
for k1=1:length(arcstmp(:,1))
    if (strncmp(cell2mat(arcstmp{k1,2}),'ST',2)==1)&&(strncmp(cell2mat(arcstmp{k1,3}),'ST',2)==1)
        len_arc=length(cell2mat(arcstmp{k1,3}));
        arc_tmp=cell2mat(arcstmp{k1,3});
                
        tmp_name=cell2mat(arcstmp{k1,3});
        if strncmp(arc_tmp(len_arc-5:len_arc),strcat('d',sprintf('%02.0f',d),'h',sprintf('%02.0f',b)),6)==1
%        if(strncmp(tmp_name(5:end),'d02h24',6)==1)
            final_c = final_c+1;
            % Claim this variable in the objective function
            fprintf(fid,' %16s',strcat(tmp_name(1:4),tmp_name));
            fprintf(fid,'              obj %25.4f',str2double(cell2mat(arcstmp{k1,5})));
            fprintf(fid,'\n');
            obj_count_ED=obj_count_ED+1;
            obj_variable_ED{obj_count_ED,1}=strcat(tmp_name(1:4),tmp_name);
            obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
            if flag_final_st==0
                obj_final_st=obj_count_ED;
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
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('Bn',cell2mat(arcstmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        if flag_obj_Bn==0
            obj_Bn=obj_count_ED;
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
                fprintf(fid,' %25f',-1);
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
% flag_UC=0;
% for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
%     if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%         % Claim the UC variable in the objective function
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         len_arc=length(cell2mat(arcstmp{k1,1})); 
%         arc_tmp=cell2mat(arcstmp{k1,1});
%         fprintf(fid,'              obj %25.4f',0);
%         fprintf(fid,'\n');
%         obj_pointer=obj_pointer+1;
%         if flag_UC==0 
%             obj_UC=obj_pointer;
%         end
%         flag_UC=1;
%         
%         % This part deal with the equation U(t)-U(t+1)=X-Y startup and shutdown cost
%         if strncmp(arc_tmp(len_arc-2:len_arc),strcat('h',sprintf('%02.0f',b)),3)~=1
%               fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%               fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'U'));   
%               fprintf(fid,'%25.4f',1);
%               fprintf(fid,'\n');
%             if strncmp(arc_tmp(len_arc-2:len_arc),'h01',3)~=1
%               fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%               fprintf(fid,' %16s',strcat(arc_tmp(1:len_arc-2),sprintf('%02.0f',str2num(arc_tmp(len_arc-1:len_arc))-1),'U'));
%               fprintf(fid,'%25.4f',-1);
%               fprintf(fid,'\n');
%               
%             end
%         else
%               fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%               fprintf(fid,' %16s',strcat(arc_tmp(1:len_arc-2),sprintf('%02.0f',str2num(arc_tmp(len_arc-1:len_arc))-1),'U')); 
%               fprintf(fid,'%25.4f',-1);
%               fprintf(fid,'\n');
%         end
%                 
%             
%         % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if reserve.indicator1 ==1)
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'L'));
%         fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,8})));
%         fprintf(fid,'\n');
%         
%         % Energy bid  - min.capacity >0
%         fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'S'));
%         fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,7})));
%         fprintf(fid,'\n');
%         
%         if str2num(cell2mat(arcstmp{k1,12}))==1
%             fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
%             fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,8})));
%             fprintf(fid,'\n');
%         end
%         if str2num(cell2mat(arcstmp{k1,13}))==1
%             fprintf(fid,' %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%             fprintf(fid,' %16s',strcat(cell2mat(arcstmp{k1,1}),'b'));
%             fprintf(fid,'%25.4f',-str2double(cell2mat(arcstmp{k1,8})));
%             fprintf(fid,'\n');
%         end
%         
%     end
% end

flag_loadunmet = 0;
obj_loadunmet = 0;
num_loadunmet = 0;
for k1=1:length(nodestmp(:,1))
    %%%%%%%%%%%%%% Objective function (Loss of load penalty part)
    %********************************************************************************************
    % important: the way mps file is changed to matrix form is decided by the sequence that
    % these variables appear in the objective function
    %********************************************************************************************
    if (strncmp(cell2mat(nodestmp{k1,1}),'EL',2)==1) % If it is a electric node
        fprintf(fid,' %16s',strcat('SL',cell2mat(nodestmp{k1,1})));
        tmp_name=cell2mat(nodestmp{k1,1});
        fprintf(fid,'              obj %25s\n',num2str(lol));    % discount rate NOT considered
        %% Apr 05 2011 - See the load shedding variable
        num_loadunmet = num_loadunmet +1;
        obj_count_ED=obj_count_ED+1;
        obj_variable_ED{obj_count_ED,1}=strcat('SL',cell2mat(nodestmp{k1,1}));
        obj_variable_ED{obj_count_ED,2}={num2str(obj_count_ED)};
        if flag_loadunmet==0
            obj_loadunmet=obj_count_ED;
            flag_loadunmet=1;
        end
        
        %%
            
        %********************************************************************************************
        % Only in the electric part
        %********************************************************************************************
        fprintf(fid,' %16s',strcat('SL',cell2mat(nodestmp{k1,1})));
        fprintf(fid,' %16s                         1\n',cell2mat(nodestmp{k1,1}));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of (Loss of load penalty part)
    end
end


% An external function named sortcell is used
% Sorted_angle=sortcell(Angle_tmp);
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
load Reg_req
% load Regh % Apr 10
% NLVhr = NLVhr(1:dpi,1);%% new
% NLVhr = NLVhravg(1+(ii-1)*dpi:dpi*ii,1);%% new
NLVhr = NLVhravg(1:dpi,1);%% new
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
uc_s = 0;
for k1=1:length(arcstmp(:,1))  % This part adds two rows for each generator arc
    if (strncmp(cell2mat(arcstmp{k1,4}),'Gen',3)==1)&&(strncmp(cell2mat(arcstmp{k1,4}),'GenW',4)~=1)
%     tmp_E=str2num(cell2mat(arcstmp{k1,8}))*str2num(cell2mat(arcstmp{k1,9}));
%     fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'E'));
%     fprintf(fid,' %25.4f',tmp_E);

% Mar 28 2011 - Update UC into ED
            uc_s = uc_s+1;



    fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'M');
    fprintf(fid,' %25s',cell2mat(arcstmp{k1,8}));

 
    
    % This part define the max.cap*U and min.cap*U in arcname+L and arcname+S
    % max.cap*U and min.cap*U are numbers now, so they are RHS
    % Energy bid + reserve 1 + reserve 2 - max.capacity <0 (only useful if reserve.indicator1 ==1)
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'L');
%        fprintf(fid,' %25.4f',str2double(cell2mat(arcstmp{k1,9}))*UC_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8}))); 
        fprintf(fid,' %25.4f',str2double(cell2mat(arcstmp{k1,9}))*str2double(cell2mat(arcstmp{k1,8}))); 
         
    % Energy bid  - min.capacity >0
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'S');
%        fprintf(fid,' %25.4f',str2double(cell2mat(arcstmp{k1,9}))*UC_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,7}))); 
        fprintf(fid,' %25.4f',str2double(cell2mat(arcstmp{k1,9}))*str2double(cell2mat(arcstmp{k1,7}))); 
        
        
           fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'a'));
%            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*UC_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8})));           
            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,12}))*str2double(cell2mat(arcstmp{k1,14}))/6);           
            fprintf(fid,'\n');
            
            %% Mar 27 Regulation            
            fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'c')); % up regu
%            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*UC_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8})));
            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*str2double(cell2mat(arcstmp{k1,14}))/12);      
            fprintf(fid,'\n');
             
            fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'d')); % down regu
%            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*UC_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8})));
            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*str2double(cell2mat(arcstmp{k1,15}))/12);            
            fprintf(fid,'\n');
         
             fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'b')); %% Non-spin requires U?
%             fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*Nsp_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8})));
            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,13}))*str2double(cell2mat(arcstmp{k1,14}))/6);            
              fprintf(fid,'\n');
              
    % Feb 18 2012 new
             fprintf(fid,'\n rhs %16s',strcat(cell2mat(arcstmp{k1,1}),'e')); %% Non-spin requires U?
%             fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,9}))*Nsp_status(uc_s,1)*str2double(cell2mat(arcstmp{k1,8})));
            fprintf(fid,'%25.4f',str2double(cell2mat(arcstmp{k1,14}))/6);            
              fprintf(fid,'\n');    

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
        
        % Oct 10 2012 new - RD comp +RU tur< rating
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,2}),'6');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
        fprintf(fid,' %25.4f',20*0.9);
        
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
            
        % Oct 10 2012 new - RD comp +RU tur< rating
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,3}),'6');  % lower bound  >=     spin reserve >= 0.5*largest gen. unit
        fprintf(fid,' %25.4f',20*0.9);
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
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,14}));
        fprintf(fid,'\n rhs %15s%s',cell2mat(arcstmp{k1,1}),'D');
        fprintf(fid,' %25s',cell2mat(arcstmp{k1,15}));
    end
end

fprintf(fid,'\n');

%   Write section BOUNDS
%% Apr 05 2011 - add new load
load loadhourly
loadhr = loadhr(1:dpi,1);% new
% load Yearly
% loadhr = Yloadh(1+(ii-1)*dpi:dpi*ii,1)*3000/(max(Yloadh));% new

loop_c=0;
loop = a*b;
numl = 0; %% load number - 24 loads

fprintf(fid,'\nBOUNDS\n');
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
        
        %% Apr 10 - ramp constraint
        fprintf(fid,' UP bnd %16s',strcat('A',cell2mat(arcstmp{k1,1})));
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,12}))));
        fprintf(fid,' UP bnd %16s',strcat('RU',cell2mat(arcstmp{k1,1})));%Mar 27 Regu up
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,32})))); % Mar 17 2012
        fprintf(fid,' UP bnd %16s',strcat('RD',cell2mat(arcstmp{k1,1})));%Mar 27 Regu down
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,32})))); % Mar 17 2012
%         fprintf(fid,' UP bnd %16s',strcat('B',cell2mat(arcstmp{k1,1})));%Oct 07 2011
%         fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,12}))));
        fprintf(fid,' UP bnd %16s',strcat('Bn',cell2mat(arcstmp{k1,1})));%Oct 07 2011
        fprintf(fid,' %25s\n',num2str(tmp_L0 * str2num(cell2mat(arcstmp{k1,13}))));
%         fprintf(fid,' BV bnd %16s',strcat('U',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,'\n');
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
%% Apr 10 - ramp for AS
%         fprintf(fid,' UP bnd %16s',strcat('A',cell2mat(arcstmp{k1,1})));
%         fprintf(fid,' %25s\n',num2str(tmp_L0));
%         fprintf(fid,' UP bnd %16s',strcat('RU',cell2mat(arcstmp{k1,1})));%Mar 30 Regu up
%         fprintf(fid,' %25s\n',num2str(tmp_L0));
%         fprintf(fid,' UP bnd %16s',strcat('RD',cell2mat(arcstmp{k1,1})));%Mar 30 Regu down
%         fprintf(fid,' %25s\n',num2str(tmp_L0));%% Apr 10 - ramp constraint       
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
%        if(strncmp(tmp_name(5:end),'d02h24',6)==1)
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
[F_ED,c_ED,A_ED,b_L_ED,b_U_ED,x_L_ED,x_U_ED,IntVars_ED] = cpx2mat('C:\tomlab\problem_ED_OptimalAllocation.mps',0);
% A=sparse(A);
F_ED=sparse(F_ED);

% Should consider saving the needed information and erase all of the other variables
% Remember: X is determined, then everything evolves around X.
% X is comprised of 1. arcs (num_arcs) 2. loss of node penalty part
% (num_bus) 3. angle (num_bus)

% [x_uc,slack_uc,v_uc,nouse_uc,objv_uc,nouse2_uc,nouse3_uc,inform_uc] = cplex(c_UC,A_UC,x_L_UC,x_U_UC,b_L_UC,b_U_UC,[],[],[],[],IntVars_UC);
