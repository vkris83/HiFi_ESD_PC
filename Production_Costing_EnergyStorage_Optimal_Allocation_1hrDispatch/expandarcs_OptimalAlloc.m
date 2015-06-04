%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
%Expands the system network flow topology data in arcsinitial.txt to multiperiods (default 48
%hours, though it can changed by changing variables a=#days and b=#hours) 
%%Re-run everytime after changing arcsinitial.txt!

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Set expansion conditions
%   C{:,2} and C{:,3} represents the time decomposition
%   If, for example:
%   C{1,1}='Elec';C{1,2}=[52 7 24];C{1,3}={'w' 'd' 'h'}
%   Then the electric system will be decomposed in 52 weeks, 7 days/week,
%   and 24 hours/day
% w=1;
a=2; % days
b=24; % hours
C{1,1}='Elec';C{1,2}=[a b];C{1,3}={'d' 'h'};
C{2,1}='Nuc';C{2,2}=[a b];C{2,3}={'d' 'h'};
C{3,1}='Coal';C{3,2}=[a b];C{3,3}={'d' 'h'};
C{4,1}='Oil';C{4,2}=[a b];C{4,3}={'d' 'h'};
C{5,1}='Wind';C{5,2}=[a b];C{5,3}={'d' 'h'};
C{6,1}='Stor';C{6,2}=[a b];C{6,3}={'d' 'h'};
C{6,1}='NG';C{6,2}=[a b];C{6,3}={'d' 'h'};
%% CT????? - new

% C{1,1}='Elec';C{1,2}=[w a b];C{1,3}={'w' 'd' 'h'};
% C{2,1}='Nuc';C{2,2}=[w a b];C{2,3}={'w' 'd' 'h'};
% C{3,1}='Coal';C{3,2}=[w a b];C{3,3}={'w' 'd' 'h'};
% C{4,1}='Oil';C{4,2}=[w a b];C{4,3}={'w' 'd' 'h'};
% C{5,1}='Wind';C{5,2}=[w a b];C{5,3}={'w' 'd' 'h'};
% C{6,1}='Stor';C{6,2}=[w a b];C{6,3}={'w' 'd' 'h'};
% C{6,1}='NG';C{6,2}=[w a b];C{6,3}={'w' 'd' 'h'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Read arcs.txt
fid=fopen('arcsinitial_OptimalAlloc.txt');
arcs=textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',-1);  % 34 colmns

fclose('all');
for k1=1:length(arcs)
    for k2=1:length(arcs{1})
        arcstmp{k2,k1}=arcs{k1}(k2);
    end;
end;
arcs=arcstmp;
clear arcstmp;
tn=length(arcs(1,:)); % total number of columns 
%   Arcs expansion
arcstmp={};

for k1=1:length(arcs(:,1))                  %   k1: Arc

    for k2=1:length(C(:,1))                 %   k2: Type of energy
        if strcmp(arcs{k1,4},C{k2,1})==1
            T=C{k2,2};
            W=C{k2,3};
        end;
    end;

    if strncmp(arcs{k1,4},'Gen',3)==1       % If it's a generator arc, find the fuel type
        T=C{1,2};
        W=C{1,3};
        for k2=1:length(C(:,1))             %   k2: Type of energy
            x=cell2mat(arcs{k1,4});
            if strcmp(x(4:5),C{k2,1}(1:2))==1
                T2=C{k2,2};                 %  Find T for the type of energy generator uses
            end;
        end;
    end;

    a={};
    a(1,:)=arcs(k1,:);
    Storage_t={};
    a1={};
    a1(1,:)=arcs(k1,:);

    if strncmp(arcs{k1,2},arcs{k1,3},4)==0      % if it's not a storage line
        for k3=1:length(T)                      %   k3: Time category
            c={};
            for k5=1:length(a(:,1))             %   k5: Element of a  % Number of arcs
                b={};
                if (strncmp(a{k5,4},'Gen',3)==0) % Not a generator arc
                    for k4=1:T(k3)              %   k4: number of repetition
                        b(k4,1:tn)=a(k5,1:tn);
                        b{k4,1}={[cell2mat(b{k4,1}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                        b{k4,2}={[cell2mat(b{k4,2}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                        b{k4,3}={[cell2mat(b{k4,3}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                    end; 
                else                            %  Generator arc
                    for k4=1:T(k3)              %   k4: number of repetition
                        b(k4,1:tn)=a(k5,1:tn);
                        b{k4,1}={[cell2mat(b{k4,1}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                        if k3<=length(T2)
                            % T2 is time range of fuel type, T is the time range of the arc, When K3>Length(T2),
                            % one fuel arc attached to multiple generator arcs
                            b{k4,2}={[cell2mat(b{k4,2}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                        end;
                        b{k4,3}={[cell2mat(b{k4,3}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                    end;
                end;
                c=[c;b];
            end;
            a=c;
        end;
        arcstmp=[arcstmp;a];  % original arcstmp={}
    else   % If it is a storage arc
        for k3=1:length(T)                      %   k3: Time category
            c1={};
            for k5=1:length(a1(:,1))             %   k5: Element of a  % Number of arcs
                b1={};
                for k4=1:T(k3)              %   k4: number of repetition
                    b1(k4,1:tn)=a1(k5,1:tn);
                    b1{k4,1}={[cell2mat(b1{k4,1}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                    b1{k4,2}={[cell2mat(b1{k4,2}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                    b1{k4,3}={[cell2mat(b1{k4,3}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
                end;
                c1=[c1;b1];
            end;
            a1=c1;
        end
        Storage_t=[Storage_t;a1];
        for count=1:(length(Storage_t(:,1))-1)
            Storage_t{count,3}=Storage_t{count+1,3};
        end
        arcstmp=[arcstmp;Storage_t(1:(length(Storage_t(:,1))-1),:)];
    end
end;
clear a b c T W T2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Writes arcstmp to a file

fid=fopen('arcstmp_OptimalAlloc.txt', 'w');
for k1=1:length(arcstmp(:,1))
    for k2=1:length(arcstmp(1,:))
        if iscell(arcstmp{k1,k2})==1
            fprintf(fid,'%s ',cell2mat(arcstmp{k1,k2}));
        else
            fprintf(fid,'%g ',arcstmp{k1,k2});
        end;
    end;
    fprintf(fid,'\n');
end;
fclose('all');

fid=fopen('arcs_OptimalAlloc.txt', 'w');
for k1=1:length(arcstmp(:,1))
    for k2=1:length(arcstmp(1,:))
        if iscell(arcstmp{k1,k2})==1
            fprintf(fid,'%s ',cell2mat(arcstmp{k1,k2}));
        else
            fprintf(fid,'%g ',arcstmp{k1,k2});
        end;
    end;
    fprintf(fid,'\n');
end;
fclose('all');