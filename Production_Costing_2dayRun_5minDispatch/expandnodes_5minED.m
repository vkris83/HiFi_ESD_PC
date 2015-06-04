%**************************************************************************
%    Production Costing Program - Energy Storage Integration
%    2010-2014 (c) Dr. Trishna Das & Dr. Venkat Krishnan
%    Iowa State University
%**************************************************************************
% Expands the system node data in nodesinitial.txt to multiperiods (default 48
% hours, though it can changed by changing variables a=#days and b=#hours)
% Used for running Slave_ED_5minED

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Set expansion conditions
%   C{:,2} and C{:,3} represents the time decomposition
%   If, for example:
%   C{1,1}='Elec';C{1,2}=[52 7 24];C{1,3}={'w' 'd' 'h'}
%   Then the electric system will be decomposed in 52 weeks, 7 days/week,
%   and 24 hours/day, 

% C{1,1}='Elec';C{1,2}=[10 12];C{1,3}={'y' 'm' };
% C{2,1}='NG';C{2,2}=[10 12];C{2,3}={'y' 'm'};
% C{3,1}='Coal';C{3,2}=[10];C{3,3}={'y' };
% C{4,1}='Oil';C{4,2}=[10 12];C{4,3}={'y' 'm'};
%w=1;
a_5=48; % days
b_5=12; % hours
C{1,1}='Elec';C{1,2}=[a_5 b_5];C{1,3}={'h' 'm'};
C{2,1}='Nuc';C{2,2}=[a_5 b_5];C{2,3}={'h' 'm'};
C{3,1}='Coal';C{3,2}=[a_5 b_5];C{3,3}={'h' 'm'};
C{4,1}='Oil';C{4,2}=[a_5 b_5];C{4,3}={'h' 'm'};
C{5,1}='Wind';C{5,2}=[a_5 b_5];C{5,3}={'h' 'm'};
C{6,1}='Stor';C{6,2}=[a_5 b_5];C{6,3}={'h' 'm'};
C{6,1}='NG';C{6,2}=[a_5 b_5];C{6,3}={'h' 'm'};

% C{1,1}='Elec';C{1,2}=[w a b];C{1,3}={'w' 'd' 'h'};
% C{2,1}='Nuc';C{2,2}=[w a b];C{2,3}={'w' 'd' 'h'};
% C{3,1}='Coal';C{3,2}=[w a b];C{3,3}={'w' 'd' 'h'};
% C{4,1}='Oil';C{4,2}=[w a b];C{4,3}={'w' 'd' 'h'};
% C{5,1}='Wind';C{5,2}=[w a b];C{5,3}={'w' 'd' 'h'};
% C{6,1}='Stor';C{6,2}=[w a b];C{6,3}={'w' 'd' 'h'};
% C{6,1}='NG';C{6,2}=[w a b];C{6,3}={'w' 'd' 'h'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Read nodes.txt
fid=fopen('nodesinitial.txt');
nodes=textscan(fid,'%s %s %s %s',-1);
fclose('all');
for k1=1:length(nodes)
    for k2=1:length(nodes{1})
        nodestmp{k2,k1}=nodes{k1}(k2);
    end;
end;
nodes=nodestmp;clear nodestmp;

%   Nodes expansion
nodestmp={};
for k1=1:length(nodes(:,1))                 %   k1: Node
    
    for k2=1:length(C(:,1))                 %   k2: Type of energy
        if strcmp(nodes{k1,2},C{k2,1})==1
            T=C{k2,2};
            W=C{k2,3};
        end;
    end;

    a={};
    a(1,:)=nodes(k1,:);

    for k3=1:length(T)                      %   k3: Time category
        c={};
        for k5=1:length(a(:,1))             %   k5: Element of a
            b={};
            for k4=1:T(k3)                  %   k4: number of repetition
                b(k4,1:4)=a(k5,1:4);        %   changed from 1:5 to 1:4
                b{k4,1}={[cell2mat(b{k4,1}) cell2mat(W(k3)) num2str(k4,'%02.0f')]};
               % b(k4,6)={[a{k5,6} k4]};
            end;
            c=[c;b];
        end;
        a=c;
    end;
    nodestmp=[nodestmp;a];
end;
clear a b c T W;


%   Writes nodestmp to a file

fid=fopen('nodestmp_5minED.txt', 'w');
for k1=1:length(nodestmp(:,1))
    for k2=1:length(nodestmp(1,:))
        if iscell(nodestmp{k1,k2})==1
            fprintf(fid,'%s ',cell2mat(nodestmp{k1,k2}));
        else
            fprintf(fid,'%g ',nodestmp{k1,k2});
        end;
    end;
    fprintf(fid,'\n');
end;
fclose('all');

fid=fopen('nodes_5minED.txt', 'w');
for k1=1:length(nodestmp(:,1))
    for k2=1:length(nodestmp(1,:))
        if iscell(nodestmp{k1,k2})==1
            fprintf(fid,'%s ',cell2mat(nodestmp{k1,k2}));
        else
            fprintf(fid,'%g ',nodestmp{k1,k2});
        end;
    end;
    fprintf(fid,'\n');
end;
fclose('all');

