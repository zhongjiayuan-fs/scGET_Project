clc;
clear;
close all;


%%

fpi=fopen('hESCs_to_DECs_data.txt');
hline = textscan(fpi, '%s', 1, 'delimiter', '\n');
field=textscan(hline{1}{1},'%s');
clear format;
format='%s';
% format=[format,' %s'];
for i=2:759
    format=[format,' %f'];
end
lines =textscan(fpi, format,1000000,'delimiter', '\t');
gene_names=lines{1};
mprofile = [];
for i = 2 :759
    mprofile = [mprofile, lines{i}];
end
fclose(fpi);

pretime=clock;
mprofile_size=size(mprofile);
cell_num=mprofile_size(2);

pretime=clock;
es=0.00000001;

for c=1:cell_num
    clear adjacent_network;
    fid=fopen(['tmp_network',num2str(c),'_adj_edges_all.txt']);
    adjacent_network={};
    j=0;
    while ~feof(fid)
        tline=fgetl(fid);
        j=j+1;
        adjacent_network{j}=regexp(tline, '\t', 'split');
    end
    fclose(fid);
    total_node_num=j;
    %%
    for na=1:total_node_num
        center=na;
        cell_gene_name{c}{na}=gene_names{center};
        gx=mprofile(center,:);
        e=0;
        rcc=0;
        clear p_rho;
        if (length(adjacent_network{na})==1)&&(str2num(adjacent_network{na}{1})==na)
            HP(c,na)=0;
            node_HP(c,na)=0;
        else
            for n=1:length(adjacent_network{na})
                nei=adjacent_network{na}{n};
                e=e+1;
                gy=mprofile(str2num(nei),:);
                curr_weight= abs(csnedge(gx,gy));
                node_expression=mprofile(str2num(nei),c);
                node_p_rho(n)=abs(curr_weight(c)*node_expression);
            end
            
            node_p_rho=node_p_rho/(sum(node_p_rho)+es);
            node_HP(c,na)=-(1/length(adjacent_network{na}))*sum(node_p_rho.*log(node_p_rho+es));
            
            
        end
    end
    currtime=clock;
    c,etime(currtime,pretime)
end
 
 
 
save HP.mat;
 for i=1:length(cell_gene_name)
    curr_cell=cell_gene_name{i};
    gene_num=length(curr_cell);
    node_HP(i,node_HP(i,:)<0)=0;
    [node_sorted_H,idx]=sort(node_HP(i,1:gene_num));
    node_SH(i)=sum(node_sorted_H(gene_num-150+1:gene_num));
 end


node_result(1)=mean(node_SH(1:92));
node_result(2)=mean(node_SH(93:194));
node_result(3)=mean(node_SH(195:260));
node_result(4)=mean(node_SH(261:432));
node_result(5)=mean(node_SH(433:570));
node_result(6)=mean(node_SH(571:758));


combineData = [node_SH(1:92),node_SH(93:194),node_SH(195:260),node_SH(261:432),node_SH(433:570),node_SH(571:758)];  
group = [2*ones(1,92),3*ones(1,102),4*ones(1,66),5*ones(1,172),6*ones(1,138),7*ones(1,188)]; 
boxplot(combineData,group);
hold on;
t=1:6;
plot(t,node_result,'r','LineWidth',3);
set(gca,'XTick',1:6);
B={'0h' '12h'  '24h' '36h' '72h' '96h'};
set(gca,'XTickLabel',B);
%set(gca,'yticklabel',[0:10]);
xlabel('Stages');
ylabel('SGE');
%plot(t,aver_comidx,'r','LineWidth',3);
title('Average SGE for  hESCs-to-DECs data');
















































