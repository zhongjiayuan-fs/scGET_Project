clc;
clear;
close all;
load HP.mat;
 for i=1:length(cell_gene_name)
    curr_cell=cell_gene_name{i};
    gene_num=length(curr_cell);
    node_HP(i,node_HP(i,:)<0)=0;
    [node_sorted_H,idx]=sort(node_HP(i,1:gene_num));
    node_SH(i)=sum(node_sorted_H(gene_num-300+1:gene_num));
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
















































