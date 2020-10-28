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

node_HP=zeros(mprofile_size(2),mprofile_size(1));


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
            node_HP(c,na)=0;
        else
            for n=1:length(adjacent_network{na})
                nei=adjacent_network{na}{n};
                e=e+1;
                gy=mprofile(str2num(nei),:);
                curr_weight= abs(edge_weight(gx,gy));
                node_p_rho(n)=abs(curr_weight(c)*mprofile(str2num(nei),c));
            end
            
            for num=1:length(node_p_rho)
                node_p_rho(num)=node_p_rho(num)/(sum(node_p_rho)+es);
            end
            
            for num=1:length(node_p_rho)
                node_HP(c,na)=node_HP(c,na)-(1/length(adjacent_network{na}))*(node_p_rho(num).*log(node_p_rho(num)+es));
            end
            
        end
    end
    currtime=clock;
    c,etime(currtime,pretime)
end


save HP.mat;
















































