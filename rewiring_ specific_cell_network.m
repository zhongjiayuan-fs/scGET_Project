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

mprofile=log(1+mprofile);

time_0h=mprofile(:,1:92);
time_12h=mprofile(:,93:194);
time_24h=mprofile(:,195:260);
time_36h=mprofile(:,261:432);
time_72h=mprofile(:,433:570);
time_96h=mprofile(:,571:758);

cell_network_0h=csnet(time_0h);
cell_network_12h=csnet(time_12h);
cell_network_24h=csnet(time_24h);
cell_network_36h=csnet(time_36h);
cell_network_72h=csnet(time_72h);
cell_network_96h=csnet(time_96h);



cell_network=cell_network_0h;
for i=93:194
    cell_network{1,i}=cell_network_12h{1,i-92};
end
for i=195:260
    cell_network{1,i}=cell_network_24h{1,i-194};
end
for i=261:432
    cell_network{1,i}=cell_network_36h{1,i-260};
end
for i=433:570
    cell_network{1,i}=cell_network_72h{1,i-432};
end
for i=571:758
    cell_network{1,i}=cell_network_96h{1,i-570};
end

pretime=clock;
for c=1:length(cell_network)
    clear cell;
    clear cell_full_matrix;
    cell=cell_network{1,c};
    cell_full_matrix=full(cell);
    fid=fopen(['tmp_network',num2str(c),'_adj_edges_all.txt'],'wt');
    cell_full_size=size(cell_full_matrix);
    for i=1:cell_full_size(1)
        clear nei_node;
        clear node_expression;
        nei_node=find(cell_full_matrix(i,:)==1);
        if length(nei_node)==0 ||mprofile(i,c)==0
            fprintf(fid,'%g\n',i);
            continue;
        end
        
        node_expression=mprofile(nei_node,c);
        node_zero=find(node_expression==0);
        nei_node(node_zero)=[];
        
        if length(nei_node)==0
            fprintf(fid,'%g\n',i);
            continue;
        end
        
        for j=1:length(nei_node)
            if  j==length(nei_node)
                fprintf(fid,'%g\n',nei_node(j));
            else
                fprintf(fid,'%g\t',nei_node(j));
            end
        end
    end
    fclose(fid);
    currtime=clock;
    c,etime(currtime,pretime)
end





 

