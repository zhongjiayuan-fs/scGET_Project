The implement of SGE on the scRNA-Seq dataset is shown as an example in this project. 
The input data can be changed with the other four datasets if necessary.


Rewiring the cell-specific networks
input data: hESCs_to_DECs_data.txt
output data: tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
running pipeline: Rewiring_ specific_cell_network.m


Calculating the single-cell graph entropy (SGE)
input data: hESCs_to_DECs_data.txt, tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
output data: HP.mat
running pipeline: SGE.m

Displaying result
input data: HP.mat
running pipeline: SGE_main.m
