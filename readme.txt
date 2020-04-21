The implement of SGE on the scRNA-Seq dataset is shown as an example in this project. 
The input data can be changed with the other four datasets if necessary.


rewiring the cell-specific networks
input data: hESCs_to_DECs_data.txt
output data: tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
running pipeline: Rewiring_ specific_cell_network.m


Calculating the single-cell graph entropy (SGE)
input data: hESCs_to_DECs_data.txt, tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
running pipeline: SGE.m
