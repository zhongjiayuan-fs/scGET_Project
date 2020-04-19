There gives a scRNA-Seq dataset (hESCs-to-DECs data) to implement the algorithm process.(Other four scRNA-Seq datasets are similar)


rewiring the cell-specific networks
input data: hESCs_to_DECs_data.txt
output data: tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
running pipeline: rewiring_ specific_cell_network.m


Calculating the single-cell graph entropy (SGE)
input data: hESCs_to_DECs_data.txt, tmp_network1_adj_edges_all, tmp_network2_adj_edges_all, ..., tmp_network758_adj_edges_all
running pipeline: SGE.m