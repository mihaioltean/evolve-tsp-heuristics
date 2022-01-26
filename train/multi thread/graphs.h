#ifndef graphs_H
#define graphs_H

//---------------------------------------------------------------------------
struct t_graph {
	int num_nodes;			// number of nodes
	double** distance;      // distance matrix
	double min_distance, max_distance, average_distance;
	double optimal_length;  // known optimal solution
};
//---------------------------------------------------------------------------
bool read_training_data(const char* path, t_graph*& training_graphs, int& num_training_graphs);
void delete_training_graphs(t_graph*& training_graphs, int num_training_graphs);
void compute_global_variables(t_graph* training_graphs, int num_training_graphs);
//---------------------------------------------------------------------------

#endif
