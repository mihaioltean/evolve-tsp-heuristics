#ifndef graphs_H
#define graphs_H

//---------------------------------------------------------------------------
class t_graph {
public:
	int num_nodes;			// number of nodes
	double** distance;      // distance matrix
	double min_distance, max_distance, average_distance;
	double optimal_length;  // known optimal solution

	t_graph();
	~t_graph();

	void compute_local_variables(int num_visited, int current_node, const int* node_visited, double* vars_values) const;
};
//---------------------------------------------------------------------------
bool read_training_data(const char* path, t_graph*& training_graphs, int& num_training_graphs);
void delete_training_graphs(t_graph*& training_graphs, int num_training_graphs);
void compute_global_variables(t_graph* training_graphs, int num_training_graphs);
//---------------------------------------------------------------------------

#endif
