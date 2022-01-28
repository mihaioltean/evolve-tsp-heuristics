#include <math.h>
#include <float.h>
//---------------------------------------------------------------------------
#include "mep_chromosome.h"
#include "tsp_params.h"
//---------------------------------------------------------------------------
t_mep_chromosome::t_mep_chromosome()
{
	prg = NULL;
	constants = NULL;
	best_index = -1;
	fitness = 0;
}
//---------------------------------------------------------------------------
t_mep_chromosome::~t_mep_chromosome()
{
	delete_chromosome();
}
//---------------------------------------------------------------------------
void t_mep_chromosome::print_to_screen(const t_mep_parameters& params)
{
	const char operators_string[NUM_OPERATORS][10] = { "+", "-", "*", "min", "max", "ifabcd" };

	printf("The t_mep_chromosome is:\n");

	for (unsigned int i = 0; i < params.num_constants; i++)
		printf("constants[%d] = %lf\n", i, constants[i]);

	for (unsigned int i = 0; i < params.code_length; i++)
		if (prg[i].op < 0) {
			if (prg[i].op == OP_IFABCD)
				printf("%d: %s %d %d %d %d\n", i, operators_string[abs(prg[i].op) - 1], prg[i].addr1, prg[i].addr2, prg[i].addr3, prg[i].addr4);
			else
				printf("%d: %s %d %d\n", i, operators_string[abs(prg[i].op) - 1], prg[i].addr1, prg[i].addr2);
		}
		else
			if (prg[i].op < num_variables)
				printf("%d: inputs[%d]\n", i, prg[i].op);
			else
				printf("%d: constants[%d]\n", i, prg[i].op - num_variables);

	printf("Fitness = %lf\n", fitness);
	printf("Best gene = %d\n", best_index);
}
//---------------------------------------------------------------------------
void t_mep_chromosome::print_to_file(FILE* f, const t_mep_parameters& params)
{
	fprintf(f, "%d %d\n", params.code_length, params.num_constants);
	for (unsigned int i = 0; i < params.code_length; i++)
		fprintf(f, "%d %d %d %d %d\n", prg[i].op, prg[i].addr1, prg[i].addr2, prg[i].addr3, prg[i].addr4);

	for (unsigned int i = 0; i < params.num_constants; i++)
		fprintf(f, "%lf ", constants[i]);
	fprintf(f, "\n%d\n", best_index);
}
//---------------------------------------------------------------------------
void t_mep_chromosome::allocate_chromosome(const t_mep_parameters& params)
{
	prg = new t_code3[params.code_length];
	if (params.num_constants)
		constants = new double[params.num_constants];
	else
		constants = NULL;
}
//---------------------------------------------------------------------------
void t_mep_chromosome::delete_chromosome(void)
{
	if (prg) {
		delete[] prg;
		prg = NULL;
	}
	if (constants) {
		delete[] constants;
		constants = NULL;
	}
}
//---------------------------------------------------------------------------
void t_mep_chromosome::copy_individual(const t_mep_chromosome& source, const t_mep_parameters& params)
{
	for (unsigned int i = 0; i < params.code_length; i++)
		prg[i] = source.prg[i];
	for (unsigned int i = 0; i < params.num_constants; i++)
		constants[i] = source.constants[i];
	fitness = source.fitness;
	best_index = source.best_index;
}
//---------------------------------------------------------------------------
double t_mep_chromosome::evaluate(double* vars_values,
	double* partial_values_array, int output_gene_index)
{
	for (int i = 0; i <= output_gene_index; i++)
		switch (prg[i].op) {
		case OP_ADD:// +
			partial_values_array[i] = partial_values_array[prg[i].addr1] + partial_values_array[prg[i].addr2];
			break;
		case OP_SUB:// -
			partial_values_array[i] = partial_values_array[prg[i].addr1] - partial_values_array[prg[i].addr2];
			break;
		case OP_MUL:// *
			partial_values_array[i] = partial_values_array[prg[i].addr1] * partial_values_array[prg[i].addr2];
			break;
		case OP_MAX:// max
			partial_values_array[i] = partial_values_array[prg[i].addr1] > partial_values_array[prg[i].addr2] ? 
									partial_values_array[prg[i].addr1] : partial_values_array[prg[i].addr2];
			break;
		case OP_MIN:// min
			partial_values_array[i] = partial_values_array[prg[i].addr1] < partial_values_array[prg[i].addr2] ? 
									partial_values_array[prg[i].addr1] : partial_values_array[prg[i].addr2];
			break;
		case OP_IFABCD:// a<b?c:d
			partial_values_array[i] = partial_values_array[prg[i].addr1] <
				partial_values_array[prg[i].addr2] ?
				partial_values_array[prg[i].addr3] :
				partial_values_array[prg[i].addr4];
			break;
		default:
			if (prg[i].op < num_variables)
				partial_values_array[i] = vars_values[prg[i].op];
			else
				partial_values_array[i] = constants[prg[i].op - num_variables];
		}

	return partial_values_array[output_gene_index];
}
//---------------------------------------------------------------------------
double t_mep_chromosome::compute_path_length_for_a_graph(int gene_index,
	const t_graph& graph, double* vars_values, double* partial_values_array)
{
	// set the value of variables
	vars_values[min_distance_in_matrix] = graph.min_distance;
	vars_values[max_distance_in_matrix] = graph.max_distance;
	vars_values[average_distance_in_matrix] = graph.average_distance;
	vars_values[num_total_nodes] = graph.num_nodes;

	// initialize memory for the path
	int* tsp_path = new int[graph.num_nodes];
	int* node_visited = new int[graph.num_nodes];
	for (int i = 0; i < graph.num_nodes; i++)
		node_visited[i] = 0;

	// init path
	int count_nodes = 0;
	tsp_path[0] = 0; // first node in the path is 0;
	node_visited[0] = 1;
	count_nodes++;
	double path_length = 0;

	while (count_nodes < graph.num_nodes) {
		// fill the path node by node
		double min_eval = DBL_MAX;
		int best_node = -1;
		// compute other statistics
		vars_values[length_so_far] = path_length;
		vars_values[num_visited_nodes] = count_nodes;
		graph.compute_local_variables(count_nodes, tsp_path[count_nodes - 1], node_visited, vars_values);

		for (int nod = 0; nod < graph.num_nodes; nod++)
			// consider each unvisited node
			if (!node_visited[nod]) {// not visited yet
				vars_values[distance_to_next_node] = graph.distance[tsp_path[count_nodes - 1]][nod];
				double eval = evaluate(vars_values, partial_values_array, gene_index);
				if (eval < min_eval) {
					best_node = nod; // keep the one with minimal evaluation
					min_eval = eval;
				}
			}
		// add the best next node to the path
		path_length += graph.distance[tsp_path[count_nodes - 1]][best_node];
		node_visited[best_node] = 1;

		tsp_path[count_nodes] = best_node;
		count_nodes++;
	}
	// connect the last with the first in the path
	path_length += graph.distance[tsp_path[count_nodes - 1]][tsp_path[0]];

	delete[] tsp_path;
	delete[] node_visited;

	return path_length;
}
//---------------------------------------------------------------------------
double t_mep_chromosome::compute_error_for_a_gene(int gene_index,
	const t_graph* training_graphs, int num_training_graphs,
	double* vars_values, double* partial_values_array)
{
	double error = 0;
	for (int k = 0; k < num_training_graphs; k++) {
		double path_length = compute_path_length_for_a_graph(gene_index,
			training_graphs[k], vars_values, partial_values_array);
		error += (path_length - training_graphs[k].optimal_length) / training_graphs[k].optimal_length * 100;
			
		// keep it in percent from the optimal solution, otherwise we have scalling problems
	}
	error /= (double)num_training_graphs; // average over the number of training graphs
	return error;
}
//---------------------------------------------------------------------------
void t_mep_chromosome::compute_fitness(int code_length,
	const t_graph* training_graphs, int num_training_graphs,
	double* vars_values, double* partial_values_array)
{
	// fitness is the sum of errors over all training graphs.
	// error is the distance from the optimum

	// we cannot send the entire graph as a parameter to the function
	// so we extract only some summarized information from it
	// we fill the "vars_values" array with this summarized information
	// a function having vars_values as a parameter will be evolved

	fitness = DBL_MAX;
	for (int g = 0; g < code_length; g++) {
		double error = compute_error_for_a_gene(g,
			training_graphs, num_training_graphs,
			vars_values, partial_values_array);

		if (fitness > error) {
			fitness = error;
			best_index = g;
		}
	}
}
//-----------------------------------------------------------------
void t_mep_chromosome::from_current_file_pos(FILE* f, t_mep_parameters& params)
{
	fscanf(f, "%d%d", &params.code_length, &params.num_constants);

	if (prg)
		delete[] prg;
	prg = new t_code3[params.code_length];

	for (unsigned int i = 0; i < params.code_length; i++) 
		fscanf(f, "%d%d%d%d%d", &prg[i].op, &prg[i].addr1, &prg[i].addr2, &prg[i].addr3, &prg[i].addr4);

	if (constants)
		delete[] constants;

	if (params.num_constants)
		constants = new double[params.num_constants];
	else
		constants = NULL;

	for (unsigned int i = 0; i < params.num_constants; i++)
		fscanf(f, "%lf", &constants[i]);
	
	fscanf(f, "%d", &best_index);
}
//-----------------------------------------------------------------
bool t_mep_chromosome::from_file(const char* file_name, t_mep_parameters& params)
{
	FILE* f = fopen(file_name, "r");
	if (!f)
		return false;

	from_current_file_pos(f, params);

	fclose(f);
	return true;
}
//------------------------------------------------------------------------------