//---------------------------------------------------------------------------
//   Evolve TSP heuristics with Genetic Programming - with multiple subpopulations and threads and MPI
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com), Virginia Niculescu (vniculescu@cs.ubbcluj.ro)
//   Version 2016.06.26.0 // year.month.day.build#

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------
//   Each subpopulation is evolved in a separate thread
//   The subpopulations have a circular structure
//   At the end of each generation we move, from each subpopulation, some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations and threads)

//   Compiled with Microsoft Visual C++ 2013
//   Also compiled with XCode 5.
//   Requires C++11 or newer (for thread support)

//   how to use it:
//   just create an empty console project and add this file to the project

//   More info at:  https://github.com/mihaioltean/evolve-tsp-heuristics

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com or to vniculescu@cs.ubbcluj.ro


//--------------------------------------------------------------------


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <time.h>


#include <float.h>

//#include <vld.h> // for detecting memory leaks in VC++



#define max_file_name 50

#define num_operators 7

// +   -1
// -   -2
// *   -3
// min   -4
// max   -5

#define O_ADDITION -1
#define O_SUBTRACTION -2
#define O_MULTIPLICATION -3
#define O_MIN -4
#define O_MAX -5
#define O_SIN -6
#define O_COS -7

char operators_string[num_operators][10] = { "+", "-", "*", "min", "max", "sin", "cos" };

// variables
// some indecses in the variables array

#define min_distance_in_matrix 0     // 0 - min distance in the distance matrix
#define max_distance_in_matrix 1     // 1 - max distance in the distance matrix
#define average_distance_in_matrix 2 // 2 - average distance in the matrix
#define length_so_far 3              // 3 - length so far
#define distance_to_next_node 4      // 4 - distance to the next node
#define distance_to_nearest_node 5   // 5 - distance to the nearest node
#define distance_to_fartest_node 6        // 6 - distance to the fartest node
#define average_distance_to_unvisited 7
#define num_visited_nodes 8          // number of visisted nodes so far
#define num_total_nodes 9            // total number of nodes of that graph


//---------------------------------------------------------------------------
// problem parameters
//---------------------------------------------------------------------------
struct t_parameters{
	int code_length;             // number of instructions in a t_chromosome
	int num_generations;
	int num_sub_populations;       // number of subpopulations
	int sub_population_size;                // subpopulation size
	double mutation_probability, crossover_probability;
	int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;
	int num_migrations;
	int num_training_graphs;
};
//---------------------------------------------------------------------------
// related to chromosomes
//---------------------------------------------------------------------------
struct t_code3{
	int op;				// either a variable, operator or constant;
	// variables are indexed from 0: 0,1,2,...;
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int adr1, adr2;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_chromosome{
	t_code3 *prg;        // the program - a string of genes
	double *constants;   // an array of constants
	int code_length;

	double fitness;      // the fitness (or the error)
	int num_constants;

	//------------------------------------------------------------------------------
	t_chromosome(void)
	{
		prg = NULL;
		constants = NULL;
		num_constants = 0;
		code_length = 0;
	}
	//------------------------------------------------------------------------------
	void from_string_simplified(char* s_source)
	{
		int num_consumed = 0;

		// read chromosome
		sscanf(s_source, "%d%n", &code_length, &num_consumed);
		s_source += num_consumed;
		if (prg)
			delete[] prg;
		prg = new t_code3[code_length];
		for (int i = 0; i < code_length; i++) {
			sscanf(s_source, "%d%d%d%n", &prg[i].op, &prg[i].adr1, &prg[i].adr2, &num_consumed);
			s_source += num_consumed;
		}
		// read constants
		sscanf(s_source, "%d%n", &num_constants, &num_consumed);
		s_source += num_consumed;

		if (constants)
			delete[] constants;
		if (num_constants)
			constants = new double[num_constants];
		else
			constants = NULL;

		for (int i = 0; i < num_constants; i++) {
			sscanf(s_source, "%lf%n", &constants[i], &num_consumed);
			s_source += num_consumed;
		}
		sscanf(s_source, "%lf", &fitness);
	}
	//------------------------------------------------------------------------------
	bool from_file_simplified(char* file_name)
	{
		FILE *f = fopen(file_name, "r");
		if (!f)
			return false;
		char s[1000];
		fgets(s, 1000, f);

		from_string_simplified(s);

		fclose(f);
		return true;
	}
	//------------------------------------------------------------------------------
};
//---------------------------------------------------------------------------
void delete_chromosome(t_chromosome &c)
{
	if (c.prg) {
		delete[] c.prg;
		c.prg = NULL;
	}

	if (c.constants) {
		delete[] c.constants;
		c.constants = NULL;
	}
}
//---------------------------------------------------------------------------
double evaluate(t_code3 *prg, int head_index, int num_variables, double *vars_values, double *partial_values_array, double *constants)
{
	for (int i = 0; i <= head_index; i++)
		switch (prg[i].op) {
		case O_ADDITION:// +
			partial_values_array[i] = partial_values_array[prg[i].adr1] + partial_values_array[prg[i].adr2];
			break;
		case O_SUBTRACTION:// -
			partial_values_array[i] = partial_values_array[prg[i].adr1] - partial_values_array[prg[i].adr2];
			break;
		case O_MULTIPLICATION:// *
			partial_values_array[i] = partial_values_array[prg[i].adr1] * partial_values_array[prg[i].adr2];
			break;
		case O_MAX:// max
			partial_values_array[i] = partial_values_array[prg[i].adr1] > partial_values_array[prg[i].adr2] ? partial_values_array[prg[i].adr1] : partial_values_array[prg[i].adr2];
			break;
		case O_MIN:// min
			partial_values_array[i] = partial_values_array[prg[i].adr1] < partial_values_array[prg[i].adr2] ? partial_values_array[prg[i].adr1] : partial_values_array[prg[i].adr2];
			break;
		case O_SIN:// sin
			partial_values_array[i] = sin(partial_values_array[prg[i].adr1]);
			break;
		case O_COS:// cos
			partial_values_array[i] = cos(partial_values_array[prg[i].adr1]);
			break;
		default:
			if (prg[i].op < num_variables)
				partial_values_array[i] = vars_values[prg[i].op];
			else
				partial_values_array[i] = constants[prg[i].op - num_variables];
	}

	return partial_values_array[head_index]; // last gene is the one providing the output
}
//--------------------------------------------------------------------
//related to test graphs
//--------------------------------------------------------------------

struct t_graph{
	int num_nodes;			// number of nodes
	double **distance;      // distance matrix
	double min_distance, max_distance, average_distance;
	double optimal_length;  // known optimal solution
};
//---------------------------------------------------------------------------
bool allocate_training_graphs(t_graph *&training_graphs, int num_training_graphs)
{
	training_graphs = new t_graph[num_training_graphs];

	return true;
}
//---------------------------------------------------------------------------

int read_graphs(t_graph *training_graphs)
{
	// training is done on 4 graphs
	//4 graphs are read
	int k = 0; // count the graphs
	{
		FILE* f = fopen("data//bayg29.tsp", "r");
		if (!f)
			return 0;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];
		// now read the data
		for (int i = 0; i < training_graphs[k].num_nodes - 1; i++)
			for (int j = i; j < training_graphs[k].num_nodes; j++)
				if (i != j) {
					fscanf(f, "%lf", &training_graphs[k].distance[i][j]);
					training_graphs[k].distance[j][i] = training_graphs[k].distance[i][j];
				}
				else
					training_graphs[k].distance[i][i] = 0;

		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);

		fclose(f);
	}
	{
		k++;
		FILE* f = fopen("data//a280.tsp", "r");
		if (!f)
			return 0;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}
	{
		k++;
		FILE* f = fopen("data//berlin52.tsp", "r");
		if (!f)
			return 0;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}
	{
		k++;
		FILE* f = fopen("data//bier127.tsp", "r");
		if (!f)
			return 0;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double*[training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double *x = new double[training_graphs[k].num_nodes];
		double *y = new double[training_graphs[k].num_nodes];
		int index;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			fscanf(f, "%d%lf%lf", &index, &x[i], &y[i]);
		// now read the length of the shortest path
		fscanf(f, "%lf", &training_graphs[k].optimal_length);
		fclose(f);

		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				training_graphs[k].distance[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));

		delete[] x;
		delete[] y;
	}	return 1;
}
//---------------------------------------------------------------------------
void delete_training_graphs(t_graph *&training_graphs, int num_training_graphs)
{
	if (training_graphs) {
		for (int i = 0; i < num_training_graphs; i++) {
			for (int j = 0; j < training_graphs[i].num_nodes; j++)
				delete[] training_graphs[i].distance[j];
			delete[] training_graphs[i].distance;
		}
		delete[] training_graphs;
	}
}

//--------------------------------------------------------------------
// here we compute the min, max and average distance in a given graph
//--------------------------------------------------------------------
void compute_global_variables(t_graph *training_graphs, int num_training_graphs)
{

	for (int k = 0; k < num_training_graphs; k++) {
		training_graphs[k].max_distance = training_graphs[k].distance[0][1];
		training_graphs[k].min_distance = training_graphs[k].distance[0][1];
		training_graphs[k].average_distance = 0;
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			for (int j = 0; j < training_graphs[k].num_nodes; j++)
				if (i != j) {
					training_graphs[k].average_distance += training_graphs[k].distance[i][j];
					if (training_graphs[k].max_distance < training_graphs[k].distance[i][j])
						training_graphs[k].max_distance = training_graphs[k].distance[i][j];
					if (training_graphs[k].min_distance > training_graphs[k].distance[i][j])
						training_graphs[k].min_distance = training_graphs[k].distance[i][j];
				}
		training_graphs[k].average_distance /=
			(training_graphs[k].num_nodes * training_graphs[k].num_nodes) - training_graphs[k].num_nodes;
	}
}

//---------------------------------------------------------------------------
void calibrate_vars_values(double * vars_values, t_graph *training_graphs, int k)
{
	vars_values[min_distance_in_matrix] /= training_graphs[k].max_distance; //in order to reduce to [0,1] range
	vars_values[max_distance_in_matrix] /= training_graphs[k].max_distance; //in order to reduce to [0,1] range
	vars_values[average_distance_in_matrix] /= training_graphs[k].max_distance;//in order to reduce to [0,1] range

	vars_values[num_total_nodes] /= training_graphs[k].num_nodes; //in order to reduce to [0,1] range

	vars_values[length_so_far] /= training_graphs[k].max_distance;
	vars_values[num_visited_nodes] /= training_graphs[k].num_nodes;
	vars_values[distance_to_next_node] /= training_graphs[k].max_distance;

}
//---------------------------------------------------------------------------
void compute_local_variables(t_graph &graph, int num_visited, int current_node, int *node_visited, double *vars_values)
{
	vars_values[distance_to_nearest_node] = 1E307; // min distance to unvisited
	vars_values[distance_to_fartest_node] = 0;     // max_distance to unvisited
	vars_values[average_distance_to_unvisited] = 0;  // average distance from current to all unvisited

	for (int i = 0; i < graph.num_nodes; i++)
		if (!node_visited[i]) {
			if (vars_values[distance_to_nearest_node] > graph.distance[current_node][i])
				vars_values[distance_to_nearest_node] = graph.distance[current_node][i];
			if (vars_values[distance_to_fartest_node] < graph.distance[current_node][i])
				vars_values[distance_to_fartest_node] = graph.distance[current_node][i];

			vars_values[average_distance_to_unvisited] += graph.distance[current_node][i];
		}

	vars_values[average_distance_to_unvisited] /= (double)(graph.num_nodes - num_visited);
}
//---------------------------------------------------------------------------

void run(t_chromosome &individual, t_graph *graphs, int num_graphs,
	int num_variables, double * vars_values, double *partial_values_array)

{
	// fitness is the sum of errors over all training graphs.
	// error is the distance from the optimum

	// we cannot send the entire graph as a parameter to the function
	// so we extract only some summarized information from it
	// we fill the "vars_values" array with this summarized information
	// a function having vars_values as a parameter will be evolved

	individual.fitness = 0;
	printf("graph #% || path length || shortest path length || error percent\n");

	for (int k = 0; k < num_graphs; k++) {

		// set the value of variables
		vars_values[min_distance_in_matrix] = graphs[k].min_distance;
		vars_values[max_distance_in_matrix] = graphs[k].max_distance;
		vars_values[average_distance_in_matrix] = graphs[k].average_distance;
		vars_values[num_total_nodes] = graphs[k].num_nodes;

		// initialize memory for the path

		int *tsp_path = new int[graphs[k].num_nodes];
		int *node_visited = new int[graphs[k].num_nodes];
		for (int i = 0; i < graphs[k].num_nodes; i++)
			node_visited[i] = 0;

		// init path
		int count_nodes = 0;
		tsp_path[0] = 0; // first node in the path is 0;
		node_visited[0] = 1;
		count_nodes++;
		double path_length = 0;

		while (count_nodes < graphs[k].num_nodes) {
			// fill the path node by node
			double min_eval = DBL_MAX;
			int best_node = -1;
			// compute other statistics
			vars_values[length_so_far] = path_length;
			vars_values[num_visited_nodes] = count_nodes;
			compute_local_variables(graphs[k], count_nodes, tsp_path[count_nodes - 1], node_visited, vars_values);

			for (int node = 0; node < graphs[k].num_nodes; node++)
				// consider each unvisited node
				if (!node_visited[node]) {// not visited yet
					vars_values[distance_to_next_node] = graphs[k].distance[tsp_path[count_nodes - 1]][node];

					calibrate_vars_values(vars_values, graphs, k);

					double eval = evaluate(individual.prg, individual.code_length - 1, num_variables, vars_values, partial_values_array, individual.constants);
					if (eval < min_eval) {
						best_node = node; // keep the one with minimal evaluation
						min_eval = eval;
					}
				}
			// add the best next node to the path
			path_length += graphs[k].distance[tsp_path[count_nodes - 1]][best_node];
			node_visited[best_node] = 1;

			tsp_path[count_nodes] = best_node;
			count_nodes++;

		}
		// connect the last with the first in the path
		path_length += graphs[k].distance[tsp_path[count_nodes - 1]][tsp_path[0]];
		individual.fitness += (path_length - graphs[k].optimal_length) / graphs[k].optimal_length * 100;
		// keep it in percent from the optimal solution, otherwise we have scalling problems
		printf("[%d]  %lf  %lf %lf\n", k, path_length, graphs[k].optimal_length, (path_length - graphs[k].optimal_length) / graphs[k].optimal_length * 100);

		delete[] tsp_path;
		delete[] node_visited;
	}
	individual.fitness /= (double)num_graphs; // average over the number of training graphs
}
//--------------------------------------------------------------------
int main(int argc, char* argv[])
{
	t_graph *graphs = NULL;
	int num_graphs = 4;

	allocate_training_graphs(graphs, num_graphs);
	if (!read_graphs(graphs)) {
		printf("Cannot find input graphs! Please specify the full path!\n");
		printf("Press Enter ...");
		getchar();
		return 1;
	}

	compute_global_variables(graphs, num_graphs);

	int num_variables = 10;

	double* vars_values = new double[num_variables];

	t_chromosome a_chromosome;
	if (!a_chromosome.from_file_simplified("best.txt")) {
		printf("Cannot find input solution! Please specify the full path!\n");
		printf("Press Enter ...");
		getchar();
		return 1;
	}
	double *partial_values_array = new double[a_chromosome.code_length];
	run(a_chromosome, graphs, num_graphs, num_variables, vars_values, partial_values_array);
	delete[] partial_values_array;



	delete[] vars_values;
	delete_training_graphs(graphs, num_graphs);

	printf("Press Enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------