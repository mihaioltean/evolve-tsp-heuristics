#ifndef mep_chromosome_H
#define mep_chromosome_H

//---------------------------------------------------------------------------
#include <stdio.h>
//---------------------------------------------------------------------------
#include "graphs.h"

//---------------------------------------------------------------------------
#define NUM_OPERATORS 6

#define OP_ADD -1
#define OP_SUB -2
#define OP_MUL -3
#define OP_MIN -4
#define OP_MAX -5
#define OP_IFABCD -6


//---------------------------------------------------------------------------
struct t_mep_parameters {
	unsigned int code_length;             // number of instructions in a t_mep_chromosome
	unsigned int num_generations;
	unsigned int num_sub_populations;       // number of subpopulations
	unsigned int sub_population_size;                // subpopulation size
	double mutation_probability, crossover_probability;
	unsigned int num_constants;
	double constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

	unsigned int num_threads; // num threads. 
	//for best performances the number of subpopulations should be multiple of num_threads.
	// num_thread should no exceed the number of processor cores.
};
//---------------------------------------------------------------------------
struct t_code3 {
	int op;				// either a variable, operator or constant; 
	// variables are indexed from 0: 0,1,2,...; 
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int addr1, addr2, addr3, addr4;    // pointers to arguments
};
//---------------------------------------------------------------------------
class t_mep_chromosome {
private:
	double evaluate(double* vars_values,
		double* partial_values_array, int output_gene_index);
	double compute_error_for_a_gene(int gene_index,
		const t_graph* training_graphs, int num_training_graphs,
		double* vars_values, double* partial_values_array);

public:
	t_code3* prg;        // the program - a string of genes
	double* constants; // an array of constants

	double fitness;        // the fitness (or the error)
	int best_index;
	t_mep_chromosome();
	~t_mep_chromosome();
	void print_to_file(FILE* f, const t_mep_parameters& param);
	void print_to_screen(const t_mep_parameters& param);
	void allocate_chromosome(const t_mep_parameters& params);
	void delete_chromosome(void);
	void copy_individual(const t_mep_chromosome& source, const t_mep_parameters& params);
	void compute_fitness(int code_length,
		const t_graph* training_graphs, int num_training_graphs,
		double* vars_values, double* partial_values_array);
	double compute_path_length_for_a_graph(int gene_index,
		const t_graph& graph, double* vars_values, double* partial_values_array);

	void from_current_file_pos(FILE* f, t_mep_parameters& params);
	bool from_file(const char* file_name, t_mep_parameters& params);
};
//---------------------------------------------------------------------------
#endif