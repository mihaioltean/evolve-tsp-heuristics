//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2022.01.27.0

//   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//---------------------------------------------------------------------------
//   Each subpopulation is evolved in a separate thread
//   The subpopulations have a circular structure
//   At the end of each generation we move, from each subpopulation, some individuals in the next one

//   I recommend to check the basic variant first (without subpopulations and threads)

//   Compiled with Microsoft Visual C++ 2019
//   Requires C++11 or newer (for thread support)

//   how to use it: 
//   just create a console project and copy-paste the content this file in the main file of the project

//   More info at:  https://www.mepx.org

//   Please reports any sugestions and/or bugs to:     mihai.oltean@gmail.com

//   Training data file must have the following format (see building1.txt and cancer1.txt):
//   building1 and cancer1 data are taken from PROBEN1

//   x11 x12 ... x1n f1
//   x21 x22 ....x2n f2
//   .............
//   xm1 xm2 ... xmn fm

//   where m is the number of training data
//   and n is the number of variables.

//--------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <thread>
#include <mutex>

#include "mep_chromosome.h"
#include "graphs.h"
#include "mep_rands.h"
#include "tsp_params.h"

// variables
// some indexes in the variables array


//---------------------------------------------------------------------------
void generate_random_chromosome(t_mep_chromosome &a, const t_mep_parameters &params, t_seed& seed) // randomly initializes the individuals
{
	// generate constants first
	for (unsigned int c = 0; c < params.num_constants; c++)
		a.constants[c] = mep_real_rand(seed, params.constants_min, params.constants_max);

	// on the first position we can have only a variable or a constant
	double sum = params.variables_probability + params.constants_probability;
	double p = mep_real_rand(seed, 0, sum);

	if (p <= params.variables_probability)
		a.prg[0].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
	else
		a.prg[0].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1);

	a.prg[0].addr1 = a.prg[0].addr2 = a.prg[0].addr3 = a.prg[0].addr4 = 0;
	// for all other genes we put either an operator, variable or constant
	for (unsigned int i = 1; i < params.code_length; i++) {
		double p = mep_real_rand(seed, 0, 1);

		if (p <= params.operators_probability)
			a.prg[i].op = -1 -mep_unsigned_int_rand(seed, 0, NUM_OPERATORS - 1);       // an operator
		else {
			if (p <= params.operators_probability + params.variables_probability)
				a.prg[i].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);     // a variable
			else
				a.prg[i].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1); // index of a constant
		}
		a.prg[i].addr1 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr2 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr3 = mep_unsigned_int_rand(seed, 0, i - 1);
		a.prg[i].addr4 = mep_unsigned_int_rand(seed, 0, i - 1);
	}
}
//---------------------------------------------------------------------------
void mutation(t_mep_chromosome &a_chromosome, const t_mep_parameters &params, t_seed& seed) // mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = mep_real_rand(seed, 0, 1);
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		double p = mep_real_rand(seed, 0, sum);

		if (p <= params.variables_probability)
			a_chromosome.prg[0].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
		else
			a_chromosome.prg[0].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1);
	}
	// other genes
	for (unsigned int i = 1; i < params.code_length; i++) {
		p = mep_real_rand(seed, 0, 1);      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = mep_real_rand(seed, 0, 1);

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -1 - mep_unsigned_int_rand(seed, 0, NUM_OPERATORS - 1);
			else
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.prg[i].op = mep_unsigned_int_rand(seed, 0, num_variables - 1);
				else
					a_chromosome.prg[i].op = num_variables + mep_unsigned_int_rand(seed, 0, params.num_constants - 1); // index of a constant
		}

		p = mep_real_rand(seed, 0, 1);      // mutate the first address  (addr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr1 = mep_unsigned_int_rand(seed, 0, i - 1);

		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr2 = mep_unsigned_int_rand(seed, 0, i - 1);
		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr3)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr3 = mep_unsigned_int_rand(seed, 0, i - 1);
		p = mep_real_rand(seed, 0, 1);      // mutate the second address   (addr4)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].addr4 = mep_unsigned_int_rand(seed, 0, i - 1);
	}
	// mutate the constants
	for (unsigned int c = 0; c < params.num_constants; c++) {
		p = mep_real_rand(seed, 0, 1);
		if (p < params.mutation_probability)
			a_chromosome.constants[c] = mep_real_rand(seed, params.constants_min, params.constants_max);
	}
}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2, t_seed& seed)
{
	unsigned int cutting_pct = 1 + mep_unsigned_int_rand(seed, 0, params.code_length - 3);
	for (unsigned int i = 0; i < cutting_pct; i++) {
		offspring1.prg[i] = parent1.prg[i];
		offspring2.prg[i] = parent2.prg[i];
	}
	for (unsigned int i = cutting_pct; i < params.code_length; i++) {
		offspring1.prg[i] = parent2.prg[i];
		offspring2.prg[i] = parent1.prg[i];
	}
	// now the constants
	if (params.num_constants) {
		cutting_pct = 1 + mep_unsigned_int_rand(seed, 0, params.num_constants - 3);
		for (unsigned int i = 0; i < cutting_pct; i++) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		for (unsigned int i = cutting_pct; i < params.num_constants; i++) {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
	}
}
//---------------------------------------------------------------------------
void uniform_crossover(const t_mep_chromosome &parent1, const t_mep_chromosome &parent2, 
	const t_mep_parameters &params, 
	t_mep_chromosome &offspring1, t_mep_chromosome &offspring2, t_seed& seed)
{
	for (unsigned int i = 0; i < params.code_length; i++)
		if (rand_int_01(seed)) {
			offspring1.prg[i] = parent1.prg[i];
			offspring2.prg[i] = parent2.prg[i];
		}
		else {
			offspring1.prg[i] = parent2.prg[i];
			offspring2.prg[i] = parent1.prg[i];
		}

	// constants
	for (unsigned int i = 0; i < params.num_constants; i++)
		if (rand_int_01(seed)) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		else {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
}
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{// comparator for quick sort
	if (((t_mep_chromosome *)a)->fitness > ((t_mep_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_mep_chromosome *)a)->fitness < ((t_mep_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
unsigned int tournament_selection(const t_mep_chromosome *a_sub_pop, int sub_pop_size, unsigned int tournament_size, t_seed& seed)     // Size is the size of the tournament
{
	unsigned int p = mep_unsigned_int_rand(seed, 0, sub_pop_size - 1);
	for (unsigned int i = 1; i < tournament_size; i++) {
		int r = mep_unsigned_int_rand(seed, 0, sub_pop_size - 1);
		p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------

void evolve_one_subpopulation(unsigned int *current_subpop_index, std::mutex* mutex,
		t_mep_chromosome ** sub_populations, unsigned int generation_index, 
	const t_mep_parameters &params, const t_graph *training_graphs, int num_training_graphs, 
	double* vars_values, t_seed* seed)
{
	unsigned int pop_index = 0;
	while (*current_subpop_index < params.num_sub_populations) {// still more subpopulations to evolve?

		while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
		mutex->unlock();

		// pop_index is the index of the subpopulation evolved by the current thread

		t_mep_chromosome *a_sub_population = sub_populations[pop_index];

		t_mep_chromosome offspring1, offspring2;
		offspring1.allocate_chromosome(params);
		offspring2.allocate_chromosome(params);

		double *partial_values_array = new double[params.code_length];

		if (generation_index == 0) {
			for (unsigned int i = 0; i < params.sub_population_size; i++) {
				generate_random_chromosome(a_sub_population[i], params, seed[pop_index]);
				
				a_sub_population[i].compute_fitness(params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

			}
			// sort ascendingly by fitness inside this population
			qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
		}
		else // next generations
			for (unsigned int k = 0; k < params.sub_population_size; k += 2) {
				// we increase by 2 because at each step we create 2 offspring

				// choose the parents using binary tournament
				unsigned int r1 = tournament_selection(a_sub_population, params.sub_population_size, 1, seed[pop_index]);
				unsigned int r2 = tournament_selection(a_sub_population, params.sub_population_size, 1, seed[pop_index]);
				// crossover
				double p_0_1 = mep_real_rand(seed[pop_index], 0, 1); // a random number between 0 and 1
				if (p_0_1 < params.crossover_probability)
					one_cut_point_crossover(a_sub_population[r1], a_sub_population[r2], params, offspring1, offspring2, seed[pop_index]);
				else {// no crossover so the offspring are a copy of the parents
					offspring1.copy_individual(a_sub_population[r1], params);
					offspring2.copy_individual(a_sub_population[r2], params);
				}
				// mutate the result and compute fitness
				mutation(offspring1, params, seed[pop_index]);
				offspring1.compute_fitness(params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

				// mutate the other offspring too
				mutation(offspring2, params, seed[pop_index]);
				offspring2.compute_fitness(params.code_length, training_graphs, num_training_graphs, vars_values, partial_values_array);

				// replace the worst in the population
				if (offspring1.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
					a_sub_population[params.sub_population_size - 1].copy_individual(offspring1, params);
					qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
				}
				if (offspring2.fitness < a_sub_population[params.sub_population_size - 1].fitness) {
					a_sub_population[params.sub_population_size - 1].copy_individual(offspring2, params);
					qsort((void *)a_sub_population, params.sub_population_size, sizeof(a_sub_population[0]), sort_function);
				}
			}

		offspring1.delete_chromosome();
		offspring2.delete_chromosome();
	}
}
//---------------------------------------------------------------------------
void start_steady_state(const t_mep_parameters &params, 
	const t_graph *training_graphs, int num_training_graphs,
	FILE *f_out)
{
	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory for all sub populations
	t_mep_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_mep_chromosome*[params.num_sub_populations];
	for (unsigned int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_mep_chromosome[params.sub_population_size];
		for (unsigned int i = 0; i < params.sub_population_size; i++)
			sub_populations[p][i].allocate_chromosome(params); // allocate each individual in the subpopulation 
	}

	// allocate memory
	double** vars_values = new double*[params.num_threads];
	for (unsigned int t = 0; t < params.num_threads; t++)
		vars_values[t] = new double[num_variables];

	t_seed* seeds = new t_seed[params.num_sub_populations];
	for (unsigned int p = 0; p < params.num_sub_populations; p++)
		seeds[p].init(p);

	// an array of threads. Each sub population is evolved by a thread
	std::thread **mep_threads = new std::thread*[params.num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads

	unsigned int best_individual_subpop_index;
	// evolve for a fixed number of generations
	for (unsigned int generation = 0; generation < params.num_generations; generation++) { // for each generation

		//current_subpop_index = 0;
		unsigned int current_subpop_index = 0;
		for (unsigned int t = 0; t < params.num_threads; t++)
			mep_threads[t] = new std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, generation, params, training_graphs, num_training_graphs, vars_values[t], seeds);

		for (unsigned int t = 0; t < params.num_threads; t++) {
			mep_threads[t]->join();
			delete mep_threads[t];
		}

		// find the best individual
		best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
		for (unsigned int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
				best_individual_subpop_index = p;
		// compute average
		double average_fitness = 0;
		for (unsigned int p = 0; p < params.num_sub_populations; p++)
			for (unsigned int i = 0; i < params.sub_population_size; i++)
				average_fitness += sub_populations[p][i].fitness;

		average_fitness /= (params.num_sub_populations * params.sub_population_size);

		printf("generation %d, best fitness = %lf; average fitness = %lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness, average_fitness);
		fprintf(f_out, "%d %lf %lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness, average_fitness);
		sub_populations[best_individual_subpop_index][0].print_to_file(f_out, params);

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (unsigned int p = 0; p < params.num_sub_populations; p++) {
			unsigned int  k = mep_unsigned_int_rand(seeds[0], 0, params.sub_population_size - 1);// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			unsigned int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				sub_populations[index_next_pop][params.sub_population_size - 1].copy_individual(sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}
	}
	delete[] seeds;
	delete[] mep_threads;

	// print best t_mep_chromosome

	sub_populations[best_individual_subpop_index][0].print_to_screen(params);
	// free memory

	for (unsigned int p = 0; p < params.num_sub_populations; p++) {
		for (unsigned int i = 0; i < params.sub_population_size; i++)
			sub_populations[p][i].delete_chromosome();
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

	for (unsigned int t = 0; t < params.num_threads; t++)
		delete[] vars_values[t];
	delete[] vars_values;
}
//--------------------------------------------------------------------

int main(void)
{
	t_mep_parameters params;
	params.num_sub_populations = 30;
	params.sub_population_size = 100;			// the number of individuals in population  (must be an even number!)
	params.code_length = 50;
	params.num_generations = 200;				// the number of generations
	params.mutation_probability = 0.01;          // mutation probability
	params.crossover_probability = 0.9;         // crossover probability

	params.variables_probability = 0.4;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 10; // use 10 constants from -1 ... +1 interval
	params.constants_min = -1;
	params.constants_max = 1;

	params.num_threads = 30;

	t_graph *training_graphs = NULL;
	int num_training_graphs = 0;

	const char* path_str = "c:/Mihai/work/evolve-tsp-heuristics/data/";

	if (!read_training_data(path_str, training_graphs, num_training_graphs)) {
		printf("Cannot find input file(s) in path %s! Please specify the correct path!", path_str);
		getchar();
		return 1;
	}

	compute_global_variables(training_graphs, num_training_graphs);

	srand(0); 

	printf("evolving...\n");
	// store start time
	std::chrono::time_point<std::chrono::system_clock> start;
	start = std::chrono::system_clock::now();

	FILE* f_out = fopen("c:/temp/tsp_output.txt", "w");

	if (!f_out) {
		printf("Cannot create file for output! Please specify a valid path!\n");
		getchar();
		return 1;
	}

	start_steady_state(params, training_graphs, num_training_graphs, f_out);

	fclose(f_out);

	std::chrono::time_point<std::chrono::system_clock> end;
	end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	printf("run time = %lf\n", elapsed_seconds.count());

	delete_training_graphs(training_graphs, num_training_graphs);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------