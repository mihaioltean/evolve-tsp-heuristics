//---------------------------------------------------------------------------
//   Multi Expression Programming Software - with multiple subpopulations and threads
//   Copyright Mihai Oltean  (mihai.oltean@gmail.com)
//   Version 2016.02.03

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
//   just create a console project and copy-paste the content this file in the main file of the project

//   More info at:  http://www.mep.cs.ubbcluj.ro

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


#include "I_operations.hpp"
#include "genetic_mep.hpp"

//#include <vld.h> // for detecting memory leaks in VC++



//---------------------------------------------------------------------------
void print_chromosome(t_chromosome& a, t_parameters &params, int num_variables)
{
    printf("The t_chromosome is:\n");
    
    for (int i = 0; i < params.num_constants; i++)
        printf("constants[%d] = %lf\n", i, a.constants[i]);
    
    for (int i = 0; i < params.code_length; i++)
        if (a.prg[i].op < 0)
            printf("%d: %s %d %d\n", i, operators_string[abs(a.prg[i].op) - 1], a.prg[i].adr1, a.prg[i].adr2);
        else
            if (a.prg[i].op < num_variables)
                printf("%d: inputs[%d]\n", i, a.prg[i].op);
            else
                printf("%d: constants[%d]\n", i, a.prg[i].op - num_variables);
    
    printf("Fitness = %lf\n", a.fitness);
}

//---------------------------------------------------------------------------
void start_steady_state(t_parameters &params, t_graph *training_graphs, int num_training_graphs, int num_variables,  int num_threads)
                        // num threads.
                        //for best performances the number of subpopulations should be multiple of num_threads.
                        // num_thread should no exceed the number of processor cores.)
{
	// a steady state model - 
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

	// allocate memory for all sub populations
	t_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation 
	}

	// allocate memory
	double** vars_values = new double*[num_threads];
	for (int t = 0; t < num_threads; t++)
		vars_values[t] = new double[num_variables];

	// an array of threads. Each sub population is evolved by a thread
	std::thread **mep_threads = new std::thread*[num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads
	
	// initial population (generation 0)
	int current_subpop_index = 0;
	for (int t = 0; t < num_threads; t++)
		mep_threads[t] = new std::thread(evolve_one_subpopulation,
                                         &current_subpop_index,
                                         &mutex, sub_populations,
                                         0, &params,
                                         training_graphs,
                                         num_training_graphs,
                                         num_variables,
                                         vars_values[t]);


	for (int t = 0; t < num_threads; t++) {
		mep_threads[t]->join(); // wait for all threads to execute
		delete mep_threads[t];
	}

	// find the best individual from the entire population
	int best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
	for (int p = 1; p < params.num_sub_populations; p++)
		if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
			best_individual_subpop_index = p;

	printf("generation %d, best fitness = %lf\n", 0, sub_populations[best_individual_subpop_index][0].fitness);
	
	// evolve for a fixed number of generations
	for (int generation = 1; generation < params.num_generations; generation++) { // for each generation

		current_subpop_index = 0;
		for (int t = 0; t < num_threads; t++)
			mep_threads[t] = new std::thread(evolve_one_subpopulation, &current_subpop_index, &mutex, sub_populations, generation, &params, training_graphs, num_training_graphs, num_variables, vars_values[t]);

		for (int t = 0; t < num_threads; t++) {
			mep_threads[t]->join();
			delete mep_threads[t];
		}

		// find the best individual
		best_individual_subpop_index = 0; // the index of the subpopulation containing the best invidual
		for (int p = 1; p < params.num_sub_populations; p++)
			if (sub_populations[p][0].fitness < sub_populations[best_individual_subpop_index][0].fitness)
				best_individual_subpop_index = p;
		printf("generation %d, best fitness = %lf\n", generation, sub_populations[best_individual_subpop_index][0].fitness);

		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (int p = 0; p < params.num_sub_populations; p++) {
			int  k = rand() % params.sub_population_size;// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}
	}

	delete[] mep_threads;

		// print best t_chromosome

	print_chromosome(sub_populations[best_individual_subpop_index][0], params, num_variables);
	// free memory

	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

	for (int t = 0; t < num_threads; t++)
		delete[] vars_values[t];
	delete[] vars_values;
}
//--------------------------------------------------------------------


void init_params(t_parameters& params){

    params.num_sub_populations = 4;
    params.sub_population_size = 50;						    // the number of individuals in population  (must be an even number!)
    params.code_length = 50;
    params.num_generations = 100;					// the number of generations
    params.mutation_probability = 0.1;              // mutation probability
    params.crossover_probability = 0.9;             // crossover probability
    
    params.variables_probability = 0.4;
    params.operators_probability = 0.5;
    params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !
    
    params.num_constants = 3; // use 3 constants from -1 ... +1 interval
    params.constants_min = -1;
    params.constants_max = 1;
    


}


bool read_graphs(t_graph *&training_graphs, int& num_training_graphs){
    
    // training is done on 4 graphs
    num_training_graphs = 0;
    
    training_graphs = new t_graph[4];
    
    char* file_name = new char[30];
    
    strcpy(file_name,"data//bayg29.tsp");
    
    if (!read_training_data_distances(training_graphs, num_training_graphs, file_name)) {
        printf("Cannot find input file(s)! Please specify the full path!");
        getchar();
        return 0;
    }

    
  
    
    strcpy(file_name,"data//a280.tsp");
    
    if (!read_training_data_coordinates(training_graphs, num_training_graphs, file_name)) {
        printf("Cannot find input file(s)! Please specify the full path!");
        getchar();
        return 0;
    }
    
    
    
    strcpy(file_name,"data//berlin52.tsp");
    
    if (!read_training_data_coordinates(training_graphs, num_training_graphs, file_name)) {
        printf("Cannot find input file(s)! Please specify the full path!");
        getchar();
        return 0;
    }
    
    
    
    strcpy(file_name,"data//bier127.tsp");
    
    if (!read_training_data_coordinates(training_graphs, num_training_graphs, file_name)) {
        printf("Cannot find input file(s)! Please specify the full path!");
        getchar();
        return 0;
    }
    
    delete [] file_name;
    
    return 1;
}

//--------------------------------------------------------------------

int main(void)
{
    t_parameters params;
    
    init_params(params);
    
    
   
    
	t_graph *training_graphs = NULL;
	int num_training_graphs = 0;
    
  
    
    if (!read_graphs(training_graphs, num_training_graphs)) return 1;
    
    
	compute_global_variables(training_graphs, num_training_graphs);


	int num_variables = 10;


	srand(1); // if you want to make a 

	printf("evolving...\n");

     int num_threads = 4;
    
	start_steady_state(params, training_graphs, num_training_graphs, num_variables, num_threads);

	delete_training_graphs(training_graphs, num_training_graphs);

	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------
