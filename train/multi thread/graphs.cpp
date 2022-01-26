#include <stdio.h>
#include <math.h>
#include <string.h>

#include "graphs.h"
//---------------------------------------------------------------------------
bool read_training_data(const char* path, t_graph*& training_graphs, int& num_training_graphs)
{
	// training is done on 4 graphs
	num_training_graphs = 4;
	training_graphs = new t_graph[num_training_graphs];

	int k = 0; // count the graphs
	{
		char* file_name = new char[strlen(path) + 100];
		strcpy(file_name, path);
		strcat(file_name, "bayg29.tsp");

		FILE* f = fopen(file_name, "r");
		delete[] file_name;
		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double* [training_graphs[k].num_nodes];
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
		char* file_name = new char[strlen(path) + 100];
		strcpy(file_name, path);
		strcat(file_name, "a280.tsp");

		FILE* f = fopen(file_name, "r");
		delete[] file_name;

		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double* [training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double* x = new double[training_graphs[k].num_nodes];
		double* y = new double[training_graphs[k].num_nodes];
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
		char* file_name = new char[strlen(path) + 100];
		strcpy(file_name, path);
		strcat(file_name, "berlin52.tsp");

		FILE* f = fopen(file_name, "r");
		delete[] file_name;

		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double* [training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double* x = new double[training_graphs[k].num_nodes];
		double* y = new double[training_graphs[k].num_nodes];
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
		char* file_name = new char[strlen(path) + 100];
		strcpy(file_name, path);
		strcat(file_name, "bier127.tsp");

		FILE* f = fopen(file_name, "r");
		delete[] file_name;

		if (!f)
			return false;

		fscanf(f, "%d", &training_graphs[k].num_nodes);
		// allocate the memory first
		training_graphs[k].distance = new double* [training_graphs[k].num_nodes];
		for (int i = 0; i < training_graphs[k].num_nodes; i++)
			training_graphs[k].distance[i] = new double[training_graphs[k].num_nodes];

		double* x = new double[training_graphs[k].num_nodes];
		double* y = new double[training_graphs[k].num_nodes];
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
	}	return true;
}
//---------------------------------------------------------------------------
void delete_training_graphs(t_graph*& training_graphs, int num_training_graphs)
{
	if (training_graphs)
		for (int i = 0; i < num_training_graphs; i++) {
			for (int j = 0; j < training_graphs[i].num_nodes; j++)
				delete[] training_graphs[i].distance[j];
			delete[] training_graphs[i].distance;
		}
	delete[] training_graphs;
}
//---------------------------------------------------------------------------
void compute_global_variables(t_graph* training_graphs, int num_training_graphs)
{
	// here we compute the min, max and average distance in a given graph
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
		training_graphs[k].average_distance /= (training_graphs[k].num_nodes * training_graphs[k].num_nodes) - training_graphs[k].num_nodes;
	}
}
//--------------------------------------------------------------------
