//
//  I_operations.cpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//

#include "I_operations.hpp"



//---------------------------------------------------------------------------

bool read_training_data_distances(t_graph *training_graphs, int &num_training_graphs, char* file_name)
{
 
   
    
    int k = num_training_graphs; // count the graphs

        FILE* f = fopen(file_name, "r");
        if (!f)
            return false;
        
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
    
    num_training_graphs++;
    return true;
 
}
//---------------------------------------------------------------------------

bool read_training_data_coordinates(t_graph *training_graphs, int &num_training_graphs, char* file_name)
{
     int k = num_training_graphs; // count the graphs
  
        FILE* f = fopen(file_name, "r");
        if (!f)
            return false;
        
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
    num_training_graphs++;
    return true;
    
}

//---------------------------------------------------------------------------
void delete_training_graphs(t_graph *&training_graphs, int num_training_graphs)
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
