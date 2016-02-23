//
//  I_operations.hpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//
// reading and deleting graphs


#ifndef I_operations_hpp
#define I_operations_hpp


#include "config.hpp"

//---------------------------------------------------------------------------
bool read_training_data_distances(t_graph *training_graphs, int &num_training_graphs, char* name);

//---------------------------------------------------------------------------
bool read_training_data_coordinates(t_graph *training_graphs, int &num_training_graphs, char* name);



//---------------------------------------------------------------------------
void delete_training_graphs(t_graph *&training_graphs, int num_training_graphs);
//---------------------------------------------------------------------------

#endif /* I_operations_hpp */
