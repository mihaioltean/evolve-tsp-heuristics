//
//  config.cpp
//  mep_evolve
//
//  Created by Virginia Niculescu on 23/02/16.
//  Copyright © 2016 Virginia Niculescu. All rights reserved.
//

#include "config.hpp"

char operators_string[num_operators][10] = { "+", "-", "*", "min", "max" };

//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params)
{
    c.prg = new t_code3[params.code_length];
    if (params.num_constants)
        c.constants = new double[params.num_constants];
    else
        c.constants = NULL;
}
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
void allocate_training_data(double **&data, double *&target, int num_training_data, int num_variables)
{
    target = new double[num_training_data];
    data = new double*[num_training_data];
    for (int i = 0; i < num_training_data; i++)
        data[i] = new double[num_variables];
}
//---------------------------------------------------------------------------
void allocate_partial_expression_values(double ***&expression_value, int num_training_data, int code_length, int num_threads)
{
    // partial values are stored in a matrix of size: code_length x num_training_data
    // for each thread we have a separate matrix
    expression_value = new double**[num_threads];
    for (int t = 0; t < num_threads; t++) {
        expression_value[t] = new double*[code_length];
        for (int i = 0; i < code_length; i++)
            expression_value[t][i] = new double[num_training_data];
    }
}
//---------------------------------------------------------------------------
void delete_partial_expression_values(double ***&expression_value, int code_length, int num_threads)
{
    if (expression_value) {
        for (int t = 0; t < num_threads; t++) {
            for (int i = 0; i < code_length; i++)
                delete[] expression_value[t][i];
            delete[] expression_value[t];
        }
        delete[] expression_value;
    }
}
//---------------------------------------------------------------------------
