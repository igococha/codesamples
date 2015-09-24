#ifndef PHYANIMAL_H
#define PHYANIMAL_H

#include "phytree.h"
#include "df.h"


typedef struct PhyAnimal *PhyAnimal;

/* Returns new Animal model, NULL if there was an error */
PhyAnimal phyanimal_new(char** y, int num_y, char** x, int num_x);
Dense* phyanimal_X(PhyAnimal model);
DenseVector* phyanimal_y(PhyAnimal model);
Dense* phyanimal_V(PhyAnimal model);

void phyanimal_print(PhyAnimal model,FILE* fp);
void phyanimal_free(PhyAnimal *model);
int phyanimal_load_data(PhyAnimal model,DataFrame *df, char *taxa_column, PhyTree *tree, int all_nodes);
void phyanimal_fit(PhyAnimal model, int comp_lambda);
double phyanimal_loglh(PhyAnimal model);

/* MCMC interface */

typedef struct mc2_animal *mc2_animal;

mc2_animal mc2_animal_new(PhyAnimal m);
void mc2_animal_free(mc2_animal *m);

int mc2_animal_run(mc2_animal m,unsigned int num_samples, unsigned int burnin, unsigned int thin);

#endif
