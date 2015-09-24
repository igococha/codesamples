#ifndef PHYPARSER_H
#define PHYPARSER_H

#include "phylexer.h"
#include "phytree.h"
#include "pm_error.h"
#include <stdio.h>

// opaque pointer
typedef struct PhyParser *PhyParser;

PhyParser phyparser_new(FILE* fp, PhyFormat format);
// PhyTrees *phyparser_parse(PhyParser parser);
void phyparser_free(PhyParser* pparser);


/* ***  Tree specific methods *** */
PhyTreeList* phyparser_parse(PhyParser parser, PmError *e);

#endif
