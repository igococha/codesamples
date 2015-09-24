#ifndef DATAPARSER_H
#define DATAPARSER_H

#include "datalexer.h"
#include "df.h"
#include "pm_error.h"
#include <stdio.h>

// opaque pointer
typedef struct DataParser *DataParser;

DataParser dataparser_new(FILE* fp);
void dataparser_free(DataParser* pparser);


/* ***  Tree specific methods *** */
DataFrame* dataparser_parse(DataParser parser, PmError *e);

#endif
