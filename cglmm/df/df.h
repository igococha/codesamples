#ifndef DF_H
#define DF_H

#include "list_all.h"
#include "list_any.h"
#include "dense.h"


#include <stdio.h>
#include <stddef.h>

typedef enum { DF_NOTYPE, DF_INT, DF_DOUBLE, DF_STRING, DF_FACTOR } DataType;

typedef struct DataFrame DataFrame;

struct DataFrame {
  int num_cols;
  int num_rows;
  char **col_names;
  DataType* col_datatype;
  void **col_data;

};

DataFrame* df_load_new(const char * filename);
DataFrame* df_new_from_strings(int num_cols, ListStr col_names, DataType* col_datatype, ListStr *col_data);
const char* df_datatype_name(DataType datatype);
void df_columns_print(DataFrame* df, FILE* fp);
void df_print(DataFrame* df, FILE* fp);
void df_free(DataFrame** data);


DataFrame* df_load_columns(ListAny* col_names, const char *filename);
int df_drop_columns(DataFrame* data, ListAny* col_names);

int df_num_rows(DataFrame* df);
int df_get_column_idx(DataFrame* df, char* col_name);
void* df_get_column(DataFrame* df, char* col_name, DataType *type);
int df_permute(DataFrame* df,int* perm);

Dense* df_gen_matrix(DataFrame* df, char** col_names, int size, int ones,int prepend);
DenseVector* df_gen_vector(DataFrame* df, char *col_name,int prepend);

#endif


