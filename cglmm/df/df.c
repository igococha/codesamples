#include "df.h"
#include "dataparser.h"

#include <stdlib.h>
#include <string.h>

/* local interface */
static DataFrame* df_new();

/* end local interface */

static DataFrame* df_new() {
  DataFrame *df;
  df = (DataFrame*)malloc(sizeof(DataFrame));
  df->num_cols = 0;
  df->col_names = NULL;
  df->col_datatype = NULL;
  
  return df;
}


/* *** Implementation global interface *** */

DataFrame* df_load_new(const char * filename) {
  FILE* fp;
  DataParser dataparser;
  DataFrame *data;

  pm_error_clear(NULL);
  fp = fopen(filename,"r");
  if (fp == NULL) {
    PM_ERROR_VAL(PM_EFILE,NULL,"Couldn't open data file: %s",filename);
  }
  dataparser = dataparser_new(fp);
  if (dataparser != NULL) {
    data = dataparser_parse(dataparser,pm_error_default());
  } else {
    fclose(fp);
    PM_ERROR_VAL(PM_EPARSER,NULL,"Couldn't create parser: %s",filename);
  }
  fclose(fp);
  dataparser_free(&dataparser);

  return data;

}

DataFrame* df_new_from_strings(int num_cols, ListStr col_names, DataType* col_datatype, ListStr *col_data) {
  int i,j,*pint,num_rows;
  double *pdouble;
  char *name, **pstring;
  DataFrame *df;
  DataType dtype;
  ListStr l;
  
  df = df_new();
  df->num_cols = num_cols;
  num_rows = df->num_rows = list_str_size(col_data[i]);
  df->col_names = (char**)malloc(sizeof(char*)*num_cols);
  df->col_datatype = (DataType*)malloc(sizeof(DataType)*num_cols);
  df->col_data = (void**)malloc(sizeof(void*)*num_cols);
  for(i=0; i < num_cols; i++) {
    dtype = df->col_datatype[i] = col_datatype[i]; /* set datatype */
    list_str_pop_front(col_names, &df->col_names[i]);
    l = col_data[i];
    switch(dtype) {
    case DF_INT:
      j=0;
      pint = (int*)malloc(sizeof(int)*num_rows);
      while(!list_str_empty(l)) {
	list_str_pop_front(l,&name);
	pint[j++] = atoi(name);
	free(name);
      }
      df->col_data[i] = (void*)pint;
      break;
    case DF_DOUBLE:
      j=0;
      pdouble = (double*)malloc(sizeof(double)*num_rows);
      while(!list_str_empty(l)) {
	list_str_pop_front(l,&name);
	pdouble[j++] = atof(name);
	free(name);
      }
      df->col_data[i] =(void*)pdouble;
      break;
    case DF_STRING:
      j=0;
      pstring = list_str_clear_toarray(l);
      df->col_data[i] =(void*)pstring;
      break;
    default:
      df_free(&df);
      PM_ERROR_VAL_NOMSG(PM_EARGUMENTS,NULL);
      return NULL;
    }
  }
  /* col_names list should be empty */
  /* all lists in col_data should be empty */
  return df;
}

const char* df_datatype_name(DataType datatype) {
  switch(datatype) {
  case DF_NOTYPE:
    return "NOTYPE";
  case DF_INT:
    return "INT";
  case DF_DOUBLE:
    return "DOUBLE";
  case DF_STRING:
    return "STRING";
  case DF_FACTOR:
    return "FACTOR";
  }
  return "UNKNOWN"; /* should be an error */
}

void df_columns_print(DataFrame* df, FILE* fp) {
  int i,j;
  fprintf(fp,"Columns:\n");
  for(i=0; i < df->num_cols; i++) {
    fprintf(fp,"%s(%s)\n",df->col_names[i],df_datatype_name(df->col_datatype[i]));
  }
}

void df_print(DataFrame* df, FILE* fp) {
  int i, j,*pint;
  double *pdouble;
  char **pstring;
  fprintf(fp,"Datafarame %d by %d\n",df->num_rows,df->num_cols);
  for(i=0; i < df->num_cols; i++) {
    fprintf(fp,"%s(%s)\n",df->col_names[i],df_datatype_name(df->col_datatype[i]));
    switch(df->col_datatype[i]) {
    case DF_INT:
      pint = (int*)df->col_data[i];
      for(j=0; j < df->num_rows; j++,pint++) {
	fprintf(fp,"%d ",*pint);
      }
      fprintf(fp,"\n"); break;
    case DF_DOUBLE:
      pdouble = (double*)df->col_data[i];
      for(j=0; j < df->num_rows; j++,pdouble++) {
	fprintf(fp,"%f ",*pdouble);
      }
      fprintf(fp,"\n"); break;
    case DF_STRING:
      pstring = (char**)df->col_data[i];
      for(j=0; j < df->num_rows; j++,pstring++) {
	fprintf(fp,"%s ",*pstring);
      }
      fprintf(fp,"\n"); break;
    default:
      printf("Incorrect data type\n");
      return;
    }
  }
  return;
}

void df_free(DataFrame** data) {
  DataFrame *df = *data;
  int i,j;
  char *str,**pstring;
  
  if (df == NULL) return;
  /* free data */
  for(i=0;i<df->num_cols;i++) {
    if (df->col_datatype[i]==DF_STRING) {
      pstring = (char**)df->col_data[i];
      for(j=0;j<df->num_rows;j++,pstring++) {
	/* printf("freeing %s\n",*pstring); */
	free(*pstring);
      }
    }
    /* printf("freeing %s\n",df->col_names[i]); */
    free(df->col_data[i]);
  }
  free(df->col_data);
  /* free col names and array storing them */
  for(i=0; i < df->num_cols; i++) {
    str = df->col_names[i];
    if (str) free(str);
  }
  free(df->col_names);
  free(df->col_datatype);
  /* free structure */
  free(df);
  *data = NULL;
}


DataFrame* df_load_columns(ListAny* col_names, const char *filename) {
  
  return NULL;
}

int df_drop_columns(DataFrame* data, ListAny* col_names) {
  return 0;
}

int df_num_rows(DataFrame* data) {
  return data->num_rows;
}

int df_get_column_idx(DataFrame* df, char* col_name) {
  int i;
  for(i=0; i < df->num_cols; i++) {
    if (strcmp(col_name,df->col_names[i])==0)
      break;
  }
  if (i == df->num_cols) {
    return -1;
  } else {
    return i;
  }
}


/* Search for col_name in dataframe
 * returns pointer to column array, NULL otherwise
 * and sets datatype to corresponding value (DF_NOTYPE if falure)
 * NOTE: No allocation is made!
 */
void* df_get_column(DataFrame* df, char* col_name, DataType *type) {
  PmError *e;
  int i;

  i = df_get_column_idx(df,col_name);
  if (i < 0) {
    *type = DF_NOTYPE;
    return NULL;
  }
  *type = df->col_datatype[i];
  return df->col_data[i];
}

int df_permute(DataFrame* df,int* perm) {
  double *new_double, *pdouble;
  int *new_int, *pint, num_rows;
  char **new_string, **pstring;
  void* temp_void;
  
  int i,j;
  new_double = NULL;
  new_int = NULL;
  new_string = NULL;
  num_rows = df->num_rows;

  for(i=0; i < df->num_cols; i++) {
    temp_void = df->col_data[i];
    switch(df->col_datatype[i]) {
    case DF_INT:
      pint = (int*)temp_void;
      if (new_int==NULL) new_int=(int*)malloc(sizeof(int)*num_rows);
      for(j=0; j < num_rows; j++,pint++) {
	new_int[perm[j]] = *pint;
      }
      df->col_data[i] = (void*)new_int;
      new_int = (int*)temp_void;
      break;
    case DF_DOUBLE:
      pdouble = (double*)temp_void;
      if (new_double==NULL) new_double=(double*)malloc(sizeof(double)*num_rows);
      for(j=0; j < num_rows; j++,pdouble++) {
	new_double[perm[j]] = *pdouble;
      }
      df->col_data[i] = (void*)new_double;
      new_double = (double*)temp_void;
      break;
    case DF_STRING:
      pstring = (char**)temp_void;
      if (new_string==NULL) new_string=(char**)malloc(sizeof(char*)*num_rows);
      for(j=0; j < num_rows; j++,pstring++) {
	new_string[perm[j]] = *pstring;
      }
      df->col_data[i] = (void*)new_string;
      new_string = (char**)temp_void;
      break;
    default:
      printf("Incorrect data type\n");
    }
  }
  if (new_double) free(new_double);
  if (new_int) free(new_int);
  if (new_string) free(new_string);
  return PM_SUCCESS;
}



/* Generate Dense matrix from column values*/
/* Prepend: number of 0 rows to be added before the df data */
Dense* df_gen_matrix(DataFrame* df, char** col_names, int size, int ones, int prepend) {
  Dense* m=NULL;
  int i, *col_index, idx,num_rows;
  DataType type;
  void* col_data;

  pm_error_clear(NULL);    
  ones = (ones!=0)?1:0;
  col_index = (int*)malloc(sizeof(int)*size);
  for(i=0; i < size; i++) {
    idx = df_get_column_idx(df,col_names[i]);
    if (idx < 0) {
      free(col_index);
      PM_ERROR_VAL(PM_EDF,NULL,"Couldn't generate matrix. Column %s not in dataframe",col_names[i]);
    } else {
      type = df->col_datatype[idx];
      if ((type!=DF_INT)&&(type!=DF_DOUBLE)) {
	free(col_index);
	PM_ERROR_VAL(PM_EDF,NULL,"Couldn't generate matrix. Column %s not of numeric type",col_names[i]);
      }
      col_index[i] = idx;
    }
  }
  printf("Generation X\n");
  num_rows = df->num_rows+prepend;
  m = dense_new(num_rows,ones+size);
  if (prepend>0) {
    dense_set_all(m,0.0);
  }

  if (ones) {
    dense_set_column_range(m,0,1.0,prepend,num_rows-1);
  }
  for(i=0; i < size; i++) { /* for each column */
    printf("%s\n",col_names[i]);
    idx = col_index[i];
    type = df->col_datatype[idx];
    if (type==DF_INT) {
      dense_set_column_intarray_range(m,i+ones,(int*)df->col_data[idx],
				prepend,num_rows-1);
    } else { /* DF_DOUBLE */
      dense_set_column_doublearray_range(m,i+ones,(double*)df->col_data[idx],
				   prepend,num_rows-1);
    }
  }

  free(col_index);
  return m;
  
}

DenseVector* df_gen_vector(DataFrame* df, char *col_name,int prepend) {
  DenseVector *v=NULL;
  void *data;
  DataType type;
  int num_rows;

  pm_error_clear(NULL);    
  data = df_get_column(df,col_name,&type);
  if (data==NULL) {
    PM_ERROR_VAL(PM_EDF,NULL,"Couldn't generate vector. Column %s not in dataframe",col_name);
  }
  if ((type!=DF_INT)&&(type!=DF_DOUBLE)) {
    PM_ERROR_VAL(PM_EDF,NULL,"Couldn't generate vector. Column %s not of numeric type",col_name);
  }
  num_rows = df->num_rows+prepend;
  v = dense_vector_new(num_rows);
  if (prepend > 0) {
    dense_vector_set_all(v,0.0);
  }
  if (type==DF_INT) {
    dense_vector_set_intarray_range(v,(int*)data,prepend,num_rows-1);
  } else { /* it's DOUBLE */
    dense_vector_set_doublearray_range(v,(double*)data,prepend,num_rows-1);
  }
  
  return v;
}
