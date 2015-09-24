#include "dataparser.h"
#include "datalexer.h"
#include "list_all.h"

#include <stdlib.h>
#include <string.h>
#include <setjmp.h>

struct DataParser {
  DataLexer lexer;
  Token token, matched_token;
  int error_no;
  jmp_buf parser_failed;
  // Parser state
  PmError *e;
  int row, num_cols;
  ListStr col_names;
  DataType *col_datatype;
  ListStr *col_data;
};



/* ***  General parsing methods *** */
static TokenType next_token(DataParser parser);
static TokenType match_token(DataParser parser, TokenType type);
static TokenType get_tokenType(DataParser parser);
static char* get_matched_buffer(DataParser parser);
/* **** Data parsing specific methods **** */
static void clear_state(DataParser parser);
static void parse_datafile(DataParser parser);
static void parse_header(DataParser parser);
static void parse_rows(DataParser parser);


DataParser dataparser_new(FILE* fp) {
  DataParser parser;
  parser = (DataParser)malloc(sizeof(struct DataParser));
  parser->lexer = datalexer_new(fp);
  /* initialize state variables */
  parser->row = 0;
  parser->col_names = list_str_new();
  parser->num_cols = 0;
  parser->col_data = NULL;
  return parser;
}

void dataparser_free(DataParser* pparser) {
  DataParser parser = *pparser;
  int i;
  clear_state(parser);
  list_str_free(&parser->col_names);
  datalexer_free(&parser->lexer);
  // data columns
  free(parser); // the actual structure: struct PhyParser
  *pparser = NULL;
}

DataFrame* dataparser_parse(DataParser parser, PmError *e) {
    // initialize state
  DataFrame* data = NULL;
  clear_state(parser);
  parser->e = e;
  if (setjmp(parser->parser_failed)) {
    // PARSER error handling code
    clear_state(parser);
    return NULL;
  }
  parse_datafile(parser);
  /* create dataframe with parsed info */
  data = df_new_from_strings(parser->num_cols,parser->col_names,parser->col_datatype,parser->col_data); /* ojo remove last e arg */
  return data;
}

/* ***  General parsing methods *** */

/* Reads/consumes token from lexer and stores it on state */
/* Usually called after a succesful match. Previous matched token is logged */
static TokenType next_token(DataParser parser) {
  int status;
  TokenType ttype;
  parser->matched_token = parser->token;
  ttype = datalexer_gettoken(parser->lexer,&parser->token);
  //printf("next token %d %d",ttype,parser->token.type);
  //token_print(&parser->token,stdout);
  if (ttype == TOKEN_ERROR) {
    PmError *e = parser->e;
    Token* token = &parser->token;
    pm_error_set_type(e,PM_EPARSER);
    pm_error_set_message(e,"Lexical analyzer: Couldn't get next token\n");
    pm_error_append_message(e,"(line %d,column %d)",token->line,token->column);
    pm_error_append_message(e, " found: %s",token->buffer);
    // throw ParserExc(this->token,"Lexical analyzer error");
    longjmp(parser->parser_failed,1);
  }
  return ttype;
}

/* Tries to match token with token stored in parser */
/* If succesful, it consumes next token */
static TokenType match_token(DataParser parser, TokenType match_type) {
  TokenType ttype;
  ttype = token_get_type(&parser->token);
  if (ttype != match_type) {
    PmError *e = parser->e;
    Token* token = &parser->token;
    pm_error_set_type(e,PM_EPARSER);
    pm_error_set_message(e,"Couldn't match: %s\n",datalexer_typetochar(parser->lexer,match_type));
    pm_error_append_message(e,"(line %d,column %d)",token->line,token->column);
    pm_error_append_message(e, " found instead: ",token->buffer);
    pm_error_append_message(e,"%s(%s)\n",
			    token->token_set->token_name[ttype],token->buffer);
    // throw ParserExc(this->token,"Matching error");
    longjmp(parser->parser_failed,1);
  }
  //printf("MATCHED: "); token_print(&parser->token,stdout);
  ttype = next_token(parser);
  return ttype;
} 

static TokenType get_tokenType(DataParser parser) { 
  return token_get_type(&parser->token); 
}

static char* get_matched_buffer(DataParser parser) { 
  // return token_get_pbuffer(&parser->token);
  // problem with linker. Maybe because its inlined?
  return &parser->matched_token.buffer[0];
}


/* **************************************** */

static void clear_state(DataParser parser) {
  int i;
  ListStr l;
  char *name;

  /* printf("Clearing state\n"); */
  /* clear column names */
  l = parser->col_names;
  while(!list_str_empty(l)) {
    list_str_pop_front(l,&name);
    /* printf(" free: %s\n",name); */
    free(name);
  }
  /* clear column data */
  for(i=0; i < parser->num_cols; i++) {
    /* printf("datatype %s\n",df_datatype_name(parser->col_datatype[i])); */
    l = parser->col_data[i];
    if (l!=NULL) {
      while(!list_str_empty(l)) {
	list_str_pop_front(l,&name);
	/* printf(" free: %s\n",name); */
	free(name);
      }
      list_str_free(&parser->col_data[i]);
    }
  }
  if (parser->num_cols > 0) {
    free(parser->col_data);
    parser->col_data=NULL;
    parser->num_cols = 0;
    free(parser->col_datatype);
  }
  
  return;
}

static void parse_datafile(DataParser parser) {
  int i;
  parse_header(parser);
  /* we could, at his point, read all new lines in order to have an idea 
     of the number or rows */
  /* initilize data columns */
  parser->col_data = (ListStr*)malloc(sizeof(ListStr)*parser->num_cols);
  parser->col_datatype = (DataType*)malloc(sizeof(DataType)*parser->num_cols);
  for(i=0; i < parser->num_cols; i++) {
    parser->col_data[i] = list_str_new();
    parser->col_datatype[i] = DF_NOTYPE;
  }
  parse_rows(parser);

}


static void parse_header(DataParser parser) {
  char *col_name,*buffer;
  TokenType ttype;
  Token *token;
  int finish = 0;
  ttype = next_token(parser);
  while(ttype != TOKEN_EOL) {  
    ttype = match_token(parser,TOKEN_STRING);
    buffer = get_matched_buffer(parser);
    col_name=(char*)malloc(sizeof(char)*(strlen(buffer)+1));
    strcpy(col_name,buffer);
    list_str_push_back(parser->col_names,col_name);
    /* ttype = get_tokenType(parser); */
  }  
  match_token(parser,TOKEN_EOL);
  parser->row++;
  parser->num_cols = list_str_size(parser->col_names);
  return;

}



/* No empty rows, no NaN */
/* row size must match cnumber of columns */
static void parse_rows(DataParser parser) {
  TokenType ttype;
  DataType dtype, *col_datatype;
  int col_num, num_cols, too_many_columns;
  Token *token;
  char *buffer,*string;
  PmError *e;

  num_cols = parser->num_cols;
  col_datatype = parser->col_datatype;
  ttype = get_tokenType(parser);

  too_many_columns = 0;
  while (ttype != TOKEN_SCANEOF) { /* main loop */
    /* read row */
    col_num = 0;
    while((ttype != TOKEN_EOL)&&(ttype != TOKEN_SCANEOF)) { /* row loop */
      if (col_num >= num_cols) { too_many_columns = 1; break; }
      switch(ttype) {
      case TOKEN_INT:
	dtype = DF_INT;
	break;
      case TOKEN_FLOAT:
	dtype = DF_DOUBLE;
	break;
      case TOKEN_STRING:
	dtype = DF_STRING;
	break;
      default: /* error */
	e = parser->e;
	token = &parser->token;
	pm_error_set_type(e,PM_EPARSER);
	pm_error_set_message(e,"Invalid data value\n");
	pm_error_append_message(e,"(line %d,column %d)",token->line,token->column);
	pm_error_append_message(e, " found: %s",token->buffer);
	longjmp(parser->parser_failed,1);
      }
      if (dtype > col_datatype[col_num]) col_datatype[col_num] = dtype;
      ttype = match_token(parser,ttype); /* consume token */
      buffer = get_matched_buffer(parser);
      string = (char*)malloc(sizeof(char)*(strlen(buffer)+1));
      strcpy(string,buffer);
      list_str_push_back(parser->col_data[col_num],string);
      col_num++;
      /* ttype = get_tokenType(parser); */
    } /* end row loop */
    if ((col_num < num_cols)|| too_many_columns) {
      e = parser->e;
      pm_error_set_type(e,PM_EPARSER);
      pm_error_set_message(e,"Incorrect number of columns in row %d\n",parser->row);
      pm_error_append_message(e, "Data must match number of headers (%d). ",num_cols);
      if (too_many_columns) {
	pm_error_append_message(e,"Too many columns.");
      } else {
	pm_error_append_message(e, "Found %d column(s)",col_num);
      }
      longjmp(parser->parser_failed,1);
    }
    ttype = match_token(parser,ttype);
    parser->row++;
  } /* end main loop */
  

}
