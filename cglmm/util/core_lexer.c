#include "core_lexer.h"
#include <stdio.h>
#include <stdlib.h>



const TokenSet default_ts = { 3, { "NOTHING", "ERROR", "SCANEOF"} };

void token_init(Token* token, const TokenSet* token_set) {
  token->type = TOKEN_NOTHING;
  token->buffer[0] = '\0';
  token->line = token->column = 0;
  token->token_set = (token_set==NULL)? &default_ts: token_set;
}

void token_print(Token* token, FILE* fp) {
  TokenType tt = token->type;
  if (tt < token->token_set->token_num) {
      fprintf(fp,"%s(%s)\n",token->token_set->token_name[tt],token->buffer);
  } else {
    fprintf(fp,"Unknown Token\n");
    //exit(0);
  }
  return;
}

inline TokenType token_get_type(Token* token) {
  return token->type;
}

void token_set_type(Token* token, TokenType type) {
  token->type = type;
}

inline char* get_pbuffer(Token* token) { 
  return &token->buffer[0]; 
}

void token_set(Token* token, TokenType t, 
	       unsigned short l, unsigned short c) {
  token->type = t;
  token->line = l;
  token->column = c;
}

inline const char* tokenset_typetochar(const TokenSet* token_set, 
					TokenType type) {
  return token_set->token_name[type];
}

/* ************ LEXER ******************** */

Lexer lexer_new(FILE* fp) {
  Lexer lexer;
  int num_read;
  lexer = (Lexer)malloc(sizeof(struct Lexer));
  lexer->fp = fp;
  lexer->second1 = &lexer->buffer1[0] + 1;
  lexer->second2 = &lexer->buffer2[0] + 1;
  lexer->last1 = lexer->second1 + LEXER_BS - 1; // end of first buffer
  lexer->last2 = lexer->second2 + LEXER_BS - 1; // end of second buffer
  /* load into buffer */
  lexer->second = lexer->second1;
  lexer->last = lexer->last1;
  num_read = fread(lexer->second,sizeof(char), LEXER_BS,fp);
  lexer->start = lexer->second;
  if (num_read < LEXER_BS) {
    lexer->peof = lexer->second + num_read;
  } else {
    lexer->peof = NULL;
  }
  lexer->is_eof = 0;
  lexer->last_token = TOKEN_NOTHING;
  lexer->line = 1;
  lexer->column=0;
  lexer->error_num = 0;
  return lexer;
}

// condition: forward == last
char* load_buffer(Lexer lexer) {
  int num_read;
  char c;
  c = *lexer->last;
  if (lexer->second == lexer->second1) { // load to buffer2
    lexer->second = lexer->second1;
    lexer->last = lexer->last1;
  } else { // load to buffer2
    lexer->second = lexer->second2;
    lexer->last = lexer->last2;
  }
  *(lexer->second - 1) = c;
  num_read = fread(lexer->second,sizeof(char), LEXER_BS,lexer->fp);
  if (num_read < LEXER_BS) {
    lexer->peof = lexer->second + num_read;
  } else {
    lexer->peof = NULL;
  }
  return lexer->second;
}

void lexer_free(Lexer* plexer) {
  Lexer lexer = *plexer;
  //fclose(lexer->fp);
  free(lexer);
  *plexer = NULL;
  return;
}
