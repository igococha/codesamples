#ifndef CORE_LEXER_H
#define CORE_LEXER_H

#include <stdio.h>

/*
Basic structures used for the implementation of simple parsers.
 */

/* ****** Token declaration ***** */
#define TOKEN_BS 128
#define TOKEN_NOTHING 0
#define TOKEN_ERROR 1
#define TOKEN_SCANEOF 2


// set of valid token types
// token_num >= 2
// token_name[0] = "NOTHING", token_name[1] = "ERROR"
typedef struct TokenSet {
  int token_num;
  char* token_name[];
} TokenSet;

typedef int TokenType;
typedef struct Token {
  TokenType type;
  char buffer[TOKEN_BS];
  unsigned short line,column;
  const TokenSet *token_set;
} Token;

const char* tokenset_typetochar(const TokenSet* ts, TokenType type);

void token_init(Token* token, const TokenSet* ts);
void token_print(Token* token, FILE* fp);
TokenType token_get_type(Token* token);
void token_set_type(Token* token, TokenType type);
//char* token_get_pbuffer(Token* token);
void token_set(Token* token, TokenType t, unsigned short l, unsigned short c);

/* ****** Lexer declaration ***** */

// typedef struct Lexer *Lexer;
#define LEXER_BS 4

typedef struct Lexer {
  FILE* fp;
  char buffer1[LEXER_BS+1];
  char buffer2[LEXER_BS+1];
  char *buffer;
  char *last1, *last2, *second1, *second2, *last, *second, *peof;
  char *start;
  unsigned short line,column;
  int error_num, is_eof;
  TokenType last_token;
} *Lexer;

Lexer lexer_new(FILE* fp);
//int lexer_gettoken(Lexer lexer, Token* token);
void lexer_free(Lexer* plexer);
char* load_buffer(Lexer lexer);

#endif
