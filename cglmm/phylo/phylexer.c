#include "phylexer.h"
#include "core_lexer.h"
#include "pm_mem.h"

#include <stdlib.h>
#include <ctype.h>

const TokenSet phy_ts = { 11, { "NOTHING", "ERROR", "SCANEOF",
				  "FLOAT", "INT", "ID", "LPAREN", "RPAREN",
				  " SEMICOLON", "COLON", "COMMA"} };

struct PhyLexer {
  Lexer lexer;
  PhyFormat format;
  const TokenSet *token_set;
};


PhyLexer phylexer_new(FILE* fp, PhyFormat format) {
  PhyLexer phylexer;
  phylexer = (PhyLexer) PM_MEM_ALLOC(sizeof(struct PhyLexer));
  if (phylexer==NULL) return NULL;
  phylexer->lexer = lexer_new(fp);
  phylexer->format = format;
  phylexer->token_set = &phy_ts;
  return phylexer;
}


void phylexer_free(PhyLexer* plexer) {
  lexer_free(&(*plexer)->lexer);
  PM_MEM_FREE(*plexer);
}


int phylexer_gettoken(PhyLexer phylexer, Token* token) {
  unsigned short state;
  unsigned int line, column, column_start;
  int backtrack,row,col;
  char c, *forward, *pbuffer;
  TokenType matched;
  Lexer lexer = phylexer->lexer;
  token_init(token,phylexer->token_set);
  if (lexer->is_eof) {
    token->type = TOKEN_SCANEOF;
    return TOKEN_SCANEOF;
  }
  if (lexer->last_token==TOKEN_ERROR) {
    token->type = TOKEN_ERROR;
    return TOKEN_ERROR;
  }
  line = lexer->line;
  column = lexer->column;
  forward = lexer->start;
  state = 0;
  backtrack = 0;
  pbuffer = &token->buffer[0];
  matched = TOKEN_NOTHING;
  while(matched == TOKEN_NOTHING) {
    if (forward == lexer->peof) {
      c = '\0'; // denotes EOF character
      lexer->is_eof = 1;
    } else {
      c = *forward;
      if (forward == lexer->last) {
	forward = load_buffer(lexer);
    } else
      forward++;
    }
    column++;
    if (state==0) {
	column_start = column;
    }
    switch(state) {
    case 0: /* consume whitespace and point to next state given first char */
      switch (c) {
      case ' ':
      case '\t':
	break;
      case '\n': /* missing CR+LF or LF+CR -check df lexer for example  */
      case '\r':
	line++;
	column = 0;
	break;
      case '(':
	matched = TOKEN_LPAREN; break;
      case ')':
	matched = TOKEN_RPAREN; break;
      case ';':
	matched = TOKEN_SEMICOLON; break;
      case ':':
	matched = TOKEN_COLON; break;
      case ',':
	matched = TOKEN_COMMA; break;
      case '\0':
	matched = TOKEN_SCANEOF; break;
      default:
	if (isalpha(c)) {
	  *(pbuffer++) = c;
	  state = 1;  // match ID
	} else if (isdigit(c)) {
	  *(pbuffer++) = c;
	  state = 2;  // match NUMBER ::= INT | FLOAT
	} else {
	  *(pbuffer++) = c;
	  matched = TOKEN_ERROR;
	}
      } // switch(c)
      break;
    case 1: // ID
      if (isalpha(c) || (isdigit(c)) || (c=='_') || (c=='_')) {
	state = 1; *(pbuffer++) = c;
      } else { // end of ID - backtrack
	matched = TOKEN_ID; backtrack = 1;
      }
      break;
    case 2: // INT
      if (isdigit(c)) {
	state = 2; *(pbuffer++) = c;
      } else if (c == '.') {
	state = 3; *(pbuffer++) = c;
      } else {
	matched = TOKEN_FLOAT; backtrack = 1;
      }
      break;
    case 3: // FLOAT
      if (isdigit(c)) {
	state = 3; *(pbuffer++) = c;
      } else {
	matched = TOKEN_FLOAT; backtrack = 1;
      }
      break;
    default:
      printf("Error: unknown state\n"); // it should terminate actually
      token->type = TOKEN_ERROR;
      return TOKEN_ERROR;
    } // switch (state)
  } // while
  // Bookkeeping after matching
  token_set(token,matched,line,column_start);
  *pbuffer = '\0';
  lexer->line = line;
  if (backtrack>0) {
    lexer->start = forward - backtrack;
    lexer->column = column - backtrack;
    // tested for backtrack=1 - may need to backtrack pbuffer as well
  } else {
    lexer->start = forward;
    lexer->column = column;
  }
  lexer->last_token = token->type;
  return matched;
}

const char* phylexer_typetochar(PhyLexer lexer, TokenType type) {
  return tokenset_typetochar(lexer->token_set,type);
}
