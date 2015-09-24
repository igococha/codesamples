#include "datalexer.h"
#include "core_lexer.h"
#include <stdlib.h>
#include <ctype.h>

static const TokenSet data_ts = { 10,
	{ "NOTHING", "ERROR", "SCANEOF","FLOAT","INT", "STRING",
	  "COLON","COMMA", "SPACE", "EOL"} };

struct DataLexer {
  Lexer lexer;
  const TokenSet *token_set;
};

DataLexer datalexer_new(FILE* fp) {
  DataLexer dlexer = (DataLexer)malloc(sizeof(struct DataLexer));
  dlexer->lexer = lexer_new(fp);
  dlexer->token_set = &data_ts;
  return dlexer;
}


void datalexer_free(DataLexer* dlexer) {
  lexer_free(&(*dlexer)->lexer);
  free(*dlexer);
  *dlexer = NULL;
}

/* this version considers spaces as separators */
/* does not produce SPACE token  */
/* TODO: a number state can change to a string state  
 * more freedom is required for strings (column-name/data)
*/

int datalexer_gettoken(DataLexer datalexer, Token* token) {
  unsigned short state;
  unsigned int line, column, column_start;
  int backtrack,row,col;
  char c, *forward, *pbuffer;
  TokenType matched;
  Lexer lexer = datalexer->lexer;
  token_init(token,datalexer->token_set);
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
      case '\n': /* <LF> 10 */
	state = 10; break;
      case '\r': /* <CR> 13 followed by 10 */
	state = 13; break;
      case ':':
	matched = TOKEN_COLON; break;
      case ',':
	matched = TOKEN_COMMA; break;
      case '\0':
	matched = TOKEN_SCANEOF; break;
      case '-':
	state = 4;
	*(pbuffer++) = c; break; /* negative number */
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
	matched = TOKEN_STRING; backtrack = 1;
      }
      break;
    case 2: // INT
      if (isdigit(c)) {
	state = 2; *(pbuffer++) = c;
      } else if (c == '.') {
	state = 3; *(pbuffer++) = c;
      } else {
	matched = TOKEN_INT; backtrack = 1;
      }
      break;
    case 3: // FLOAT
      if (isdigit(c)) {
	state = 3; *(pbuffer++) = c;
      } else {
	matched = TOKEN_FLOAT; backtrack = 1;
      }
      break;
    case 4: /* Read '-', Negative number*/
      if (isdigit(c)) {
	state = 2; *(pbuffer++) = c;
      } else {
	*(pbuffer++) = c;
	matched = TOKEN_ERROR;
      }
      break;
    case 10: /* EOL read '\n' */
      if (c != '\r') {
	backtrack = 1;
      }
      line++; column = 0;  matched = TOKEN_EOL;
    case 13: /* EOL read '\r' */
      if (c != '\n') { /* <LF> 10 */
	backtrack = 1;
      }
      line++; column = 0;  matched = TOKEN_EOL;
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

const char* datalexer_typetochar(DataLexer lexer, TokenType type) {
  return tokenset_typetochar(lexer->token_set,type);
}
