#ifndef PHYLEXER_H
#define PHYLEXER_H

#include "core_lexer.h"
#include "phytree.h"
#include <stdio.h>

/* ****** Tokens ******* */

// already defined TOKEN_NOTHING, TOKEN_ERROR, TOKEN_SCANEOF
// start from 3
#define TOKEN_FLOAT 3
#define TOKEN_INT 4
#define TOKEN_ID 5
#define TOKEN_LPAREN 6
#define TOKEN_RPAREN 7
#define TOKEN_SEMICOLON 8
#define TOKEN_COLON 9
#define TOKEN_COMMA 10


/* ********* Lexer ********** */

// opaque pointer
typedef struct PhyLexer *PhyLexer;

PhyLexer phylexer_new(FILE* fp, PhyFormat format);
void phylexer_free(PhyLexer* plexer);
TokenType phylexer_gettoken(PhyLexer lexer, Token* token);
const char* phylexer_typetochar(PhyLexer lexer, TokenType type);


#endif
