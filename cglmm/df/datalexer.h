#ifndef DATALEXER_H
#define DATALEXER_H

#include "core_lexer.h"
#include "df.h"
#include <stdio.h>

/* ****** Tokens ******* */

// already defined TOKEN_NOTHING, TOKEN_ERROR, TOKEN_SCANEOF
// start from 3
#define TOKEN_FLOAT 3
#define TOKEN_INT 4
#define TOKEN_STRING 5
#define TOKEN_COLON 6
#define TOKEN_COMMA 7
#define TOKEN_SPACE 8
#define TOKEN_EOL 9

/* ********* Lexer ********** */

// opaque pointer
typedef struct DataLexer *DataLexer;

DataLexer datalexer_new(FILE* fp);
void datalexer_free(DataLexer* lexer);
TokenType datalexer_gettoken(DataLexer lexer, Token* token);
const char* datalexer_typetochar(DataLexer lexer, TokenType type);


#endif
