#include "phyparser.h"
#include "phylexer.h"
#include "list_any.h"
#include "pm_mem.h"
#include <stdlib.h>
#include <setjmp.h>

struct PhyParser {
  PhyLexer lexer;
  Token token, matched_token;
  int error_no;
  jmp_buf parser_failed;
  // Parser state
  PhyNode* root;
  ListAny tree_list; // PhyTree*
  ListAny internal_list; // PhyNode*
  ListAny tip_list; // PhyNode*
  ListAny children_stack;
  NodeId tip_counter, internal_counter;
  NodeId node_counter;    
  PmError *e;
};


/* ***  General parsing methods *** */
static TokenType next_token(PhyParser parser);
static void match_token(PhyParser parser, TokenType type);
static TokenType get_tokenType(PhyParser parser);
static char* get_matched_buffer(PhyParser parser);
/* **** Tree specific methods **** */
static void parse_trees(PhyParser parser);
static PhyTree* parse_tree(PhyParser parser);
static void parse_nodelist(PhyParser parser,PhyNode* anc, ListAny nodes);
static PhyNode* parse_node(PhyParser parser, PhyNode* anc);
static PhyNode* parse_internal(PhyParser parser, PhyNode* anc, int distance_must);
static PhyNode* parse_tip(PhyParser parser, PhyNode* anc);
static void clear_state(PhyParser parser);

PhyParser phyparser_new(FILE* fp,PhyFormat format) {
  PhyParser parser;
  parser = (PhyParser)PM_MEM_ALLOC(sizeof(struct PhyParser));
  if (parser==NULL) return NULL;
  parser->lexer = phylexer_new(fp,format);
  parser->tree_list = list_any_new();
  parser->internal_list = list_any_new();
  parser->tip_list = list_any_new();
  parser->children_stack = list_any_new();
  return parser;
}

void phyparser_free(PhyParser* pparser) {
  PhyParser parser = *pparser;
  clear_state(parser);
  list_any_free(&parser->children_stack);
  list_any_free(&parser->internal_list);
  list_any_free(&parser->tip_list);
  list_any_free(&parser->tree_list);
  phylexer_free(&parser->lexer);
  PM_MEM_FREE(*pparser); // the actual structure: struct PhyParser
}

PhyTreeList* phyparser_parse(PhyParser parser, PmError *e) {
    // initialize state
  PhyTreeList* tree_list = NULL;

  clear_state(parser);
  parser->e = e;
  parser->node_counter = 0;
  parser->tip_counter = parser->internal_counter = 0;
  if (setjmp(parser->parser_failed)) {
    // PARSER error handling code
    clear_state(parser);
    return NULL;
  }
  parse_trees(parser);
  tree_list = phytreelist_new(parser->tree_list);
  return tree_list;
}


/* ***  General parsing methods *** */

/* Reads/consumes token from lexer and sotores it on state */
static TokenType next_token(PhyParser parser) {
  int status;
  TokenType ttype;
  parser->matched_token = parser->token;
  ttype = phylexer_gettoken(parser->lexer,&parser->token);
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

static void match_token(PhyParser parser, TokenType match_type) {
  TokenType ttype;
  PmError *e;
  Token *token;
  ttype = token_get_type(&parser->token);
  if (ttype != match_type) {
    e = parser->e;
    token = &parser->token;
    pm_error_set_type(e,PM_EPARSER);
    pm_error_set_message(e,"Couldn't match: %s\n",phylexer_typetochar(parser->lexer,match_type));
    pm_error_append_message(e,"(line %d,column %d)",token->line,token->column);
    pm_error_append_message(e, " found instead: ",token->buffer);
    pm_error_append_message(e,"%s(%s)\n",
			    token->token_set->token_name[ttype],token->buffer);
    // throw ParserExc(this->token,"Matching error");
    longjmp(parser->parser_failed,1);
  }
  //printf("MATCHED: "); token_print(&parser->token,stdout);
  ttype = next_token(parser);
  return;
} 

static TokenType get_tokenType(PhyParser parser) { 
  return token_get_type(&parser->token); 
}

static char* get_matched_buffer(PhyParser parser) { 
  // return token_get_pbuffer(&parser->token);
  // problem with linker. Maybe because its inlined?
  return &parser->matched_token.buffer[0];
}

/* **** Tree specific methods **** */



static void parse_trees(PhyParser parser) {
  TokenType ttype;
  PhyTree* tree;
  // we allow only one tree for the time being
  ttype = next_token(parser); // start
  tree = parse_tree(parser);
  list_any_push_back(parser->tree_list,tree);
  match_token(parser,TOKEN_SCANEOF);
  return;
}


static PhyTree* parse_tree(PhyParser parser) {
  TokenType ttype;
  PhyTree* tree = NULL;
  PhyNode* root = NULL;

  root = parse_internal(parser,NULL,0);
  match_token(parser,TOKEN_SEMICOLON);
  tree = phytree_new(root,parser->internal_list,parser->tip_list);
  return tree;
}

static void parse_nodelist(PhyParser parser,PhyNode* anc, ListAny nodes) {
  TokenType ttype;
  PhyNode* node;
  // token in memory 
  node = parse_node(parser,anc);
  list_any_push_back(nodes, node);
  ttype = get_tokenType(parser);
  while (ttype == TOKEN_COMMA) {
    match_token(parser,TOKEN_COMMA);
    // read next token
    node = parse_node(parser,anc);
    list_any_push_back(nodes,node);
    ttype = get_tokenType(parser);
  }
  return;
}

static PhyNode* parse_node(PhyParser parser, PhyNode* anc) {
  TokenType ttype;
  PhyNode* node;
  // token in memory
  ttype = get_tokenType(parser);
  if (ttype == TOKEN_ID) {
    node = parse_tip(parser,anc);
  } else if (ttype == TOKEN_LPAREN) {
    node = parse_internal(parser,anc,1);
  } else {
    printf("Tree parse error. Expecting new tree.\nFound: ");
    token_print(&parser->token,stdout);
    // this->token.print(std::cerr);
    //throw ParserExc(this->token,"Expecting new tree");
    longjmp(parser->parser_failed,1);
  }
  return node;
}

static PhyNode* parse_tip(PhyParser parser, PhyNode* anc) {
  PhyNode* node;
  float l = 0.0;

  parser->tip_counter++;
  parser->node_counter++;
  node = phynode_new(0,anc,1);
  list_any_push_back(parser->tip_list, node);
  match_token(parser,TOKEN_ID);
  //printf("TIP matched %s\n",get_matched_buffer(parser));
  phynode_set_taxa(node, get_matched_buffer(parser));
  match_token(parser, TOKEN_COLON);
  match_token(parser, TOKEN_FLOAT);

  l = atof(get_matched_buffer(parser));
  node->edge_length = l;

  return node;
}

PhyNode* parse_internal(PhyParser parser, PhyNode* anc, int distance_must) {
  float l = 0.0;
  PhyNode* node;
  TokenType ttype;
  ListAny children_list = list_any_new();
  ListAny lc;
  void *lvoid;

  parser->internal_counter++;
  parser->node_counter++;
  node = phynode_new(0,anc,0); 
  list_any_push_back(parser->internal_list,node);
  match_token(parser, TOKEN_LPAREN);
  list_any_push_front(parser->children_stack,children_list);
  parse_nodelist(parser,node,children_list);
  list_any_pop_front(parser->children_stack,&lvoid);
  lc = (ListAny)lvoid;
  match_token(parser,TOKEN_RPAREN);
  if (distance_must) {
    match_token(parser, TOKEN_COLON);
    match_token(parser, TOKEN_FLOAT);
    l = atof(get_matched_buffer(parser));
  } else {
    ttype = get_tokenType(parser);
    if (ttype == TOKEN_COLON) {  // optional
      match_token(parser, TOKEN_COLON);
      match_token(parser, TOKEN_FLOAT);
    l = atof(get_matched_buffer(parser));
    }
  }
  phynode_set_descendants(node,children_list);
  node->edge_length = l;
  list_any_free(&children_list);
  return node;
}

static void clear_state(PhyParser parser) {
  PhyTree *tree;
  PhyNode *node;
  ListAny l, children;
  void *pvoid;

  // Print statuc - debug
  //printf("Clearing state\n");
  //printf("stack %d ",list_any_size(parser->children_stack));
  //printf("internal %d ",list_any_size(parser->internal_list));
  //printf("tip %d ",list_any_size(parser->tip_list));
  //printf("tlist %d\n",list_any_size(parser->tree_list));

  // Unravel children_stack and clean list structures
  // the nodes are already part of 
  l = parser->children_stack;
  while (!list_any_empty(l)) {
    list_any_pop_front(l,&pvoid);
    children = (ListAny)pvoid;
    list_any_free(&children); // frees list_nodes and list
  } 
  // children_stack is empty
  // Nodes must be freed properly
  l = parser->internal_list;
  while (!list_any_empty(l)) {
    list_any_pop_front(l,&pvoid);
    node = (PhyNode*)pvoid;
    phynode_free(&node);
  }
  l = parser->tip_list;
  while (!list_any_empty(l)) {
    list_any_pop_front(l,&pvoid);
    node = (PhyNode*)pvoid;
    phynode_free(&node);
  }

  // same with tree_list
  l = parser->tree_list;
  while (!list_any_empty(l)) {
    list_any_pop_front(l,&pvoid);
    tree = (PhyTree*)pvoid;
    phytree_free(&tree);
  }

  return;
}
