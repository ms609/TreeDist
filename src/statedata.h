//
//  statedata.h
//  MorPhy2
//
//  Created by mbrazeau on 26/04/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//

#ifndef statedata_h
#define statedata_h

#include <limits.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>


#define VALID_NEXMAT_PUNC   "{}();"
#define VALID_XREAD_MATPUNC "[];"
#define VALID_WILDCAR       "-?"
#define VALID_STATESYMB     "+0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
#define VALID_WS            "\n\t "
#define VALIDSYMB   VALID_NEXMAT_PUNC VALID_XREAD_MATPUNC VALID_WILDCAR \
                    VALID_STATESYMB VALID_WS


/* Function prototypes */
int         mpl_init_symbolset(Morphyp m);
int         mpl_set_numsymbols(int numsymb, Morphyp handl);
int         mpl_get_numsymbols(Morphyp handl);
void        mpl_destroy_symbolset(Morphyp m);
char*       mpl_skip_closure(const char *closure, const char openc, const char closec);
int         mpl_compare_symbol_lists(const char* sym1, const char* sym2);
int         mpl_assign_symbol_list_from_matrix(const char *symbs, MPLsymbols* symlist);
char*       mpl_query_symbols_from_matrix(Morphyp m);
int         mpl_get_states_from_rawdata(Morphyp handl);
int         mpl_copy_raw_matrix(const char* rawmatrix, Morphyp handl);
int         mpl_check_nexus_matrix_dimensions(char *input_matrix, int input_num_taxa, int input_num_chars);
char*       mpl_get_preprocessed_matrix(Morphyp handl);
int         mpl_write_input_rawchars_to_cells(Morphyp handl);
int         mpl_create_state_dictionary(Morphyp handl);
int         mpl_convert_cells(Morphyp handl);
int         mpl_preproc_rawdata(Morphyp handl);
MPLmatrix*  mpl_new_mpl_matrix(const int ntaxa, const int nchar, const int nstates);
int         mpl_delete_mpl_matrix(MPLmatrix* m);
MPLmatrix*  mpl_get_mpl_matrix(Morphyp m);
char*       mpl_translate_state2char(MPLstate cstates, Morphyp handl);
int         mpl_init_charac_info(Morphyp handl);
void        mpl_delete_charac_info(Morphyp handl);

#endif /* statedata_h */
