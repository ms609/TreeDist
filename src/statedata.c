/*
//  statedata.c
//  MorPhy2
//
//  Created by mbrazeau on 26/04/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
*/
#include "mpl.h"
#include "morphydefs.h"
#include "morphy.h"
#include "mplerror.h"
#include "statedata.h"

char* mpl_skip_closure(const char *closure, const char openc, const char closec)
{
    if (*closure != openc) {
        return (char*)ERR_BAD_PARAM;
    }
    char *ret = (char*)closure;
    
    do {
        ++ret;
        if (*ret == closec) {
            return ret;
        }
    } while (*ret);
    
    return NULL;
}


int compare_char_t_states(const void *ptr1, const void *ptr2)
{
    return *(char*)ptr1 - *(char*)ptr2;
}


int mpl_compare_symbol_lists(const char* sym1, const char* sym2)
{
    int i = 0;
    
    for (i = 0; sym1[i]; ++i) {
        if (!strchr(sym2, sym1[i])) {
            if (!isspace(sym1[i])) {
                return 1;
            }
        }
    }
    
    for (i = 0; sym2[i]; ++i) {
        if (!strchr(sym1, sym2[i])) {
            if (!isspace(sym2[i])) {
                return 1;
            }
        }
    }
    
    return 0;
}


int mpl_assign_symbol_list_from_matrix
(const char *symbs, MPLsymbols* symlist)
{
    assert(symbs && symlist);
    
    int nsymbs = (int)strlen(symbs);
    
    ++nsymbs;
    symlist->symbolsinmatrix = (char*)calloc(nsymbs, sizeof(char));
    
    if (!symlist->symbolsinmatrix) {
        return ERR_BAD_MALLOC;
    }
    
    strcpy(symlist->symbolsinmatrix, symbs);
    
    return ERR_NO_ERROR;
}


char *mpl_query_symbols_from_matrix(Morphyp m)
{
    return m->symbols.symbolsinmatrix;
}


int mpl_get_states_from_rawdata(Morphyp handl)
{
    
    assert(handl);
    
    int count = 0;
    char *rawmatrix = handl->char_t_matrix;
    char *current = NULL;
    int listmax = MAXSTATES + 1; // +1 for terminal null.
    char* statesymbols = (char*)calloc(listmax, sizeof(char));//[listmax];
//    int dbg_loopcount = 0;
    
    statesymbols[0] = '\0';
    current = rawmatrix;
    
    do {
        if (strchr(VALIDSYMB, *current)) {
            
            if (strchr(VALID_NEXMAT_PUNC, *current)) {
                ++current;
            }
            if (!strchr(statesymbols, *current) &&
                strchr(VALID_STATESYMB, *current)) {
                // Put in list
                statesymbols[count] = *current;
                ++count;
                statesymbols[count] = '\0';
            }
        }
        else {
            return ERR_INVALID_SYMBOL;
        }
    
        ++current;
    
    } while (*current);
    
    // Sort alphanumerically
    qsort(statesymbols, strlen(statesymbols), sizeof(char),
          compare_char_t_states);
    
    int numstates = (int)strlen(statesymbols);
    mpl_set_numsymbols(numstates, handl);
    mpl_assign_symbol_list_from_matrix(statesymbols, &handl->symbols);
    free(statesymbols);
    return count-1;
}


int mpl_set_numsymbols(int numsymb, Morphyp handl)
{
    assert(handl);
    handl->symbols.numstates = numsymb;
    return ERR_NO_ERROR;
}


int mpl_get_numsymbols(Morphyp handl)
{
    assert(handl);
    return handl->symbols.numstates;
}


int mpl_create_state_dictionary(Morphyp handl)
{
    int i           = 0;
    int gappush     = 0;
    int numsymbs    = handl->symbols.numstates;
    mpl_get_symbols((Morphy)handl);
    
    if (!handl->symbols.packed) {
        
        handl->symbols.packed = (MPLstate*)calloc(handl->symbols.numstates,
                                                  sizeof(MPLstate));
        if (!handl->symbols.packed) {
            return ERR_BAD_MALLOC;
        }
    }

    if (handl->gaphandl == GAP_INAPPLIC || handl->gaphandl == GAP_NEWSTATE) {
        gappush = 1;
    }
    
    for (i = 0; i < numsymbs; ++i) {
        handl->symbols.packed[i] = 1 << (i + gappush);
    }
    
    return ERR_NO_ERROR;
}


MPLstate mpl_convert_gap_symbol(Morphyp handl, bool over_cutoff)
{
    if (handl->gaphandl == GAP_INAPPLIC) {
        if (over_cutoff) {
            return NA;
        }
        else {
            return MISSING;
        }
    }
    else if (handl->gaphandl == GAP_NEWSTATE) {
        return (MPLstate)1;
    }
    else if (handl->gaphandl == GAP_MISSING) {
        return MISSING;
    }
    
    return ERR_NO_DATA;
}


MPLstate mpl_convert_char_to_MPLstate(const char* celldata, Morphyp handl)
{
    int i = 0;
    MPLstate result = 0;
    
    do {
        i = 0;
        do {
            if (*celldata == handl->symbols.statesymbols[i]) {
                result |= handl->symbols.packed[i];
            }
            ++i;
        } while (handl->symbols.statesymbols[i]);
        ++celldata;
    } while (*celldata);
    
    return result;
}

// TODO: Call this function during the 'apply data' routine.
int mpl_convert_cells(Morphyp handl)
{
    
    int i = 0;
    int j = 0;
    int ncols = mpl_get_num_charac((Morphy)handl);
    int nrows = mpl_get_numtaxa((Morphy)handl);
    MPLmatrix *inmatrix = &handl->inmatrix;
    MPLcharinfo* chinfo = handl->charinfo;
    MPLcell *cell;
   
    if (handl->gaphandl == GAP_INAPPLIC) {
        mpl_count_gaps_in_columns(handl);
    }
    
    char *celldata = NULL;
    
    for (i = 0; i < ncols; ++i) {
        
        for (j = 0; j < nrows; ++j) {
            
            cell = &inmatrix->cells[j * ncols + i];
            celldata = cell->asstr;
            
            if (*celldata == handl->symbols.gap) {
                
                bool over_cutoff = false;
                
                if (chinfo[i].ninapplics > NACUTOFF) {
                    over_cutoff = true;
                }
                
                cell->asint = mpl_convert_gap_symbol(handl, over_cutoff);
            }
            else if (*celldata == handl->symbols.missing) {
                cell->asint = MISSING;
            }
            else {
                cell->asint = mpl_convert_char_to_MPLstate(celldata, handl);
            }
            
        }
    
    }
    
    return ERR_NO_ERROR;
}


void mpl_destroy_symbolset(Morphyp m)
{
    assert(m);
    if (m->symbols.statesymbols) {
        if (m->symbols.statesymbols == m->symbols.symbolsinmatrix) {
            free(m->symbols.statesymbols);
            m->symbols.statesymbols = NULL;
            m->symbols.symbolsinmatrix = NULL;
        }
        else {
            free(m->symbols.statesymbols);
            m->symbols.statesymbols = NULL;
            if (m->symbols.symbolsinmatrix) {
                free(m->symbols.symbolsinmatrix);
                m->symbols.symbolsinmatrix = NULL;
            }
        }
    }
    if (m->symbols.packed) {
        free(m->symbols.packed);
        m->symbols.packed = NULL;
    }
}


bool mpl_is_valid_matrix_symbol(const char c)
{
    if (strchr(VALID_STATESYMB, c)) {
        return true;
    }
    else if (strchr(VALID_WILDCAR, c)) {
        return true;
    }
    else if (strchr(VALID_NEXMAT_PUNC, c)) {
        return true;
    }
    
    return false;
}


unsigned long mpl_get_valid_matrix_length(const char* rawmatrix)
{
    unsigned long len = 0;
    char* matptr = (char*)rawmatrix;
    
    do {
        if (mpl_is_valid_matrix_symbol(*matptr)) {
            ++len;
        }
        else if (*matptr == '[') {
            matptr = mpl_skip_closure(matptr, '[', ']');
            assert(!(matptr < 0));
        }
        ++matptr;
    } while (*matptr);
    
    return len;
}


void mpl_copy_valid_matrix_data(char *copy, const char* rawmatrix)
{
    int i = 0;
    char* matptr = (char*)rawmatrix;
    
    do {
        if (mpl_is_valid_matrix_symbol(*matptr)) {
            copy[i] = *matptr;
            ++i;
        }
        else if (*matptr == '[') {
            matptr = mpl_skip_closure(matptr, '[', ']');
            assert(!(matptr < 0));
        }
        ++matptr;
    } while (*matptr);
    
    copy[i-1] = '\0';
}


// Copy the raw matrix, take out whitespace and comments
int mpl_copy_raw_matrix(const char* rawmatrix, Morphyp handl)
{
    unsigned long len = mpl_get_valid_matrix_length(rawmatrix);
    
    char *matcpy = (char*)calloc(len + 1, sizeof(char));
   
    if (!matcpy) {
        return ERR_BAD_MALLOC;
    }
    mpl_copy_valid_matrix_data(matcpy, rawmatrix);
    handl->char_t_matrix = matcpy;
    return ERR_NO_ERROR;
}


int mpl_check_nexus_matrix_dimensions
(char *preproc_matrix, int input_num_taxa, int input_num_chars)
{
    /* An input matrix should not have inline taxon names. This function
     * scans each row of the input matrix to determine whether or not the
     * number of places in the row corresponds to input number of
     * of characters. If the number exceeds the expected number of data
     * columns (num_input_chars), then it is inferred that taxon names or
     * other extraneous info are included in the matrix. */
    
    char* current = NULL;
    int matrix_size = 0;
    int expected_size = 0;
    
    expected_size = input_num_chars * input_num_taxa;
    
    current = preproc_matrix;
    assert(current);
    
    do {
        if (strchr(VALID_STATESYMB, *current)
            || strchr(VALID_WILDCAR, *current)) {
            ++matrix_size;
        }
        else if (*current == '(' || *current == '{') {
            
            char* err = 0;
            
            if (*current == '(') {
                err = mpl_skip_closure(current, '(', ')');
            }
            else {
                err = mpl_skip_closure(current, '{', '}');
            }
            if (*err <= 0) {
                return ERR_MATCHING_PARENTHS;
            }
            
            current = err;
            assert(!(current < 0));
            ++matrix_size;
        }

        ++current;
    } while (*current);
    
    if (matrix_size > expected_size) {
        return ERR_DIMENS_UNDER;
    }
    else if (matrix_size < expected_size) {
        return ERR_DIMENS_OVER;
    }

    return ERR_NO_ERROR;
}


char* mpl_get_preprocessed_matrix(Morphyp handl)
{
    assert(handl);
    return handl->char_t_matrix;
}


MPLstate mpl_gap_value(Morphyp handl)
{
    switch (mpl_get_gaphandl(handl)) {
        case GAP_INAPPLIC:
            return NA;
        case GAP_MISSING:
            return MISSING;
        case GAP_NEWSTATE:
            return (MPLstate)1;
        case GAP_MAX:
            return -1;
        default:
            break;
    }
    
    return -2;
}

int mpl_init_inmatrix(Morphyp handl)
{
    assert(handl);
    MPLmatrix* mat = &handl->inmatrix;
    int ntaxa = mpl_get_numtaxa((Morphyp)handl);
    int nchar = mpl_get_num_charac((Morphyp)handl);
    int nstates = mpl_get_numsymbols(handl);
    
//    mat->chtypes = (MPLchtype*)calloc(nchar, sizeof(MPLchtype));
//    if (!mat->chtypes) {
//        return ERR_BAD_MALLOC;
//    }
//    
//    mat->intweights = (int*)calloc(nchar, sizeof(int));
//    if (!mat->intweights) {
//        mpl_delete_mpl_matrix(mat);
//        return ERR_BAD_MALLOC;
//    }
//    
//    mat->fltweights = (Mflt*)calloc(nchar, sizeof(Mflt));
//    if (!mat->fltweights) {
//        mpl_delete_mpl_matrix(mat);
//        return ERR_BAD_MALLOC;
//    }
    
    mat->cells = (MPLcell*)calloc(ntaxa * nchar, sizeof(MPLcell));
    if (!mat->cells) {
        mpl_delete_mpl_matrix(mat);
        return ERR_BAD_MALLOC;
    }
    
    mat->ncells = ntaxa * nchar;
    int i = 0;
    
    for (i = 0; i < mat->ncells; ++i) {
        mat->cells[i].asstr = (char*)calloc(nstates + 1, sizeof(char));
        if (!mat->cells[i].asstr) {
            int j = 0;
            for (j = 0; j < i; ++j) {
                free(mat->cells[i].asstr);
                mat->cells[i].asstr = NULL;
            }
            mpl_delete_mpl_matrix(mat);
            return ERR_BAD_MALLOC;
        }
    }
    
    return ERR_NO_ERROR;
}


int mpl_delete_mpl_matrix(MPLmatrix* m)
{
    if (!m) {
        return ERR_BAD_PARAM;
    }
    
    int i = 0;
    
    if (m->cells) {
        for (i = 0; i < m->ncells; ++i) {
            if (m->cells[i].asstr) {
                free(m->cells[i].asstr);
                m->cells[i].asstr = NULL;
            }
        }
        free(m->cells);
        m->cells = NULL;
    }
    
//    if (m->chtypes) {
//        free(m->chtypes);
//        m->chtypes = NULL;
//    }
//    
//    if (m->fltweights) {
//        free(m->fltweights);
//        m->fltweights = NULL;
//    }
//    
//    if (m->intweights) {
//        free(m->intweights);
//        m->intweights = NULL;
//    }
    
    return ERR_NO_ERROR;
}


MPLmatrix* mpl_get_mpl_matrix(Morphyp m)
{
    return &m->inmatrix;
}


int mpl_set_gap_push(Morphyp handl)
{
    MPLgap_t gt = mpl_get_gaphandl(handl);
    
    if (gt == GAP_INAPPLIC || gt == GAP_NEWSTATE) {
        return 1;
    }
    else if (gt == GAP_MISSING) {
        return 0;
    }
    
    return -1;
}


int mpl_get_uncorrected_shift_value(char symb, Morphyp handl)
{
    // Gets the raw shift value as determined by the order in the symbols list
    assert(symb != DEFAULTGAP && symb != DEFAULTMISSING);
    int shift = 0;
    char* symbols = mpl_get_symbols((Morphy)handl);
    
    while (*symbols != symb && *symbols) {
        ++symbols;
        ++shift;
    }
    
    return shift;
}


void mpl_use_symbols_from_matrix(Morphyp handl)
{
    handl->symbols.statesymbols = handl->symbols.symbolsinmatrix;
}


int mpl_write_input_rawchars_to_cells(Morphyp handl)
{
    assert(handl);
    int i = 0;
    int j = 0;
//    int rows = mpl_get_numtaxa((Morphyp)handl);
//    int cols = mpl_get_num_charac((Morphyp)handl);
//    int length = rows * cols;
    
    char* prpdata = mpl_get_preprocessed_matrix(handl);
    
    while (*prpdata) {
        
        if (!strchr(VALID_NEXMAT_PUNC, *prpdata)) {
            handl->inmatrix.cells[i].asstr[0] = *prpdata;
            handl->inmatrix.cells[i].asstr[1] = '\0';
        }
        else {
            if (*prpdata == '(') {
                j = 0;
                ++prpdata;
                do {
                    
                    handl->inmatrix.cells[i].asstr[j] = *prpdata;
                    ++j;
                    ++prpdata;
                } while (*prpdata != ')');
                handl->inmatrix.cells[i].asstr[j] = '\0';
            }
            
            if (*prpdata == '{') {
                j = 0;
                ++prpdata;
                do {
                    handl->inmatrix.cells[i].asstr[j] = *prpdata;
                    ++j;
                    ++prpdata;
                } while (*prpdata != '}');
                handl->inmatrix.cells[i].asstr[j] = '\0';
            }
            if (*prpdata == ';') {
                break;
            }
        }

        ++i;
        ++prpdata;
    };
    
    //prpdata = mpl_get_preprocessed_matrix(handl);
    
    return ERR_NO_ERROR;
}

// TODO: Rename this.
int mpl_preproc_rawdata(Morphyp handl)
{
    int ret = ERR_NO_ERROR;
    
    mpl_get_states_from_rawdata(handl);
    
    // If no symbols supplied by the user
    if (!mpl_get_symbols((Morphy)handl)) {
        // Assign internal symbols to list
        mpl_use_symbols_from_matrix(handl);
    }
    else {
        
        char *frommatrix    = mpl_query_symbols_from_matrix(handl);
        char *user          = mpl_get_symbols((Morphyp)handl);
        
        if (mpl_compare_symbol_lists(frommatrix, user)) {
            return ERR_SYMBOL_MISMATCH;
        }
    }
    
    mpl_init_inmatrix(handl);
    
    // Now safe to write characters into cells.
    mpl_write_input_rawchars_to_cells(handl);
    
    return ret;
}

// TODO: This probably needs to be more memory-safe
char *mpl_translate_state2char(MPLstate cstates, Morphyp handl)
{
    int i = 0;
    int shift = 0;
    int gapshift = 0;
   
    MPLgap_t gaphandl = mpl_query_gaphandl((Morphyp)handl);
    if (gaphandl == GAP_INAPPLIC || gaphandl == GAP_NEWSTATE) {
        gapshift = 1;
    }
    char *res = calloc(MAXSTATES+1, sizeof(char));
    if (!res) {
        return NULL;
    }
    char* symbols = mpl_get_symbols((Morphy)handl);
    
    if (cstates < (MISSING-NA)) {
        while (cstates) {
            if (1 & cstates) {
                if (shift == 0 && gapshift) {
                    res[i] = mpl_get_gap_symbol(handl);
                }
                else {
                    res[i] = symbols[shift - gapshift];
                }
                ++i;
            }
            cstates = cstates >> 1;
            ++shift;
        } 
    }
    else {
        res[0] = '?';
    }
    
    return res;
}

int mpl_init_charac_info(Morphyp handl)
{
    int nchar = mpl_get_num_charac((Morphy)handl);
    
    if (handl->charinfo) {
        free(handl->charinfo);
    }
    
    handl->charinfo = (MPLcharinfo*)calloc(nchar, sizeof(MPLcharinfo));
    if (!handl->charinfo) {
        return ERR_BAD_MALLOC;
    }
    
    int i = 0;
    for (i = 0; i < nchar; ++i) {
        handl->charinfo[i].charindex    = i;
        handl->charinfo[i].chtype       = DEFAULCHARTYPE;
        handl->charinfo[i].realweight   = 1.0;
        handl->charinfo[i].basewt       = 1;
        handl->charinfo[i].intwt        = 1;
    }
    
    return ERR_NO_ERROR;
}

void mpl_delete_charac_info(Morphyp handl)
{
    assert(handl);
    if (!handl->charinfo) {
        return;
    }
    free(handl->charinfo);
}
