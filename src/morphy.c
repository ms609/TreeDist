/*
//  morphy.c
//  MorPhy2
//
//  Created by mbrazeau on 23/04/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
*/
#include "mpl.h"
#include "morphydefs.h"
#include "morphy.h"
#include "mplerror.h"
#include "statedata.h"
#include "fitch.h"
#include "wagner.h"

void *mpl_alloc(size_t size, int setval)
{
    void *ret = malloc(size);
    if (ret) {
        memset(ret, setval, size);
    }
    return ret;
}


Morphyp mpl_new_Morphy_t(void)
{
    Morphyp new = (Morphyp)calloc(1, sizeof(Morphy_t));
    
    mpl_set_gaphandl(GAP_INAPPLIC, (Morphy)new);
    new->symbols.gap        = DEFAULTGAP;
    new->symbols.missing    = DEFAULTMISSING;
    new->nthreads           = 1; 
    new->usrwtbase          = 0;
    new->wtbase             = 1;
    
    return new;
}


/*void* mpl_get_from_matrix
(const int      row,
 const int      col,
 const int      ncol,
 const size_t   size,
 const void*    data)
{
    return (void*)(data + (row * ncol * size + (size * col)));
}*/


int mpl_get_gaphandl(Morphyp handl)
{
    assert(handl);
    return handl->gaphandl;
}


int mpl_set_num_charac(const int nchar, Morphyp m)
{
    if (!m) {
        return ERR_BAD_PARAM;
    }
    
    m->numcharacters = nchar;
    
    return ERR_NO_ERROR;
}


int mpl_set_numtaxa(const int ntax, Morphyp m)
{
    if (!m) {
        return ERR_BAD_PARAM;
    }
    
    m->numtaxa = ntax;
    
    return ERR_NO_ERROR;
}


int mpl_check_data_loaded(Morphyp m)
{
    if (m->char_t_matrix) {
        return 1;
    }
    
    return 0;
}


char mpl_get_gap_symbol(Morphyp handl)
{
    return handl->symbols.gap;
}

void mpl_flt_rational_approx
(unsigned long *a, unsigned long *b, const double fval)
{
   /* Using David Eppstein's method 
    http://www.ics.uci.edu/~eppstein/numth/frap.c*/
    
    
    long m[2][2];
    long ai;
    long maxden = 100;
    double x = fval;
    
    m[0][0] = m[1][1] = 1;
    m[0][1] = m[1][0] = 0;
    
    while (m[1][0] * (ai = (long)x) + m[1][1] <= maxden) {
        long t;
        t = m[0][0] * ai + m[0][1];
        m[0][1] = m[0][0];
        m[0][0] = t;
        t = m[1][0] * ai + m[1][1];
        m[1][1] = m[1][0];
        m[1][0] = t;
        if (x == (double)ai) {
            break;
        }
        x = 1 /(x - (double) ai);
        if (x > (double)INT_MAX) {
            break;
        }
    }
    
    *a = m[0][0];
    *b = m[1][0];
}


bool mpl_almost_equal(double a, double b)
{
    double diff = fabs(a - b);
    double largest = 0.0;
    
    a = fabs(a);
    b = fabs(b);
   
    largest = (b > a) ? b : a;
    
    if (diff <= largest * MPL_EPSILON) {
        return true;
    }
    
    return false;
}

/*!
 @brief Checks whether or not a value corresponds to a real number or is whole
 @discussion Should be used only to check values from external calling functions
 that are not expected to be the product of a calculation (i.e. are a user-
 supplied value, such as an input weight).
 @param n The value to be tested.
 @return A true/false value: true if value is fractional, false if value is
 whole.
 */
bool mpl_isreal(const double n)
{
    assert(!(n > (double)LONG_MAX));
    long i = (long)n;
    if (n == (double)i) {
        return false;
    }
    return true;
}


int mpl_change_weight_base(const unsigned long wtbase, Morphyp handl)
{
    if (handl->usrwtbase) {
        return 1;
    }
    
    handl->wtbase = wtbase;
    
    return 0;
}

static inline unsigned long mpl_greatest_common_denom
(unsigned long a, unsigned long b)
{
    unsigned long t = 0;
    while (b) {
        t = b;
        b = a % b;
        a = t;
    }
    
    return a;
}

unsigned long mpl_least_common_multiple(unsigned long a, unsigned long b)
{
    return (a * b) / mpl_greatest_common_denom(a, b);
}

void mpl_set_new_weight_public
(const double wt, const int char_id, Morphyp handl)
{
    bool wtisreal = mpl_isreal(wt);
  
    
    if (wtisreal) {
        
        if (!mpl_isreal(handl->charinfo[char_id].realweight) ||
            handl->charinfo[char_id].realweight == 0.0)
        {
            ++handl->numrealwts;
        }
    }
    else {
        
        if (mpl_isreal(handl->charinfo[char_id].realweight)) {
            --handl->numrealwts;
        }
    }
    
    handl->charinfo[char_id].realweight = wt;
}

void mpl_scale_all_intweights(Morphyp handl)
{
    int i = 0;
    int nchar = mpl_get_num_charac((Morphy)handl);
    
   
    if (!handl->numrealwts) {
        if (handl->usrwtbase != 0) {
            handl->wtbase = handl->usrwtbase;
        } else {
            handl->wtbase = DEFAULTWTBASE;
        }
    }
    
    for (i = 0; i < nchar; ++i) {
        mpl_flt_rational_approx(&handl->charinfo[i].intwt,
                                &handl->charinfo[i].basewt,
                                handl->charinfo[i].realweight);
    }
    
    for (i = 0; i < nchar; ++i) {
        handl->wtbase = mpl_least_common_multiple(handl->wtbase,
                                                  handl->charinfo[i].basewt);
    }
    
    for (i = 0; i < nchar; ++i) {
        handl->charinfo[i].intwt *= (handl->wtbase / handl->charinfo[i].basewt);
        handl->charinfo[i].basewt = handl->wtbase;
    }
        
}

void mpl_assign_fitch_fxns(MPLpartition* part)
{
    assert(part);
    
    if (part->isNAtype) {
        part->inappdownfxn      = mpl_NA_fitch_second_downpass;
        part->inappupfxn        = mpl_NA_fitch_second_uppass;
        part->prelimfxn         = mpl_NA_fitch_first_downpass;
        part->finalfxn          = mpl_NA_fitch_first_uppass;
        part->tipupdate         = mpl_fitch_NA_tip_update;
        part->tipfinalize       = mpl_fitch_NA_tip_finalize;
        part->tiproot           = mpl_fitch_NA_first_one_branch;
        part->tiprootfinal      = mpl_fitch_NA_second_one_branch;
        part->loclfxn           = mpl_fitch_NA_local_reopt;
        part->downrecalc1       = mpl_NA_fitch_first_update_downpass;
        part->uprecalc1         = mpl_NA_fitch_first_update_uppass;
        part->inappdownrecalc2  = mpl_NA_fitch_second_update_downpass;
        part->inapuprecalc2     = mpl_NA_fitch_second_update_uppass;
        part->tipupdaterecalc   = mpl_fitch_NA_tip_recalc_update;
        part->tiprootrecalc     = mpl_fitch_NA_first_one_branch;
        part->tiprootupdaterecalc = mpl_fitch_NA_second_one_branch_recalc;
    }
    else {
        part->prelimfxn         = mpl_fitch_downpass;
        part->finalfxn          = mpl_fitch_uppass;
        part->tipupdate         = mpl_fitch_tip_update;
        part->tiproot           = mpl_fitch_one_branch;
        part->tipfinalize       = NULL;
        part->inappdownfxn      = NULL;
        part->inappupfxn        = NULL;
        part->loclfxn           = mpl_fitch_local_reopt;
        part->downrecalc1       = NULL;
        part->uprecalc1         = NULL;
        part->inappdownrecalc2  = NULL;
        part->inapuprecalc2     = NULL;
    }
}

void mpl_assign_wagner_fxns(MPLpartition* part)
{
    assert(part);
    
        part->prelimfxn     = mpl_wagner_downpass;
        part->finalfxn      = mpl_wagner_uppass;
        part->tipupdate     = mpl_wagner_tip_update;
        part->tipfinalize   = NULL;
        part->inappdownfxn  = NULL; 
        part->inappupfxn    = NULL;
        part->loclfxn       = NULL;

}



int mpl_fetch_parsim_fxn_setter
(void(**pars_assign)(MPLpartition*), MPLchtype chtype)
{
    int err = ERR_NO_ERROR;
    
    switch (chtype) {
        case FITCH_T:
            if (pars_assign) {
                *pars_assign = mpl_assign_fitch_fxns;
            }
            break;
        case WAGNER_T:
            if (pars_assign) {
                *pars_assign = mpl_assign_wagner_fxns;
            }
            break;
            
            
        default:
            err = ERR_CASE_NOT_IMPL;
            break;
    }
    
    return err;
}


int mpl_assign_partition_fxns(MPLpartition* part)
{
    assert(part);
    int err = ERR_NO_ERROR;

    void (*pars_assign)(MPLpartition*) = NULL;
    
    err = mpl_fetch_parsim_fxn_setter(&pars_assign, part->chtype);
    
    if (!err && pars_assign) {
        pars_assign(part);
    }
    
    return err;
}


int mpl_extend_intarray(int** array, size_t size)
{
    int* temp = (int*)realloc(*array, size);
    if (!temp) {
        return ERR_BAD_MALLOC;
    }
    
    *array = temp;
    
    return ERR_NO_ERROR;
}


int mpl_part_push_index(int newint, MPLpartition* part)
{
    int err = ERR_NO_ERROR;
    
    if (part->ncharsinpart < part->maxnchars) {
        part->charindices[part->ncharsinpart] = newint;
        ++part->ncharsinpart;
    }
    else {
        err = mpl_extend_intarray(&part->charindices,
                                 (part->maxnchars + 1) * sizeof(int));
        if (!err) {
            part->charindices[part->ncharsinpart] = newint;
            ++part->ncharsinpart;
            ++part->maxnchars;
        }
    }
    
    return err;
}


int mpl_part_remove_index(int index, MPLpartition* part)
{
    if (!part->ncharsinpart) {
        return 1;
    }
    
    --part->ncharsinpart;
    assert(part->ncharsinpart >= 0);
    
    int i = 0;
    for (i = 0; i < part->ncharsinpart; ++i) {
        part->charindices[i] = part->charindices[i + 1];
    }
    part->charindices[i] = MPLCHARMAX; 
    
    return 0;
}


int mpl_delete_partition(MPLpartition* part)
{
    int err = ERR_UNEXP_NULLPTR;
    
    if (part) {
        
        if (part->charindices) {
            free(part->charindices);
            part->charindices   = NULL;
        }
        if (part->intwts) {
            free(part->intwts);
            part->intwts = NULL;
        }
        part->maxnchars     = 0;
        part->ncharsinpart  = 0;
        part->chtype        = NONE_T;
        part->tipupdate     = NULL;
        part->tipfinalize   = NULL;
        part->inappdownfxn  = NULL;
        part->inappupfxn    = NULL;
        part->prelimfxn     = NULL;
        part->finalfxn      = NULL;
        part->next          = NULL;
        free(part);
        err = ERR_NO_ERROR;
    }
    
    return err;
}

int mpl_delete_all_partitions(Morphyp handl)
{
    assert(handl);
    int i = 0;
    
    if (handl->numparts) {
        
        MPLpartition* p = handl->partstack;
        MPLpartition* q = NULL;
        while (p) {
            q = p->next;
            mpl_delete_partition(p);
            p = q;
        }
        
        for (i = 0; i < handl->numparts; ++i) {
            handl->partitions[i] = NULL;
        }
        free(handl->partitions);
        handl->partitions = NULL;
        
        return ERR_NO_ERROR;
    }
    return ERR_UNEXP_NULLPTR;
}


MPLpartition* mpl_new_partition(const MPLchtype chtype, const bool hasNA)
{
    assert(chtype);
    
    MPLpartition *new = (MPLpartition*)calloc(1, sizeof(MPLpartition));
    
    if (!new) {
        return NULL;
    }
    
    new->chtype     = chtype;
    new->isNAtype   = hasNA;

    new->charindices = (int*)calloc(1, sizeof(int));
    if (!new->charindices) {
        free(new);
        return NULL;
    }

    new->maxnchars      = 1;
    new->ncharsinpart   = 0;
    
    mpl_assign_partition_fxns(new);
    
    return new;
}


int mpl_count_gaps_in_columns(Morphyp handl)
{
    int i = 0;
    int j = 0;
    char gap            = mpl_get_gap_symbol(handl);
    int numchar         = mpl_get_num_charac((Morphy)handl);
    int numtax          = mpl_get_numtaxa((Morphy)handl);
    MPLmatrix* matrix   = mpl_get_mpl_matrix(handl);
    MPLcharinfo* chinfo = handl->charinfo;
    int numna = 0;
    
    for (i = 0; i < numchar; ++i) {
        chinfo[i].ninapplics = 0;
        for (j = 0; j < numtax; ++j) {
            
            MPLcell* cell = &matrix->cells[j * numchar + i];
            
            if (strchr(cell->asstr, gap)) {
                ++chinfo[i].ninapplics;
            }
            
            if (chinfo[i].ninapplics > NACUTOFF) {
                ++numna;
                break;
            }
        }
    }
    
    return numna;
}


int mpl_compare_partition_with_char_info
(const MPLcharinfo *chinfo, const MPLpartition* part, const MPLgap_t gaphandl)
{
    int ret = 0;
    
    if (chinfo->chtype != part->chtype) {
        ++ret;
    }
    
    if (gaphandl == GAP_INAPPLIC) {
        if (chinfo->ninapplics <= NACUTOFF) {
            if (part->isNAtype) {
                ++ret;
            }
        }
        else {
            if (!part->isNAtype) {
                ++ret;
            }
        }
    }
    
    return ret;
}


/*!
 @brief Searches the partition list for a partition matching the supplied info
 @discussion Traverses a linked list of partitions, looking for a partition 
 matching the supplied information. If this function returns NULL, then the
 supplied info does not match a character in the list. A new partition will
 need to be created.
 @param chinfo MPLchtype providing data on a character in the matrix.
 @param part A data partition; should be the first partition in the list.
 @return A pointer to the partition corresponding to the supplied character
 information.
 */
MPLpartition* mpl_search_partitions
(MPLcharinfo *chinfo, MPLpartition* part, MPLgap_t gaphandl)
{
    assert(chinfo);
    MPLpartition* p = part;
    
    while (p) {
        if (!mpl_compare_partition_with_char_info(chinfo, p, gaphandl)) {
            return p;
        }
        p = p->next;
    }
    
    return p;
}


int mpl_compare_partitions(const void* ptr1, const void* ptr2)
{
    MPLpartition* part1 = *(MPLpartition**)ptr1;
    MPLpartition* part2 = *(MPLpartition**)ptr2;
    
    int ret;
    MPLchtype cdiff = NONE_T;
    
    cdiff = part1->chtype - part2->chtype;
    ret = (int)cdiff;
    
    if (!cdiff) {
        if (part2->isNAtype) {
            ret = 1;
        }
        else {
            ret = 0;
        }
    }
    
    return ret;
}


int mpl_put_partitions_in_handle(MPLpartition* first, Morphyp handl)
{
    assert(handl);
    if (!handl->numparts) {
        return ERR_NO_DATA;
    }
    handl->partitions = (MPLpartition**)calloc(handl->numparts,
                                                sizeof(MPLpartition*));
    if (!handl->partitions) {
        return ERR_BAD_MALLOC;
    }
    
    int i = 0;
    MPLpartition* p = first;
    while (p) {
        handl->partitions[i] = p;
        ++i;
        p = p->next;
    }
    assert(i == handl->numparts);
    
    qsort(handl->partitions, handl->numparts, sizeof(MPLpartition*), mpl_compare_partitions);
    handl->partstack = first;
    
    return ERR_NO_ERROR;
}


void mpl_delete_all_update_buffers(Morphyp handl)
{
    /* Erases all the update buffers */
    int i = 0;
    for (i = 0; i < handl->numparts; ++i) {
        MPLpartition* p = handl->partitions[i];
        if (p->update_indices) {
            free(p->update_indices);
            p->update_indices = NULL;
        }
        if (p->update_NA_indices) {
            free(p->update_NA_indices);
            p->update_NA_indices = NULL;
        }
    }
}


int mpl_allocate_update_buffers(Morphyp handl)
{
    /* Allocates memory for an array of indices needed update on the tree */
    int i = 0;
    for (i = 0; i < handl->numparts; ++i) {
        MPLpartition* p = handl->partitions[i];
        p->update_indices = (int*)calloc(p->ncharsinpart, sizeof(int));
        if (!p->update_indices) {
            mpl_delete_all_update_buffers(handl);
            return ERR_BAD_MALLOC;
        }
        p->update_NA_indices = (int*)calloc(p->ncharsinpart, sizeof(int));
        if (!p->update_NA_indices) {
            mpl_delete_all_update_buffers(handl);
            return ERR_BAD_MALLOC;
        }
        
        p->ntoupdate = 0;
        p->nNAtoupdate = 0;
    }
    
    return ERR_NO_ERROR;
}

int mpl_setup_partitions(Morphyp handl)
{
    assert(handl);
    
    int err = ERR_NO_ERROR;
    
    int i = 0;
    int nchar = mpl_get_num_charac((Morphyp)handl);
    
    MPLcharinfo* chinfo = NULL;
    MPLpartition* first = NULL;
    MPLpartition* last  = NULL;
    MPLpartition* p     = NULL;
    int numparts        = 0;
    
    if (handl->partitions) {
        mpl_delete_all_partitions(handl);
    }
    
    for (i = 0; i < nchar; ++i) {
        chinfo = &handl->charinfo[i];
        
        p = mpl_search_partitions(chinfo, first, mpl_get_gaphandl(handl));
        
        if (p) {
            mpl_part_push_index(i, p);
        }
        else {
            bool hasNA = false;
            if (handl->gaphandl == GAP_INAPPLIC) {
                if (chinfo->ninapplics > NACUTOFF) {
                    hasNA = true;
                }
            }
            p = mpl_new_partition(chinfo->chtype, hasNA);
            mpl_part_push_index(i, p);
            if (!first) {
                first = p;
                last = p;
            }
            else {
                last->next = p;
                last = p;
            }
            
            ++numparts;
        }
    }
    
    handl->numparts = numparts;
    err = mpl_put_partitions_in_handle(first, handl);
    mpl_allocate_update_buffers(handl);
    
    return err;
}


int mpl_get_numparts(Morphyp handl)
{
    return handl->numparts;
}


void mpl_delete_nodal_strings(const int nchars, MPLndsets* set)
{
    int i = 0;
    
    for (i = 0; i < nchars; ++i) {
        if (set->downp1str) {
            free(set->downp1str[i]);
            set->downp1str[i] = NULL;
        }
        if (set->upp1str) {
            free(set->upp1str[i]);
            set->upp1str[i] = NULL;
        }
        if (set->downp2str) {
            free(set->downp2str[i]);
            set->downp2str[i] = NULL;
        }
        if (set->upp2str) {
            free(set->upp2str[i]);
            set->upp2str[i] = NULL;
        }
    }
}


int mpl_allocate_stset_stringptrs(const int nchars, MPLndsets* set)
{
    if (!set->downp1str) {
        set->downp1str = (char**)calloc(nchars, sizeof(char*));
        if (!set->downp1str) {
            mpl_delete_nodal_strings(nchars, set);
            return ERR_BAD_MALLOC;
        }
    }
    
    if (!set->upp1str) {
        set->upp1str = (char**)calloc(nchars, sizeof(char*));
        if (!set->upp1str) {
            mpl_delete_nodal_strings(nchars, set);
            return ERR_BAD_MALLOC;
        }
    }
    
    if (!set->downp2str) {
        set->downp2str = (char**)calloc(nchars, sizeof(char*));
        if (!set->downp2str) {
            mpl_delete_nodal_strings(nchars, set);
            return ERR_BAD_MALLOC;
        }
    }
    
    if (!set->upp2str) {
        set->upp2str = (char**)calloc(nchars, sizeof(char*));
        if (!set->upp2str) {
            mpl_delete_nodal_strings(nchars, set);
            return ERR_BAD_MALLOC;
        }
    }
    
    
    return ERR_NO_ERROR;
}


MPLndsets* mpl_alloc_stateset(int numchars)
{
    MPLndsets* new = (MPLndsets*)calloc(1, sizeof(MPLndsets));
    if (!new) {
        return NULL;
    }
    
    new->downpass1 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->downpass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->uppass1 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->uppass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->downpass2 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->downpass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->uppass2 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->uppass2) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->subtree_actives = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->subtree_actives) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->temp_subtr_actives = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->temp_subtr_actives) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->temp_downpass1 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->temp_downpass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->temp_uppass1 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->temp_uppass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->temp_downpass2 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->temp_downpass1) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    
    new->temp_uppass2 = (MPLstate*)calloc(numchars, sizeof(MPLstate));
    if (!new->temp_uppass2) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    new->changes = (bool*)calloc(numchars, sizeof(bool));
    if (!new->changes) {
        mpl_free_stateset(numchars, new);
        return NULL;
    }
    mpl_allocate_stset_stringptrs(numchars, new);
    
    return new;
}


void mpl_free_stateset(const int nchars, MPLndsets* statesets)
{
    if (!statesets) {
        return;
    }
    if (statesets->downpass1) {
        free(statesets->downpass1);
        statesets->downpass1 = NULL;
    }
    if (statesets->uppass1) {
        free(statesets->uppass1);
        statesets->uppass1 = NULL;
    }
    if (statesets->downpass2) {
        free(statesets->downpass2);
        statesets->downpass2 = NULL;
    }
    if (statesets->uppass2) {
        free(statesets->uppass2);
        statesets->uppass2 = NULL;
    }
    if (statesets->subtree_actives) {
        free(statesets->subtree_actives);
        statesets->subtree_actives = NULL;
    }
    if (statesets->temp_subtr_actives) {
        free(statesets->temp_subtr_actives);
        statesets->temp_subtr_actives = NULL;
    }
    if (statesets->temp_downpass1) {
        free(statesets->temp_downpass1);
        statesets->temp_downpass1 = NULL;
    }
    if (statesets->temp_uppass1) {
        free(statesets->temp_uppass1);
        statesets->temp_uppass1 = NULL;
    }
    if (statesets->temp_downpass2) {
        free(statesets->temp_downpass2);
        statesets->temp_downpass2 = NULL;
    }
    if (statesets->temp_uppass2) {
        free(statesets->temp_uppass2);
        statesets->temp_uppass2 = NULL;
    }
    if (statesets->changes) {
        free(statesets->changes);
        statesets->changes = NULL;
    }
    
    mpl_delete_nodal_strings(nchars, statesets);
    
    if (statesets->downp1str) {
        free(statesets->downp1str);
        statesets->downp1str = NULL;
    }
    if (statesets->upp1str) {
        free(statesets->upp1str);
        statesets->upp1str = NULL;
    }
    if (statesets->downp2str) {
        free(statesets->downp2str);
        statesets->downp2str = NULL;
    }
    if (statesets->upp2str) {
        free(statesets->upp2str);
        statesets->upp2str = NULL;
    }
    
    free(statesets);
}


int mpl_setup_statesets(Morphyp handl)
{
    MPL_ERR_T err = ERR_NO_ERROR;
    
    int numnodes = handl->numnodes;
    
    if (handl->statesets) {
        return 1;
    }
    
    handl->statesets = (MPLndsets**)calloc(numnodes, sizeof(MPLndsets*));
        
    if (!handl->statesets) {
        return ERR_BAD_MALLOC;
    }
    
    int i = 0;
    int nchars = mpl_get_num_charac((Morphyp)handl);
    
    for (i = 0; i < numnodes; ++i) {
        if (!(handl->statesets[i] = mpl_alloc_stateset(nchars))) {
            err = ERR_BAD_MALLOC;
            mpl_destroy_statesets(handl);
            break;
        }
    }
    
    return err; 
}


int mpl_destroy_statesets(Morphyp handl)
{
    int i = 0;
    int numnodes = handl->numnodes;
    
    if (handl->statesets) {
        
        for (i = 0; i < numnodes; ++i) {
            mpl_free_stateset(mpl_get_num_charac((Morphy)handl), handl->statesets[i]);
        }
        
        free(handl->statesets);
        handl->statesets = NULL;
    }
    
    return ERR_NO_ERROR;
}


int mpl_copy_data_into_tips(Morphyp handl)
{
    int i = 0;
    int j = 0;
    int ntax = mpl_get_numtaxa((Morphy)handl);
    int nchar = mpl_get_num_charac((Morphy)handl);
    MPLndsets** nsets = handl->statesets;
    
    for (i = 0; i < ntax; ++i) {
        for (j = 0; j < nchar; ++j) {
            nsets[i]->downpass1[j] =
            handl->inmatrix.cells[i * nchar + j].asint;
            nsets[i]->uppass1[j] = nsets[i]->downpass1[j];
            nsets[i]->uppass2[j] = nsets[i]->downpass1[j];
            nsets[i]->downpass2[j] = nsets[i]->downpass1[j];
        }
    }
    
    return ERR_NO_ERROR;
}

int mpl_assign_intwts_to_partitions(Morphyp handl)
{
    int i = 0;
    int j = 0;
    int numparts = mpl_get_numparts(handl);
    
    if (!numparts) {
        return ERR_NO_DATA;
    }
    
    for (i = 0; i < numparts; ++i) {
        
        if (handl->partitions[i]->intwts) {
            free(handl->partitions[i]->intwts);
            handl->partitions[i]->intwts = NULL;
        }
        
        handl->partitions[i]->intwts = (unsigned long*)calloc
                                        (handl->partitions[i]->ncharsinpart,
                                         sizeof(unsigned long));
        
        for (j = 0; j < handl->partitions[i]->ncharsinpart; ++j) {
            int charindex = handl->partitions[i]->charindices[j];
            handl->partitions[i]->intwts[j] = handl->charinfo[charindex].intwt;
        }
    }
    
    return 0;
}

int mpl_update_root(MPLndsets* lower, MPLndsets* upper, MPLpartition* part)
{
    int i = 0;
    int j = 0;
    int nchar = part->ncharsinpart;
    int *indices = part->charindices;
    

    for (i = 0; i < nchar; ++i) {
        j = indices[i];
        lower->downpass1[j] = upper->downpass1[j];
        lower->uppass1[j]  = upper->downpass1[j];
    }
    
    return 0;
}


int mpl_update_NA_root(MPLndsets* lower, MPLndsets* upper, MPLpartition* part)
{
    int i = 0;
    int j = 0;
    int nchar = part->ncharsinpart;
    int *indices = part->charindices;
    
    for (i = 0; i < nchar; ++i) {
        j = indices[i];
        
        if (upper->downpass1[j] & ISAPPLIC) {
            lower->downpass1[j] = upper->downpass1[j] & ISAPPLIC;
        }
        else {
            lower->downpass1[j] = NA;
        }
        
        lower->uppass2[j]   = upper->downpass2[j];
        lower->downpass1[j] = lower->downpass1[j];
        lower->uppass1[j]   = lower->downpass1[j];
        
        lower->temp_downpass1[j] = lower->downpass1[j];
        lower->temp_uppass1[j] = lower->uppass1[j];
        lower->temp_downpass2[j] = lower->downpass2[j];
        lower->temp_uppass2[j] = lower->uppass2[j];
    }
    
    return 0;
}

int mpl_update_NA_root_recalculation(MPLndsets* lower, MPLndsets* upper, MPLpartition* part)
{
    int i = 0;
    int j = 0;
    int nchar = part->nNAtoupdate;
    int *indices = part->update_NA_indices;
    
    for (i = 0; i < nchar; ++i) {
        j = indices[i];
        
        if (upper->downpass1[j] & ISAPPLIC) {
            lower->downpass1[j] = upper->downpass1[j] & ISAPPLIC;
        }
        else {
            lower->downpass1[j] = NA;
        }
        
        lower->uppass2[j]   = upper->downpass2[j];
        lower->downpass1[j] = lower->downpass1[j];
        lower->uppass1[j]   = lower->downpass1[j];
    }
    
    return 0;
}
