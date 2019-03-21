/*
//  wagner.c
//  morphylib
//
//  Created by mbrazeau on 21/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
*/
#include "mpl.h"
#include "morphydefs.h"
#include "morphy.h"
#include "wagner.h"

static inline unsigned mpl_closed_interval(MPLstate* res, MPLstate a, MPLstate b)
{
    unsigned steps = 0;
    MPLstate c = 0;
    
    if (b > a) {
        c = b;
        b = a;
        a = c;
    }
    
    *res = a & (-a);
    
    while(!(*res & b)) {
        ++steps;
        *res |= a >> steps;
    }
    
    return steps;
}

int mpl_wagner_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int steps = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* n     = nset->downpass1;
    
    unsigned long* weights = part->intwts;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if (left[j] & right[j]) {
            n[j] = left[j] & right[j];
        }
        else {
            
            n[j] = 0;
            steps += weights[i] * mpl_closed_interval(&n[j], left[j], right[j]);
        }
    }
    
    return steps;
}

int mpl_wagner_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    MPLstate* left  = lset->downpass1;
    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nfin  = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if ((anc[j] & npre[j]) == anc[j]) {
            nfin[j] = anc[j] & npre[j];
        }
        else {
            MPLstate res = 0;
            mpl_closed_interval(&res, left[j], right[j]);
            nfin[j] = (res & anc[j]) | npre[j];
        }
       
        assert(nfin[j]);
    }
    
    return 0;
}

int mpl_wagner_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    
    MPLstate* tprelim = tset->downpass1;
    MPLstate* tfinal  = tset->uppass1;
    MPLstate* astates = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        j = indices[i];
        if (tprelim[j] & astates[j]) {
            tfinal[j] = tprelim[j] & astates[j];
        }
        else {
            tfinal[j] = tprelim[j];
        }
        assert(tfinal[j]);
    }
    return 0;
}
