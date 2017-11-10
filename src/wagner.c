//
//  wagner.c
//  morphylib
//
//  Created by mbrazeau on 21/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//
#include "morphydefs.h"
#include "morphy.h"
#include "wagner.h"



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
            MPLstate temp = 0;
            MPLstate min = 0;
            MPLstate max = 0;
            if (left[j] > right[j]) {
                max = left[j];
                min = right[j];
            }
            else {
                min = left[j];
                max = right[j];
            }
            
            temp = max & ~(max-1);
            n[j] = temp;
            while (!(n[j] & min)) {
                steps += weights[i];
                n[j] |= temp >> steps;
            }
            // TODO: Multiply by weights BEFORE next j
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
//    MPLstate* left  = lset->downpass1;
//    MPLstate* right = rset->downpass1;
    MPLstate* npre  = nset->downpass1;
    MPLstate* nfin  = nset->uppass1;
    MPLstate* anc   = ancset->uppass1;
    
    for (i = 0; i < nchars; ++i) {
        
        j = indices[i];
        
        if ((anc[j] & npre[j]) == anc[j]) {
            nfin[j] = anc[j] & npre[j];
        }
        else {
//            if (left[j] & right[j]) {
//                nfin[j] = (npre[j] | (anc[j] & (left[j] | right[j])));
//            }
//            else {
//                nfin[j] = npre[j] | anc[j];
//            }
        }
       
        assert(nfin[j]);
    }
    
    return 0;
}

int mpl_wagner_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part)
{
    // TODO: Must change this to a Wagner-specific option
    int i     = 0;
    int j     = 0;
    int* indices    = part->charindices;
    int nchars      = part->ncharsinpart;
    // TODO: Check these!!!!!!!
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
