//
//  fitch.h
//  MorPhy2
//
//  Created by mbrazeau on 02/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//

#ifndef fitch_h
#define fitch_h

int mpl_fitch_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_fitch_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_local_reopt(MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part, int maxlen, bool domaxlen);

int mpl_NA_fitch_first_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_first_update_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_first_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_NA_fitch_first_update_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_NA_fitch_second_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_second_update_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_second_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_NA_fitch_second_update_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_NA_local_reopt (MPLndsets* srcset, MPLndsets* tgt1set, MPLndsets* tgt2set, MPLpartition* part, int maxlen, bool domaxlen);

int mpl_fitch_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_one_branch(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part);

int mpl_fitch_NA_first_one_branch(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part);

int mpl_fitch_NA_second_one_branch(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part);

int mpl_fitch_NA_second_one_branch_recalc(MPLndsets* tipanc, MPLndsets* node, MPLpartition* part);

int mpl_fitch_NA_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_NA_tip_recalc_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_NA_tip_finalize(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);
#endif /* fitch_h */
