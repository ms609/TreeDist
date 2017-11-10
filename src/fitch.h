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

int mpl_NA_fitch_first_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_first_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_NA_fitch_second_downpass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);

int mpl_NA_fitch_second_uppass(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_NA_tip_update(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);

int mpl_fitch_NA_tip_finalize(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);
#endif /* fitch_h */
