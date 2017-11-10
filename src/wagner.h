//
//  wagner.h
//  morphylib
//
//  Created by mbrazeau on 21/05/2017.
//  Copyright © 2017 brazeaulab. All rights reserved.
//

#ifndef wagner_h
#define wagner_h

int mpl_wagner_downpass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLpartition* part);
int mpl_wagner_uppass
(MPLndsets* lset, MPLndsets* rset, MPLndsets* nset, MPLndsets* ancset,
 MPLpartition* part);
int mpl_wagner_tip_update
(MPLndsets* tset, MPLndsets* ancset, MPLpartition* part);
#endif /* wagner_h */
