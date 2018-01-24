	/*!
	 @file mpl.h
	 
	 @brief Defines the Morphy Phylogenetic Library API: a library for phylogenetic
	 computation accommodating morphological character hierarchies.
	 
	 Copyright (C) 2017  Martin D. Brazeau
	 
	 This program is free software: you can redistribute it and/or modify
	 it under the terms of the GNU General Public License as published by
	 the Free Software Foundation, either version 3 of the License, or
	 (at your option) any later version.
	 
	 This program is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.
	 
	 You should have received a copy of the GNU General Public License
	 along with this program.  If not, see <http://www.gnu.org/licenses/>.
	 
	 @discussion This header includes all the externally exported definitions and
	 function prototypes. A calling program creates an instance of a Morphy object 
	 and interacts with its elements through the functions described in this 
	 interface. The Morphy object contains no tree objects, but requires a 
	 pre-specified list of indices (integers) corresponding to the node indices in 
	 the calling program. Morphy will not keep track of the relationships between
	 the nodes, and it is up to the caller to keep track of these. Each character
	 *must* be assigned a type, and Morphy will make no default assumptions. Once 
	 one or more characters are assigned a function type (which creates internal
	 partitions), and a postorder list of nodes is known, then the library functions
	 can be called to reconstruct state sets and deliver length estimates for each
	 node.
	 
	 Morphy will provide functions for local reoptimisation, partial reoptimisation
	 and optimisation of subtrees.
	 
	 */

#ifndef mpl_h
#define mpl_h


#ifdef __cplusplus
	extern "C" {
#endif /*__cplusplus */
		
#include <stdbool.h>
#include "mplerror.h"
		
typedef void* Morphy;

typedef enum {
    
    NONE_T          = 0,
    FITCH_T         = 1,
    WAGNER_T        = 2,
    DOLLO_T         = 3,
    IRREVERSIBLE_T  = 4,
    USERTYPE_T      = 5,
    
    MAX_CTYPE,
    
} MPLchtype;

typedef enum {
    
    GAP_INAPPLIC,
    GAP_MISSING,
    GAP_NEWSTATE,
    
    GAP_MAX,
    
} MPLgap_t;

	// Public functions

	/*!
	 
	 @brief Creates a new instance of a Morphy object

	 @discussion Creates a new empty Morphy object. All fields are unpopulated and
	 uninitialised.
	 
	 @return A void pointer to the Morphy instance. NULL if unsuccessful.
	 
	 */
Morphy  mpl_new_Morphy
	
		(void);


	/*!
	 
	 @brief Destroys an instance of a Morphy object.

	 @discussion Destroys an instance of the Morphy object, calling all destructors
	 for internal object completely returning the memory to the system.
	 
 @param m A Morphy object to be destroyed.
 
 @return A Morphy error code.
 
 */
int     mpl_delete_Morphy
    
        (Morphy m);


/*!

 @brief Sets up the dimensions of the dataset.

 @discussion Provides initial dimensions for the dataset, which will constrain
 any input matrix supplied to Morphy.

 @param ntax The number of taxa (or tips/terminals).

 @param nchar The number of characters (i.e. transformation series) in the
 data set.
 
 @param m An instance of the Morphy object.

 @return Morphy error code.

 */
int     mpl_init_Morphy
    
        (const int  ntax,
         const int  nchar,
         Morphy     m);
    

/*!

 @brief Retrieve the number of taxa (rows) in the dataset.
 
 @discussion Retrieves the number of taxa (rows) in the dataset.

 @param m An instance of the Morphy object.

 @return The number of taxa if success, otherwise an error code.

 */
int     mpl_get_numtaxa
    
        (Morphy m);


/*!

 @brief Retrieve the number of taxa (rows) in the dataset.

 @discussion Retrieves the number of taxa (rows) in the dataset.

 @param m An instance of the Morphy object.

 @return The number of taxa if success, otherwise an error code.

 */
int     mpl_get_num_charac
    
        (Morphy     m);
    

/*!
 
 @brief Sets the number of internal nodes in the dataset
 
 @discussion This specifies the number of internal nodes over which
 reconstruction sets need to be made. It is up to the caller to ensure the 
 correct number of nodes and the relationships between them.
 
 @param nnodes The desired number of internal nodes.
 
 @param m An instance of the Morphy object
 
 @return A Morphy error code.
 
*/
int     mpl_set_num_internal_nodes
    
        (const int  nnodes,
         Morphy     m);

/*!
 
 @brief Gets the number of internal nodal reconstruction sets being used by
 MorphyLib.

 @discussion Gets the number of internal nodal reconstruction sets being used by
 MorphyLib.
 
 @param m An instance of the Morphy object.
 
 @return The number of internal nodes.
 
*/
int     mpl_get_num_internal_nodes
    
        (Morphy     m);

/*!

 @brief Attach a caller-specified list of symbols.

 @discussion Allows the caller to specify a list of symbols in the data matrix,
 otherwise, the symbols list used by Morphy will be extracted from the matrix. 
 The symbols list must match the symbols provided in the matrix. When Morphy 
 extracts symbols from the matrix, their ordering is alphanumeric, according to
 their ASCII codes (i.e. "+0123...ABCD...abcd..."). Loading a user-specified
 symbols list will override this ordering. Symbols loaded in either the list or
 the matrix must be valid Morphy character state symbols as defined in the 
 statedata.h header file.

 @param symbols A C-style (i.e. NULL-terminated) string of valid state symbols.

 @param m An instance of the Morphy object.

 @return Morphy error code.

 */
int     mpl_attach_symbols
    
        (const char*    symbols,
         Morphy         m);
    

/*!
 
 @brief Retrieves the current list of symbols.
 
 @discussion Returns a pointer to the string of character state symbols
 currently being used by Morphy (i.e. either the list of symbols extracted
 from the matrix, or the caller-specified values).
 
 @param m An instance of the Morphy object.
 
 @return A C-style (null-terminated) string of the character state symbols being
 used. NULL if failure.
 
 */
char*   mpl_get_symbols

        (const Morphy   m);


/*!

 @brief Attach raw character state data (i.e. tip data).
 

 @discussion Attaches a raw data character state matrix in the form of a C-style
 (i.e. NULL-terminated string). This can be the matrix block extracted from a
 Nexus file or an xread table format. The matrix should contain no terminal or
 tip labels.

 @param rawmatrix C-style string corresponding to the tip data.

 @param m An instance of the Morphy object.

 @return Morphy error code.

 */
int     mpl_attach_rawdata

        (const char*    rawmatrix,
         Morphy         m);


/*!
 
 @brief Deletes the caller-input data
 
 @discussion Deletes all of the user-input data and restores all parameters to
 their original values, except for the dimensions of the matrix.
 
 @param m An instance of the Morphy object.
 
 @return Morphy error code.
*/
int     mpl_delete_rawdata

        (Morphy     m);

    
int     mpl_set_gap_symbol

        (const char gapsymb,
         Morphy     m);


int     mpl_set_missing_symbol
    
        (const char missymb,
         Morphy     m);


/*!
 
 @brief Commits parameters prior to nodal set calculations.

 @discussion Once the caller is satisfied with the setup of types, weights, and
 partitioning, this function must be called, thereby committing the parameters
 until any changes are made. If no character types have been assigned, the 
 function will fail with an error code.
 
 @param m An instance of the Morphy object.
 
 @return A Morphy error code.
*/
int     mpl_apply_tipdata
    
        (Morphy m);


int     mpl_incl_charac

        (const int  charID,
         Morphy     m);


int     mpl_excl_charac

        (const int  charID,
         Morphy     m);


/*!
    
 @brief Sets the weight for a specified character.
 
 @discussion Sets a weight for a specified character. The function takes a 
 floating point value. However, in the current implementation, fractional values 
 will be interpreted and estimated using an approximation of the rational 
 factors. In the current version, MorphyLib uses this method to circumvent any
 floating point calculations that aren't absolutely necessary.
 
 @param charID The index of the character to be weighted.
 
 @param weight The weight requested for the character.
 
 @param m An instance of the Morphy object.
 
 @return A Morphy error code.
 
 */
int     mpl_set_charac_weight

        (const int      charID,
         const double   weight,
         Morphy         m);

// TODO: Document
unsigned long mpl_get_charac_weight
    
        (double*    weight,
         const int char_id,
         const Morphy   m);

/*!

 @brief Sets a character's parsimony function type
 
 @discussion Set the parsimony function type to one defined in the morphydefs.h
 header file. Setting the character to type NONE_T will also cause it to be 
 excluded from any further calculations.
 
 @param charID The index of the character (transformation series) as defined in
 the input matrix.
 
 @param chtype The parsimony function type as defined in morphydefs.h
 
 @param m An instance of the Morphy object.
 
 @return A Morphy error code.

 */
int     mpl_set_parsim_t

        (const int charID,
         const MPLchtype chtype,
         Morphy           m);


/*!

 @brief Tells MorphyLib how to treat the gap symbol.

 @discussion The caller can specify the type of gap handling to use before the
 tipdata are applied. The options are documented in the morphydefs.h header file
 but include at least GAP_INAPPLIC, GAP_MISSING, and GAP_NEWSTATE to specify
 inapplicable, missing, or new state values respectively. These values are
 applied to all characters for which they are appropriate.
 
 @param gaptype The type of gap treatment to be applied (documented in
 morphydefs.h).
 
 @param m An instance of the Morphy object.

 @return A Morphy error code.
 
*/
int     mpl_set_gaphandl

        (const MPLgap_t gaptype,
         Morphy      m);

/*!
 
 @brief Returns the type of gap handling method currently in effect.
 
 @discussion Returns the type of gap handling method currently in effect. The 
 methods are defined in the morphydefs.h file.
 
 @param m An instance of the Morphy object.
 
 @return A Morphy error code.
 
*/
int     mpl_query_gaphandl

        (Morphy     m);


/*!
 
 @brief Reconstructs the first (downpass) nodal reconstructions
 
 @discussion Reconstructs the preliminary nodal set for all characters for a 
 particular node. This function is called over a postorder sequence of internal 
 nodes where left and right descendants are known.
 
 Because this function needs to be fairly high-performance, it does not do much 
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate 
 parameters have been set before use.
 
 @param node_id The index of the node being reconstructed.
 
 @param left_id The index of the left descendant.
 
 @param right_id The index of the right descendant.
 
 @param m An instance of the Morphy object.
 
 @return The integral parsimony length (right now)
 
 */
int     mpl_first_down_recon

        (const int  node_id,
         const int  left_id,
         const int  right_id,
         Morphy     m);

    
/*!

 @brief Reconstructs the second (uppass) nodal reconstructions.
 
 @discussion Reconstructs second-pass nodal sets. For normal (all-applicable)
 characters, this is the final pass. This function is called over a preorder
 sequence of nodes where left, right, and ancestral nodes are known.
 
 Because this function needs to be fairly high-performance, it does not do much
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate
 parameters have been set before use.
 
 @param node_id The index of the node being reconstructed.
 
 @param left_id The index of the left descendant.
 
 @param right_id The index of the right descendant.

 @param anc_id The index of the immediate ancestor of the node.
 
 @param m An instance of the Morphy object.
 
 @return A null value (for now).
 */
int     mpl_first_up_recon
        (const int node_id,
         const int left_id,
         const int right_id,
         const int anc_id,
         Morphy m);

    
/*!
 
 @brief Performs the second nodal reconstructions for characters with
 inapplicability.
 
 @discussion Updates the nodal sets that had ambiguous unions with the 
 inapplicable state and calculates steps involving applicable states after 
 the update.
 
 Because this function needs to be fairly high-performance, it does not do much
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate
 parameters have been set before use.
 
 @param node_id The index of the node being reconstructed.
 
 @param left_id The index of the left descendant.
 
 @param right_id The index of the right descendant.
 
 @param m An instance of the Morphy object.
 
 @return The integral parsimony length (right now)
 */
int     mpl_second_down_recon
    
        (const int  node_id,
         const int  left_id,
         const int  right_id,
         Morphy     m);

    
/*!
 
 @brief Finalises the ancestral state reconstructions for characters with 
 inapplicable values.
 
 @discussion Finalises the nodal sets for any characters that may have involved
 the inapplicable token and counts excess regions of applicability at nodes
 having at least two descendant subtrees that possess any applicable characters.
 
 Because this function needs to be fairly high-performance, it does not do much
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate
 parameters have been set before use.
 
 @param node_id The index of the node being reconstructed.
 
 @param left_id The index of the left descendant.
 
 @param right_id The index of the right descendant.
 
 @param anc_id The index of the immediate ancestor of the node.
 
 @param m An instance of the Morphy object.
 
 @return The integral parsimony length (for now)
 */
int     mpl_second_up_recon

        (const int node_id,
         const int left_id,
         const int right_id,
         const int anc_id,
         Morphy m);

    
/*!
 
 @brief Initial update of tip values following uppass reconstruction
 
 @discussion Polymorphic terminal state sets need to be resolved after the 
 uppass based on descendant state values in order for local reoptimisation 
 procedures to be accurate and for inapplicable step counting to proceed 
 accurately. This function calls updaters for the records of states active on 
 the subtrees, thereby allowing the second downpass to accurately reconstruct 
 subtree state activity. Missing values are left as-is in characters with 
 inapplicability, otherwise, final ancestral state reconstructions may be 
 inaccurate.
 
 Because this function needs to be fairly high-performance, it does not do much
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate
 parameters have been set before use.
 
 @param tip_id The index of the tip being updated.
 
 @param anc_id The index of the tip's immediate ancestor.
 
 @param m An instance of the Morphy object.
 
 @return A null value (for now).
 */
int     mpl_update_tip
    
        (const int tip_id,
         const int anc_id,
         Morphy m);

/*!
 
 @brief Finalizes ambiguous or missing values in the tips. 
 
 @discussion Ambiguous terminal state sets need to be resolved after the uppass
 based on descendant state values in order for local reoptimisation procedures 
 to be accurate and for inapplicable step counting to proceed accurately. This
 function calls updaters for the records of states active on the subtrees, 
 thereby allowing local reoptimization functions to accurately predict length
 increases when a subtree is added near a tip.
 
 Because this function needs to be fairly high-performance, it does not do much
 checking for parameter validity, thus unsafe usage of this function might not
 be caught. It is up to calling functions to ensure that the appropriate
 parameters have been set before use.
 
 @param tip_id The index of the tip being updated.
 
 @param anc_id The index of the tip's immediate ancestor.
 
 @param m An instance of the Morphy object.
 
 @return A null value (for now).
 */
int     mpl_finalize_tip

        (const int  tip_id,
         const int  anc_id,
         Morphy     m);

/*!
 @brief Used to update a root-like tip in an unrooted tree and length added.
 
 @discussion If using an unrooted tree structure, a tip is commonly used as an
 entry point for traversals on the tree. This tip is jointed either as an extra
 descendant or the ancestor of the calculation root node in the tree. In these
 circumstances, a binary traversal on the tree will not give complete 
 reconstructions or length counts for the tree. This function is called at the 
 end of the optimization process and is required for the complete tree length of
 an unrooted tree. This function will update the state sets of both nodes
 reciprocally.
 
 @param tip_id An index corresponding to the tip number being updated.
 
 @param node_id An index of the tip's neighboring internal node.
 
 @param m An instance of the Morphylib object.
 
 @return The weighted number of steps (positive) or a negative number 
 corresponding to a morphylib error code.
 */
        
int     mpl_do_tiproot
        
        (const int  tip_id,
         const int  node_id,
         Morphy     m);
        
int     mpl_finalize_tiproot
        
        (const int  tip_id,
         const int  node_id,
         Morphy     m);
/*!
 
 @brief Updates the nodal sets for a lower ('dummy') root node
 
 @discussion If trees are rooted, then Morphy uppass functions
 require a lower or 'dummy' root in order to function properly. This
 function should be called to set the nodal state sets to the dummy
 root. The nodal set will be equal to the set of the root node, unless
 there is an ambiguous union of applicable and gap tokens when gaps are 
 treated as in applicable. In which case, the set union is resolved in 
 favour of any applicable tokens in the set.
 
 @param l_root_id The index of the lower root.
 
 @param root_id The index of the upper root node.
 
 @return A Morphy error code.
 
 */
int     mpl_update_lower_root
    
        (const int  l_root_id,
         const int  root_id,
         Morphy     m);
    
    
int		mpl_na_first_down_recalculation

		(const int  node_id,
		 const int  left_id,
		 const int  right_id,
		 Morphy     m);
    

int		mpl_na_first_up_recalculation

		(const int  node_id,
		 const int  left_id,
		 const int  right_id,
		 const int  anc_id,
		 Morphy     m);


// Returns number of steps to add
int		mpl_na_second_down_recalculation

		(const int  node_id,
		 const int  left_id,
		 const int  right_id,
		 Morphy     m);

// Returns number of steps to add
int		mpl_na_second_up_recalculation

		(const int  node_id,
		 const int  left_id,
		 const int  right_id,
		 const int  anc_id,
		 Morphy     m);
        
int     mpl_lower_root_recalculation
        
        (const int  l_root_id,
         const int  root_id,
         Morphy     m);

int     mpl_na_tiproot_recalculation
        
        (const int  tip_id,
         const int  node_id,
         Morphy     m);

int     mpl_na_tiproot_final_recalculation
        
        (const int  tip_id,
         const int  node_id,
         Morphy     m);
        
int     mpl_get_insertcost

        (const int  srcID,
         const int  tgt1ID,
         const int  tgt2ID,
         const bool max,
         const int  cutoff,
         Morphy     m);
        
int     mpl_na_update_tip
        
        (const int  tip_id,
         const int  anc_id,
         Morphy     m);

    
int     mpl_get_step_recall
        
        (const int  node_id,
         Morphy     m);
        
        
// Indicates whether or not partitions with inapplicable characters need partial
// reoptimisation on the target subtree. SHOULD RETURN: Number of characters
// needing partial reoptimisation on the subtree.
int     mpl_check_reopt_inapplics
    
        (Morphy m);
        
bool    mpl_check_updated
        
        (const int  node_id,
         Morphy     m);
        
/*!
 
 @brief Restores original state sets at a node.
 
 @discussion This function restores the state sets in the tree to the ones 
 calculated using the initial fullpass optimizations (first, and second down and
 up functions). It is used after partial reoptimization of a tree and if the
 client program needs to restore the state sets to their original values before
 continuing to evaluate proposed insertions.
 @param node_id The index value of the node having its state sets restored.
 
 @param m An instance of the morphylib object.
 
 @return 0 if success, a morphylib error code if there has been an error.
 
 */
int     mpl_restore_original_sets
        
        (const int  node_id,
         Morphy     m);
/*!

 @brief Returns the state set for a character at a given node as set bits in an
 unsigned integer.
 
 @discussion If the caller requires the internal state representation of a nodal
 set used by MorphyLib, this function can be called to retrieve it. The caller
 needs to specify the node index, the character number, and the pass number 
 (1-based, because these are not indices in a C array).

 @param nodeID The index of the node set required.

 @param character The character number to be queried.

 @param passnum The traversal iteration corresponding to the set required. These
 range from 1 to 4 and represent first downpass, first uppass, second downpass
 and second uppass respectively.

 @return An unsigned integer with bits set corresponding to values used by
 MorphyLib.
 
 */
unsigned
int     mpl_get_packed_states

        (const int  nodeID,
         const int  character,
         const int  passnum,
         Morphy     m);


const char*   mpl_get_stateset
    
        (const int  nodeID,
         const int  character,
         const int  passnum,
         Morphy     m);

#ifdef __cplusplus
}
#endif /*__cplusplus */

#endif /* mpl_h */
