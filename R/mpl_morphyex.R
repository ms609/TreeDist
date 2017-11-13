#' Converts a numeric error code to human-readable format
#' 
#' @param errorCode Non-positive integer to be converted
#'
#' @return A character string corresponding to the provided error code
#'
#' @examples mpl_translate_error(-1) # "ERR_INVALID_SYMBOL"
#'
#' @author Martin R. Smith
#' @keywords internal
#' @export

mpl_translate_error <- function (errorCode) {
  mplErrorCodes <- rev(c(
      "ERR_EX_DATA_CONF",
      "ERR_OUT_OF_BOUNDS",
      "ERR_CASE_NOT_IMPL",
      "ERR_UNKNOWN_CHTYPE",
      "ERR_SYMBOL_MISMATCH",
      "ERR_MATCHING_PARENTHS",
      "ERR_ATTEMPT_OVERWRITE",
      "ERR_NO_DIMENSIONS",
      "ERR_DIMENS_UNDER",
      "ERR_DIMENS_OVER",
      "ERR_NO_DATA",
      "ERR_BAD_MALLOC",
      "ERR_BAD_PARAM",
      "ERR_UNEXP_NULLPTR",
      "ERR_INVALID_SYMBOL",
      "ERR_NO_ERROR"))
  return (mplErrorCodes[1-errorCode])
}

#' Creates a new instance of a Morphy object
#'
#' Creates a new empty Morphy object. All fields are unpopulated
#' and uninitialised.
#'
#' 
#' @return A void pointer to the Morphy instance. NULL if unsuccessful.
#' 
#' @examples morphyObj <- mpl_new_Morphy() # Create new object
#' ## Do some stuff ... ##
#' mpl_delete_Morphy(morphyObj) # Delete when done
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_new_Morphy <- function() {
    .Call("_R_wrap_mpl_new_Morphy")
}

#' Destroys an instance of a Morphy object.
#'
#' Destroys an instance of the Morphy object, calling all
#' destructor for internal object completely returning the memory to the system.
#'
#' @param morphyobj A Morphy object to be destroyed.
#' 
#' @return A Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_delete_Morphy <- function(morphyobj) {
    .Call("_R_wrap_mpl_delete_Morphy", morphyobj)
}

#' Sets up the dimensions of the dataset.
#'
#' Provides initial dimensions for the dataset, which will
#' constrain any input matrix supplied to Morphy.
#' 
#' @param morphyobj An instance of the Morphy object.
#' @param ntax The number of taxa (or tips/terminals).
#' @param nchar The number of characters (i.e. transformation series) in the
#' data set.
#' 
#' @return Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_init_Morphy <- function(numtaxa, numchars, morphyobj) {
    .Call("_R_wrap_mpl_init_Morphy", as.integer(numtaxa), as.integer(numchars), morphyobj)
}


#' Retrieve the number of taxa (rows) in the dataset.
#'
#' Retrieves the number of taxa (rows) in the dataset.
#'
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The number of taxa if success, otherwise an error code.
#'
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_get_numtaxa <- function(morphyobj) {
    .Call("_R_wrap_mpl_get_numtaxa", morphyobj)
}


#' Set the weight of a character in the dataset
#'
#' Sets the weight of a character in the dataset.
#'
#' @param charID    Number of the character (i.e. first character is number 1)
#' @param weight    Weight to assign
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return An error code.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
mpl_set_charac_weight <- function (charID, weight, morphyobj) {
  return (.Call('_R_wrap_mpl_set_charac_weight', as.integer(charID - 1L), as.numeric(weight),
                morphyobj))
}

#' Retrieve the weight of a character in the dataset
#'
#' Gets the weights of a character in the dataset.
#'
#' @param charID    Number of the character (i.e. first character is number 1)
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A list, detailing (item 1) the exact weight of the character; (item 2) the integer
#'         approximation used by Morphy.
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
mpl_get_charac_weight <- function (charID, morphyobj) {
  return (.Call('_R_wrap_mpl_get_charac_weight', as.integer(charID - 1L), 
                morphyobj))
}


#' Retrieve the number of character (columns) in the dataset.
#'
#' Retrieves the number of character (columns) in the dataset.
#'
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The number of internal nodes.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_get_num_charac <- function(morphyobj) {
    .Call("_R_wrap_mpl_get_num_charac", morphyobj)
}


#' Attach a caller-specified list of symbols.
#'
#' Allows the caller to specify a list of symbols in the data matrix,
#' otherwise, the symbols list used by Morphy will be extracted from the matrix. 
#' The symbols list must match the symbols provided in the matrix. When Morphy 
#' extracts symbols from the matrix, their ordering is alphanumeric, according to
#' their ASCII codes (i.e. "+0123...ABCD...abcd..."). Loading a user-specified
#' symbols list will override this ordering. Symbols loaded in either the list or
#' the matrix must be valid Morphy character state symbols as defined in the 
#' statedata.h header file.  The list must end with a semicolon.
#'
#' @param symbols A C-style (i.e. NULL-terminated) string of valid state symbols.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_attach_symbols <- function(symbols, morphyobj) {
    .Call("_R_wrap_mpl_attach_symbols", symbols, morphyobj)
}


#' Attach raw character state data (i.e. tip data).
#'
#' Attaches a raw data character state matrix in the form of a C-style
#' (i.e. NULL-terminated string). This can be the matrix block extracted from a
#' Nexus file or an xread table format. The matrix should contain no terminal or
#' tip labels.
#'
#' @param rawmatrix C-style string corresponding to the tip data.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_attach_rawdata <- function(rawdata, morphyobj) {
    .Call("_R_wrap_mpl_attach_rawdata", rawdata, morphyobj)
}

#' Retrieves the current list of symbols.
#'
#' Returns a pointer to the string of character state symbols
#' currently being used by Morphy (i.e. either the list of symbols extracted
#' from the matrix, or the caller-specified values).
#'
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A C-style (null-terminated) string of the character state symbols
#' being used. NULL if failure.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export

mpl_get_symbols <- function(morphyobj) {
  .Call("_R_wrap_mpl_get_symbols", morphyobj)  
}


#' Sets a character's parsimony function type
#'
#' Set the parsimony function type to one defined in the
#' morphydefs.h header file. Setting the character to type NONE_T will also
#' cause it to be excluded from any further calculations.
#'
#' @param char_id   The number of the character (transformation series) as defined
#'                  in the input matrix.  The first character is numbered 1 (one).
#' @param tname     The parsimony function type as defined in morphydefs.h
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_set_parsim_t <- function(char_id, tname = "typename", morphyobj) {
    .Call("_R_wrap_mpl_set_parsim_t", as.integer(char_id - 1L), tname, morphyobj)
}


#' Sets the number of internal nodes in the dataset
#'
#' This specifies the number of internal nodes over which
#' reconstruction sets need to be made. It is up to the caller to ensure the 
#' correct number of nodes and the relationships between them.
#'
#' @param nnodes The desired number of internal nodes.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_set_num_internal_nodes <- function(numnodes, morphyobj) {
    .Call("_R_wrap_mpl_set_num_internal_nodes", as.integer(numnodes), morphyobj)
}



#' Gets the number of internal nodal reconstruction sets being used by
#' MorphyLib.
#'
#' Gets the number of internal nodal reconstruction sets being used
#' by MorphyLib.
#'
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The number of internal nodes.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export

mpl_get_num_internal_nodes <- function(morphyobj) {
    .Call("_R_wrap_mpl_get_num_internal_nodes", morphyobj)
}


#' Commits parameters prior to nodal set calculations.
#'
#' Once the caller is satisfied with the setup of types, weights,
#' and partitioning, this function must be called, thereby committing the
#' parameters until any changes are made. If no character types have been
#' assigned, the function will fail with an error code.
#'
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A Morphy error code.
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_apply_tipdata <- function(morphyobj) {
    .Call("_R_wrap_mpl_apply_tipdata", morphyobj)
}

#' Reconstructs the first (downpass) nodal reconstructions
#'
#' Reconstructs the preliminary nodal set for all characters for a 
#' particular node. This function is called over a postorder sequence of internal 
#' nodes where left and right descendants are known.
#' Because this function needs to be fairly high-performance, it does not do much 
#' checking for parameter validity, thus unsafe usage of this function might not
#' be caught. It is up to calling functions to ensure that the appropriate 
#' parameters have been set before use.
#' 
#' @param node_id The index of the node being reconstructed.
#' @param left_id The index of the left descendant.
#' @param right_id The index of the right descendant.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The integral parsimony length (right now)
#' 
#' @author Martin Brazeau
#' @keywords internal
#' @export
mpl_first_down_recon <- function(node_id, left_id, right_id, morphyobj) {
    .Call("_R_wrap_mpl_first_down_recon", as.integer(node_id), as.integer(left_id), as.integer(right_id), morphyobj)
}


#' Deletes the caller-input data.
#'
#' Deletes all of the user-input data and restores all parameters
#' to their original values, except for the dimensions of the matrix.
#' 
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return Morphy error code.
#' 
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_delete_rawdata <- function(morphyobj) {
    .Call("_R_wrap_mpl_delete_rawdata", morphyobj)
}


#' Reconstructs the second (uppass) nodal reconstructions.
#'
#' Reconstructs second-pass nodal sets. For normal (all-applicable)
#' characters, this is the final pass. This function is called over a preorder
#' sequence of nodes where left, right, and ancestral nodes are known.
#' Because this function needs to be fairly high-performance, it does not do much 
#' checking for parameter validity, thus unsafe usage of this function might not
#' be caught. It is up to calling functions to ensure that the appropriate 
#' parameters have been set before use.
#' 
#' @param node_id The index of the node being reconstructed.
#' @param left_id The index of the left descendant.
#' @param right_id The index of the right descendant.
#' @param anc_id The index of the immediate ancestor of the node. 
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A null value (for now).
#' 
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_first_up_recon <- function(node_id, left_id, right_id, anc_id, morphyobj) {
    .Call("_R_wrap_mpl_first_up_recon", as.integer(node_id), as.integer(left_id), 
                 as.integer(right_id), as.integer(anc_id), morphyobj)
}


#' Performs the second nodal reconstructions for characters with
#' inapplicability.
#'
#' Updates the nodal sets that had ambiguous unions with the 
#' inapplicable state and calculates steps involving applicable states after 
#' the update.
#' Because this function needs to be fairly high-performance, it does not do much 
#' checking for parameter validity, thus unsafe usage of this function might not
#' be caught. It is up to calling functions to ensure that the appropriate 
#' parameters have been set before use.
#' 
#' @param node_id The index of the node being reconstructed.
#' @param left_id The index of the left descendant.
#' @param right_id The index of the right descendant.
#' @param anc_id The index of the immediate ancestor of the node. 
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The integral parsimony length (right now)
#' 
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_second_down_recon <- function(node_id, left_id, right_id, morphyobj) {
    .Call("_R_wrap_mpl_second_down_recon", as.integer(node_id), as.integer(left_id),
                 as.integer(right_id), morphyobj)
}


#' Finalises the ancestral state reconstructions for characters with 
#' inapplicable values.
#'
#' Finalises the nodal sets for any characters that may have involved
#' the inapplicable token and counts excess regions of applicability at nodes
#' having at least two descendant subtrees that possess any applicable characters.
#' Because this function needs to be fairly high-performance, it does not do much 
#' checking for parameter validity, thus unsafe usage of this function might not
#' be caught. It is up to calling functions to ensure that the appropriate 
#' parameters have been set before use.
#' 
#' @param node_id The index of the node being reconstructed.
#' @param left_id The index of the left descendant.
#' @param right_id The index of the right descendant.
#' @param anc_id The index of the immediate ancestor of the node. 
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The integral parsimony length (right now)
#'  
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_second_up_recon <- function(node_id, left_id, right_id, anc_id, morphyobj) {
    .Call("_R_wrap_mpl_second_up_recon", as.integer(node_id), as.integer(left_id), 
                 as.integer(right_id), as.integer(anc_id), morphyobj)
}

#' Initial update of tip values following uppass reconstruction.
#'
#' Ambiguous terminal state sets need to be resolved after the first uppass
#' based on descendant state values in order for local reoptimisation procedures 
#' to be accurate and for inapplicable step counting to proceed accurately. This
#' function calls updaters for the records of states active on the subtrees, 
#' thereby allowing the second downpass to accurately reconstruct subtree state
#' activity.
#' Because this function needs to be fairly high-performance, it does not do much 
#' checking for parameter validity, thus unsafe usage of this function might not
#' be caught. It is up to calling functions to ensure that the appropriate 
#' parameters have been set before use.
#' 
#' @param tip_id The index of the tip being updated.
#' @param anc_id The index of the tip's immediate ancestor.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return The integral parsimony length (right now)
#' 
#' @seealso A null value (for now).
#' 
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_update_tip <- function(tip_id, anc_id, morphyobj) {
    .Call("_R_wrap_mpl_update_tip", as.integer(tip_id), as.integer(anc_id), morphyobj)
}

#' Updates the nodal sets for a lower ('dummy') root node
#'
#' If trees are rooted, then Morphy uppass functions
#' require a lower or 'dummy' root in order to function properly. This
#' function should be called to set the nodal state sets to the dummy
#' root. The nodal set will be equal to the set of the root node, unless
#' there is an ambiguous union of applicable and gap tokens when gaps are 
#' treated as in applicable. In which case, the set union is resolved in 
#' favour of any applicable tokens in the set.
#' 
#' @param l_root_id The index of the lower root.
#' @param root_id The index of the upper root node.
#' @param morphyobj An instance of the Morphy object.
#' 
#' @return A Morphy error code.
#' 
#' @author Thomas Guillerme
#' @keywords internal
#' @export
mpl_update_lower_root <- function(l_root_id, root_id, morphyobj) {
    .Call("_R_wrap_mpl_update_lower_root", as.integer(l_root_id), as.integer(root_id),
                 morphyobj)
}