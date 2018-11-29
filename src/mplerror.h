/*!
 @file mplerror.h
 @brief Error codes and descriptions for MorphyLib.
 */

#ifndef mplerror_h
#define mplerror_h

/*!
 @typedef MPL_ERR_T
 @brief List of error codes. Each error is a negative value.
 @discussion These error codes are returned by library interface functions (and
 used by some internal functions) to report errors back to the caller.
 */
typedef enum {
    
    ERR_EX_DATA_CONF        = -15,  /*! Input conflicts with existing dataset */
    ERR_OUT_OF_BOUNDS       = -14,  /*! Attempt to index out of bounds of an 
                                        array */
    
    ERR_CASE_NOT_IMPL       = -13,  /*! Case not implemented. */
    
    ERR_UNKNOWN_CHTYPE      = -12,  /*! Character type is unknown. It either 
                                        exceeds the list of character types or
                                        a user type matrix has not yet been 
                                        supplied. */
    
    ERR_SYMBOL_MISMATCH     = -11,  /*! Symbols list and matrix have a mismatch
                                       (i.e. symbol not found).*/
    ERR_MATCHING_PARENTHS   = -10,  /*! Data input has unexpected non-matching 
                                        parentheses.*/
    ERR_ATTEMPT_OVERWRITE   = -9,   /*! Caller attempted to overwrite a loaded 
                                        dataset.*/
    ERR_NO_DIMENSIONS       = -8,   /*! Function requires pre-specified 
                                        dimensions to function properly.*/
    ERR_DIMENS_UNDER        = -7,   /*! Supplied dimensions underestimate size 
                                        of dataset.*/
    ERR_DIMENS_OVER         = -6,   /*! Supplied dimensions overestimate size of 
                                        dataset.*/
    ERR_NO_DATA             = -5,   /*! No dataset supplied.*/
    
    ERR_BAD_MALLOC          = -4,   /*! Memory allocation failure.*/
    
    ERR_BAD_PARAM           = -3,   /*! Unexpected parameter value passed to 
                                        function.*/
    ERR_UNEXP_NULLPTR       = -2,   /*! Unexpected NULL pointer passed to 
                                        function.*/
    ERR_INVALID_SYMBOL      = -1,   /*! Symbol in dataset or symbol list is not 
                                        allowed by Morphy.*/
    ERR_NO_ERROR            =  0    /*! No error. Everything went OK.*/
    
} MPL_ERR_T;

#endif /* mplerror_h */
