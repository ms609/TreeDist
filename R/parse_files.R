#' Parse TNT Tree
#' 
#' Reads a tree from TNT's paranthetical output.
#' 
#' @param filename character string specifying path to TNT `.tre` file
#' @param relativePath (optional) character string specifying location of the matrix 
#'                     file used to generate the TNT results, 
#'                     relative to the current working directory, for portability.
#' @param keepEnd (optional, default 1) integer specifying how many elements of the file
#'                path to conserve when creating relative path (see examples)
#'                     
#' 
#' @return a tree of class \code{phylo}.
#' 
#' @examples {
#'   \dontrun{
#'   # TNT read a matrix from c:/myproject/tnt/coding1/dataset.nex
#'   # The results of an analysis were written to c:/myproject/tnt/output/results1.tnt
#'   # results1.tnt will contain a hard-coded reference to 
#'   # "c:/myproject/tnt/coding1/dataset.nex"
#'   
#'   getwd() # Gives the current working directory
#'   
#'   # Say that working directory is c:/myproject, which perhaps corresponds to a
#'   # Git repository.
#'   # This directory may be saved into another location by collaborators, or on a 
#'   # different filesystem by a continuous integration platform.
#'   
#'   # Works on local machine but not elsewhere:
#'   ReadTntTree('tnt/output/results1.tnt')
#'   
#'   # Takes only the filename from the results
#'   ReadTntTree('tnt/output.results1.tnt', 'tnt/coding1')
#'   
#'   # Uses the last three elements of c:/myproject/tnt/coding1/dataset.nex
#'   #                                               3     2       1
#'   # '.' means "relative to the current directory", which is c:/myproject
#'   ReadTntTree('tnt/output/results1.tnt', '.', 3)
#'   
#'   # If the current working directory was c:/myproject/rscripts/testing,
#'   # you could navigate up the directory path with '..':
#'   ReadTntTree('../../tnt/output/results1.tnt', '../..', 3)
#'   
#'   }
#' }
#' 
#' @author Martin R. Smith
#' @importFrom ape read.tree
#' @export
ReadTntTree <- function (filename, relativePath = NULL, keepEnd = 1L) {
  fileText <- readLines(filename)
  trees <- lapply(fileText[2:(length(fileText)-1)], TNTText2Tree)
  
  taxonFile <- gsub("tread 'tree(s) from TNT, for data in ", '', fileText[1], fixed=TRUE)
  taxonFile <- gsub("'", '', gsub('\\', '/', taxonFile, fixed=TRUE), fixed=TRUE)
  if (!is.null(relativePath)) {
    taxonFileParts <- strsplit(taxonFile, '/')[[1]]
    nParts <- length(taxonFileParts)
    if (nParts < keepEnd) {
      stop("Taxon file path (", taxonFile, ") contains fewer than keepEnd (",
           keepEnd, ") components.")
    }
    taxonFile <- paste0(c(relativePath, taxonFileParts[(nParts + 1L - keepEnd):nParts]),
                        collapse='/')
  }
  if (!file.exists(taxonFile)) {
    warning("Cannot find linked data file ", taxonFile)
  } else {
    tipNames <- rownames(ReadTntCharacters(taxonFile, 1))
    trees <- lapply(trees, function (tree) {
      tree$tip.label <- tipNames[as.integer(tree$tip.label) + 1]
      tree
    })
  }
  
  # Return:
  if (length(trees) == 1) {
    trees[[1]]
  } else if (length(trees) == 0) {
    NULL
  } else {
    class(trees) <- 'multiPhylo'
    trees
  }
  
}

#' @documentIn ReadTntTree Converts text representation of a tree in TNT to an object of class `phylo`
#' @author Martin R. Smith
#' @export
TNTText2Tree <- function (treeText) {
  treeText <- gsub("(\\d+)", "\\1,", treeText, perl=TRUE)
  treeText <- gsub(")(", "),(", treeText, fixed=TRUE)
  treeText <- gsub("*", ";", treeText, fixed=TRUE)
  # Return:
  read.tree(text=gsub(", )", ")", treeText, fixed=TRUE))
}

#' Extract taxa from a matrix block
#' 
#' @param matrixLines lines of a file containing a phylogenetic matrix
#'  (see ReadCharacters for expected format)
#' @template characterNumParam
#' @template sessionParam
#' 
#' @return Matrix with n rows, each named for the relevant taxon, and c columns,
#' each corresponding to the respective character specified in `character_num`
#' 
#' @keywords internal
#' @export
ExtractTaxa <- function (matrixLines, character_num=NULL, session=NULL) {
  taxonLine.pattern <- "('([^']+)'|\"([^\"+])\"|(\\S+))\\s+(.+)$"
  
  taxonLines <- regexpr(taxonLine.pattern, matrixLines, perl=TRUE) > -1
  # If a line does not start with a taxon name, join it to the preceding line
  taxonLineNumber <- which(taxonLines)
  previousTaxon <- vapply(which(!taxonLines), function (x) {
    max(taxonLineNumber[taxonLineNumber < x])
  }, integer(1))
  
  
  taxa <- sub(taxonLine.pattern, "\\2\\3\\4", matrixLines, perl=TRUE)
  taxa <- gsub(" ", "_", taxa, fixed=TRUE)
  taxa[!taxonLines] <- taxa[previousTaxon]
  uniqueTaxa <- unique(taxa)
  
  tokens <- sub(taxonLine.pattern, "\\5", matrixLines, perl=TRUE)
  tokens <- gsub("\t", "", gsub(" ", "", tokens, fixed=TRUE), fixed=TRUE)
  tokens <- vapply(uniqueTaxa, function (taxon) paste0(tokens[taxa==taxon], collapse=''), character(1))
  
  tokens.pattern <- "\\([^\\)]+\\)|\\[[^\\]]+\\]|\\{[^\\}]+\\}|."
  matches <- gregexpr(tokens.pattern, tokens, perl=TRUE)
  
  n_char <- length(matches[[1]])
  
  if (!is.null(session)) {
    shiny::updateNumericInput(session, 'character_num', max=n_char)
  }
  
  if (!exists("character_num") || is.null(character_num)) {
    character_num <- seq_len(n_char)
  } else if (any(character_num > n_char) || any(character_num < 1)) {
    return(list("Character number must be between 1 and ", n_char, "."))
    character_num[character_num < 1] <- 1
    character_num[character_num > n_char] <- n_char
  }
  
  tokens <- t(vapply(regmatches(tokens, matches),
                     function (x) x[character_num, drop=FALSE],
                     character(length(character_num))))
  if (length(character_num) == 1) {
    tokens <- t(tokens)
  } else if (length(character_num) == 0) {
    stop("No characters selected")
  }
  rownames(tokens) <- uniqueTaxa
  
  # Return:
  tokens
}

#' Read characters from Nexus file
#'
#' Parses Nexus file, reading character states and names
#'
#' Tested with nexus files downloaded from MorphoBank with the "no notes"
#' option, but should also work more generally.
#'
#' Do [report](https://github.org/ms609/TreeSearch/issues "New GitHub Issue")
#' incorrectly parsed files.
#'
#' @param filepath character string specifying location of file
#' @template characterNumParam
#' @template sessionParam
#'
#' @return A matrix whose row names correspond to tip labels, and column names
#'         correspond to character labels, with the attribute `state.labels`
#'         listing the state labels for each character; or a character string
#'         explaining why the character cannot be returned.
#'
#' @author Martin R. Smith
#' @references
#'   Maddison, D. R., Swofford, D. L. and Maddison, W. P. (1997)
#'   NEXUS: an extensible file format for systematic information.
#'   Systematic Biology, 46, 590-621.
#' @export
#'
ReadCharacters <- function (filepath, character_num=NULL, session=NULL) {
  
  lines <- readLines(filepath, warn=FALSE) # Missing EOL is quite common, so 
                                           # warning not helpful
  nexusComment.pattern <- "\\[[^\\]*\\]"
  lines <- gsub(nexusComment.pattern, "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  semicolons <- which(RightmostCharacter(lines) == ';')
  upperLines <- toupper(lines)
  
  matrixStart <- which(upperLines == 'MATRIX')
  if (length(matrixStart) == 0) {
    return(list("MATRIX block not found in Nexus file."))
  } else if (length (matrixStart) > 1) {
    return(list("Multiple MATRIX blocks found in Nexus file."))
  } else {
    matrixEnd <- semicolons[semicolons > matrixStart][1]
    if (lines[matrixEnd] == ';') matrixEnd <- matrixEnd - 1
    
    matrixLines <- lines[(matrixStart + 1):matrixEnd]
    tokens <- ExtractTaxa(matrixLines, character_num, session)
    if (is.null(character_num)) character_num <- seq_len(ncol(tokens))
    
    ## Written with MorphoBank format in mind: each label on separate line,
    ## each character introduced by integer and terminated with comma.
    labelStart <- which(upperLines == 'CHARLABELS')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      #attr(dat, 'char.labels')
      colnames(tokens) <- lines[labelStart + character_num]
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple CharLabels blocks in Nexus file."))
    }
    
    stateStart <- which(upperLines == 'STATELABELS')
    if (length(stateStart) == 1) {
      stateEnd <- semicolons[semicolons > stateStart][1]
      stateLines <- lines[stateStart:stateEnd]
      stateStarts <- grep("^\\d+", stateLines)
      stateEnds <- grep("[,;]$", stateLines)
      if (length(stateStarts) != length(stateEnds)) {
        warning("Could not parse character states.")
      } else {
        attr(tokens, 'state.labels') <-
          lapply(character_num, function (i)
            stateLines[(stateStarts[i] + 1):(stateEnds[i] - 1)]
          )
      }
    } else {
      if (length(labelStart) > 1) {
        return(list("Multiple StateLabels blocks in Nexus file."))
      }
    }
  }
  
  # Return:
  tokens
}


#' @describeIn ReadCharacters Read characters from TNT file
ReadTntCharacters <- function (filepath, character_num=NULL, session=NULL) {
  
  lines <- readLines(filepath, warn=FALSE) # Missing EOL might occur in user-
                                           # generated file, so warning not helpful
  tntComment.pattern <- "'[^']*']"
  lines <- gsub(tntComment.pattern, "", lines)
  lines <- trimws(lines)
  lines <- lines[lines != ""]
  
  semicolons <- which(RightmostCharacter(lines) == ';')
  upperLines <- toupper(lines)
  
  matrixStart <- which(upperLines == '&[NUM]')
  if (length(matrixStart) == 0) {
    return(list("&[num] entry not found in TNT file."))
  } else if (length (matrixStart) > 1) {
    return(list("Multiple &[num] entries found in TNT file."))
  } else {
    matrixEnd <- semicolons[semicolons > matrixStart][1]
    if (lines[matrixEnd] == ';') matrixEnd <- matrixEnd - 1
    
    matrixLines <- lines[(matrixStart + 1):matrixEnd]
    tokens <- ExtractTaxa(matrixLines, character_num, session)
    labelStart <- which(upperLines == 'CHARLABELS')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      #attr(dat, 'char.labels')
      colnames(tokens) <- lines[labelStart + character_num]
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple CharLabels blocks in Nexus file."))
    }
    
    
    labelStart <- which(upperLines == 'CNAMES')
    if (length(labelStart) == 1) {
      labelEnd <- semicolons[semicolons > labelStart][1]
      if (lines[labelEnd] == ';') labelEnd <- labelEnd - 1
      charLines <- lines[labelStart + character_num]
      charLine.pattern <- "^\\S+\\s\\d+\\s(\\w+)(.*)\\s*;\\s*$"
      
      # Character labels
      charNames <- gsub(charLine.pattern, "\\1", charLines, perl=TRUE)
      colnames(tokens) <- gsub("_", " ", charNames, fixed=TRUE)
      
      # State labels
      stateNames <- gsub(charLine.pattern, "\\2", charLines, perl=TRUE)
      attr(tokens, 'state.labels') <- lapply(stateNames, function (line) {
        states <- strsplit(trimws(line), "\\s+", perl=TRUE)[[1]]
        trimws(gsub("_", " ", states, fixed=TRUE))
      })
    } else {
      if (length(labelStart) > 1)
        return(list("Multiple cnames entries in TNT file."))
    }
  }
  
  # Return:
  tokens
}

#' Matrix to phydat
#' 
#' Converts a matrix of tokens to a phyDat object
#' 
#' @param tokens matrix of tokens, probably created with [ReadCharacters] 
#'               or [ReadTntCharacters].
#' @return an object of class \code{phyDat}
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
#' 
MatrixToPhyDat <- function (tokens) {
  allTokens <- unique(as.character(tokens))
  tokenNumbers <- seq_along(allTokens)
  names(tokenNumbers) <- allTokens
  matches <- gregexpr("[\\d\\-\\w]", allTokens, perl=TRUE)
  whichTokens <- regmatches(allTokens, matches)
  levels <- sort(unique(unlist(whichTokens)))
  whichTokens[allTokens == '?'] <- list(levels)
  contrast <- 1 * t(vapply(whichTokens, function (x) levels %in% x, logical(length(levels))))
  rownames(contrast) <- allTokens
  colnames(contrast) <- levels
  dat <- phangorn::phyDat(tokens, type='USER', contrast=contrast)
  
  # Return:
  dat
}

#' @describeIn ReadCharacters Read nexus characters as phyDat object
#' @author Martin R. Smith
#' @importFrom phangorn phyDat
#' @export
ReadAsPhyDat <- function (filepath) {
  MatrixToPhyDat(ReadCharacters(filepath))
}


#' @describeIn ReadCharacters Read TNT characters as phyDat object
#' @author Martin R. Smith
#' @importFrom phangorn phyDat
#' @export
ReadTntAsPhyDat <- function (filepath) {
  MatrixToPhyDat(ReadTntCharacters(filepath))
}


#' @describeIn ReadCharacters A convenient wrapper for \pkg{phangorn}'s \code{phyDat},
#' which converts a list of morphological characters into a phyDat object.
#'
#' @param dataset list of taxa and characters, in the format produced by [read.nexus.data].
#'
#' @export
PhyDat <- function (dataset) {
  nChar <- length(dataset[[1]])
  if (nChar == 1) {
    mat <- matrix(unlist(dataset), dimnames=list(names(dataset), NULL))
  } else {
    mat <- t(vapply(dataset, I, dataset[[1]]))
  }
  MatrixToPhyDat(mat)
}

#' Rightmost character of string
#'
#' @author Martin R. Smith
#' @export
#' @keywords internal
RightmostCharacter <- function (string, len=nchar(string)) {
  substr(string, len, len)
}

