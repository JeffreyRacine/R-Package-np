#### single line console management
newLineConsole <- function(parentConsole = NULL){
  newConsole <- list(lineList = list(),
                     numLineElements = 0,
                     tabLen = 8,
                     lineLen = 0,
                     startCol = 1,
                     lineMsg = '')

  if(!is.null(parentConsole))
    newConsole$startCol <- parentConsole$startCol + parentConsole$lineLen
  return(newConsole)
}


charRep <- function(x, times, ...){
  if(times != 0)
    return(rep(x = x, times = times, ...))
  else
    return('')
}

toMsg <- function(msg, console = stop("no console provided")){
  ## attempt to provide minimal support for other data types
  msg <- paste(msg, collapse=' ')

  ## explode the message
  tMsg <- unlist(strsplit(msg, '\\\t'))

  ## add trailing tab if necessary
  if(substr(msg, nchar(msg), nchar(msg)) == '\t')
    tMsg[length(tMsg)+1] = ""

  if(length(tMsg) > 1){
    tmpPos <- console$lineLen + console$startCol
  
    ## do tabbing
    msg <- paste(tMsg, sapply(tMsg, function(x) {
      tmpPos <<- tmpPos + nchar(x)
      tmpPos <<- tmpPos +
        nchar(tFill <- paste(charRep(' ', console$tabLen-(tmpPos %% console$tabLen)),
                             collapse = ''))
      tFill
    }), sep='', collapse='')
  }
  list(msg = msg,
       len = nchar(msg))
}


printPush <- function(msg, console = stop("no console provided")){
  if(options('np.messages')$np.messages){
    console$numLineElements <- console$numLineElements + 1
    console$lineList[console$numLineElements] <- NA
    console$lineList[[console$numLineElements]] <- toMsg(msg = msg, console = console)
    console$lineLen <- console$lineLen + console$lineList[[console$numLineElements]]$len
    console$lineMsg <- paste(console$lineMsg, console$lineList[[console$numLineElements]]$msg, sep = '')
    cat(console$lineList[[console$numLineElements]]$msg)
    flush.console()
  }
  return(console)
}

printPop <- function(console = stop("no console provided")){
  if(console$numLineElements > 0 & options('np.messages')$np.messages) {
    cat(paste(rep('\b',console$lineList[[console$numLineElements]]$len), collapse=''))
    flush.console()
    console$lineLen <- console$lineLen - console$lineList[[console$numLineElements]]$len
    console$lineMsg <- substr(console$lineMsg, 1, console$lineLen)
    stopifnot(console$lineLen >= 0)
    console$lineList[console$numLineElements] <- NULL
    console$numLineElements <- console$numLineElements - 1
  }
  return(console)
}

printClear <- function(console = stop("no console provided")){
  if(console$numLineElements > 0 & options('np.messages')$np.messages){
    cat(paste(rep('\b', nchar(console$lineMsg)), collapse=''))
    cat(paste(rep(' ', nchar(console$lineMsg)), collapse=''))
    cat(paste(rep('\b', nchar(console$lineMsg)), collapse=''))
    flush.console()
    console$lineLen <- 0
    console <- newLineConsole(console)
  }
  return(console)
}

