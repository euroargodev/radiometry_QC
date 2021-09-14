my_menu <- function(choices, title=NULL) { 
  # equivalent to menu() but compatible with non interactive sessions
  
  choice_lines = paste0(as.character(1:length(choices)), rep(": "), choices)
  cat(title, "", choice_lines, sep="\n")
  
  repeat {
    cat("\nSelection: ")
    
    if (interactive()) {
      input <- as.numeric(readLines(stdin(), 1))
    } else {
      input <- as.numeric(readLines(file("stdin"), 1))
    }

    if ( !is.na(input) ){
      if ( input>=0 & input<=length(choices) ) {
        closeAllConnections()
        return( input )
      }
    }
    
    cat("Choose one of the menu items, or 0 to exit")
  }
}


scientific_formatter <- function(x) { 
  # format data into scientific notation with 10^ (instead of e)
  
  x = format(x, scientific=TRUE)

  x = sub("\\+0", "+", x)
  x = sub("\\-0", "-", x)

  return(TeX(paste0(sub("e", " 10^{", x),"}")))
}
