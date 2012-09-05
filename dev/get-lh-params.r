





if (0) {
  

library(XML)
#library(rfishbase)


pages <- function(fish, type = "LW") {
 
  gtype <- switch(type,
                  "Linf" = "graphLengthFM01DL",
                  "LW"   = "LW_graph01DL")
  
  paste0("http://fishbase.org/graph/", gtype, ".cfm?",
         "ID=",fish $ ID,
         "&fc=",fish $ fc,
         "&genusname=", fish $ genusname,
         "&speciesname=", fish $ speciesname)
}

getTables <- function(fish, type) {
  page <- pages(fish, type)
  tabs <- try(readHTMLTable(page), TRUE)
  i <- 1
  while((length(tabs) == 0 | class(tabs) == "try-error") & i < 10) {
    i <- i + 1
    cat("try", i, "\n"); flush.console()
    tabs <- try(readHTMLTable(page), TRUE)  
  }
  if (length(tabs) == 0 | class(tabs) == "try-error") stop("unable to access web page")
  
  n <- c( as.character(tabs[[1]][[1]][1]), 
           as.character(tabs[[1]][[2]][1]),
           as.character(tabs[[1]][[1]][2]) )
  n <- as.numeric(gsub("[a-z|A-Z|(\r)|(\n)|(\t)|[|]|=| |]", "", n))
  
  sp <- tabs[1:n[1] + 2]
  others <- tabs[1:n[2] + n[1] + 4]
  misc <- tabs[1:n[3] + n[1] + n[2] + 6]
  names(sp) <- names(others) <- names(misc) <- NULL
  
  sp <- do.call(rbind, sp)
  others <- do.call(rbind, others)
  misc <- do.call(rbind, misc)
  
  out <- rbind(sp, others, misc)
  out $ which <- rep(c("sp", "other", "misc"), n)
  
  xynames <- strsplit(gsub("graph| ", "", as.character(tabs[[1]][3,1])), "vs")[[1]]
  names(out) <- c("species", xynames, "which")
  
  out
}


dat <- data.frame(ID = c(69, 1381)
    fc = 183,
    genusname = c("Gadus", "Melanogrammus)",
    speciesname = c("morhua", "aeglefinus")



codLW   <- getTables(dat[1,], type = "LW")
codLinf <- getTables(dat[1,], type = "Linf")

hadLW   <- getTables(dat[2,], type = "LW")
hadLinf <- getTables(dat[2,], type = "Linf")

}



