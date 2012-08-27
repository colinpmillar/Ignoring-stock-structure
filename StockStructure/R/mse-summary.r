#' Methods to summaries Stock Simulations
#'
#' @param object the output from the mse function
#'
#' @return a list
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname mseSumm-methods
#'
setGeneric("mseSumm", function(object, ...) standardGeneric("mseSumm"))

#' @rdname mseSumm-methods
#' @aliases mseSumm, list
setMethod("mseSumm", c("list"), 
		function(object, scn, statSumm = TRUE, probs = c(0.2, 0.5, 0.8)) {
			lst <- as.vector("list", nrow(scn))
			
			for( i in seq_along(nrow(scn)) ) {
				x <- object[[i]]
				df0     <- as.data.frame(fbar(x))
				df0[,1] <- "fbar"
				df1     <- as.data.frame(ssb(x))
				df1[,1] <- "ssb"
				df2     <- as.data.frame(rec(x))
				df2[,1] <- "rec"
				df3     <- as.data.frame(catch(x))
				df3[,1] <- "catch"
				df4     <- as.data.frame(attr(x, "PARs"))
				out           <- rbind(df0, df1, df2, df3, df4)
				out           <- cbind(out, scn[i,][rep(1,nrow(out)),])
				out[,"runid"] <- i
				out[,"scn"]   <- paste(paste(names(scn), scn[i,], sep = "="), collapse = "; ")
				names(out)[1] <- "par"
				if(statSumm == TRUE) {
					out <- lapply(split(out, out[,c("par","year")]), 
							function(x) {
								xx <- x[1:3,]
								xx $ qtl <- probs
								xx $ data <- quantile(x $ data, probs = probs, na.rm = TRUE)
								xx
							})
					out <- do.call("rbind", out)
				}
				lst[[i]] <- out
			}
			lst <- do.call("rbind", lst)
			row.names(lst) <- NULL
			lst
		})

