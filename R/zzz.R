# #c".sm.Options" <- list(display = "none")
#
# op <- options()
# op_sm <- list(display = "none")
#
# .onLoad <- function(library, pkg)
# {
#   #unlockBinding(".sm.Options", asNamespace("sm"))
#   toset <- !(names(op_sm) %in% names(op))
#   if(any(toset)) options(op_sm[toset])
#
#
#   invisible()
# }
