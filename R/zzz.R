## family$initialize and the EM machinery eval() code that assigns
## these in the calling frame; R CMD check cannot see through eval().
utils::globalVariables(c("n", "mustart"))
