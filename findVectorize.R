e = parse("functions_orig.R")
sapply(e, function(x) is.call(x) && class(x) %in% c("<-", "=") && is.call(x[[3]]) && is.symbol(x[[3]][[1]]) && x[[3]][[1]] == "Vectorize")
sapply(e[w], function(x) as.character(x[[2]]))

# OR

library(rstatic)
v = find_nodes(to_ast(e), function(x) is(x, "Call") && is(x$fn, "Symbol") && x$fn$value == "Vectorize")
sapply(v, function(x) x$parent$write$value)
