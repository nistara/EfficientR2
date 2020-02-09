# ** Parsing R code
p = parse("TOY.R")

p[[1]]
p[[3]]
p[[4]]
p[[2]][[1]]
class(p[[2]][[1]]) # it's a symbolic reference to a function

# How many function calls are there?
w = sapply(p, is.call)
table(w)

w = sapply(p, function(x) is.call(x) && x[[1]] == "setwd")
table(w)

as.list(p[w])

w = sapply(p, function(x) class(x) == "=") # or <- depending on assignment type

e = p[[10]]
length(e)
class(e)
typeof(e)
e[[1]] # assignment function
e[[2]] # variable name
e[[3]] # values

w = sapply(p, function(x) class(x) == "=" && is.call(x[[3]]) && x[[3]][[1]] == "read.csv")
table(w)
as.list(p[w]) # NOTE: you need to say as.list
p[w]

f = p[[4]][[3]]
class(f) # it's a function call - not yet evaluated
f(10, 2)
f2 = eval(f)
f2(10, 2)
class(f2) # it's a "function" coz it's been evaluated
formals(f2)
body(f2)
class(body(f2)) # means it's a block of code = "{"
body(f2)[[1]]
body(f2)[[2]]


# Removing rm(list = ls()) bit before starting to evaluate
# This way you're not modifying the actual code
p[[1]]
p = p[-1] # to remove rm(list = ls()

e = new.env() # create new env
lapply(p, eval, e) # eval p's components in new env e

# NOTE: p is an "abstract syntax tree". e.g. ref: https://ruslanspivak.com/lsbasi-part7/

codetools::findGlobals(f2, FALSE) # returns global functions and variables
f2

# Do a diff on BIG.R and TOY.R to see how different they are (there are 3 differences)

# ** Find functions
## NOTE: Vectorize() is a bad function! The bane of our existence
## It takes a function and returns a function

## NOTE: Duncan's manually moved the functions to another file
## Then e = new.env
## Source functions and tell me how many functions? there are in it?

w = sapply(p, function(x) class(x) %in% c("=", "<-") && is.name(x[[2]]) && is.call(x[[3]]) && x[[3]][[1]] == "function")
table(w)

# names of functions in the script
sapply(p[w], function(x) as.character(x[[2]]))

# find the Vectorize functions
w = sapply(p, function(x) class(x) %in% c("=", "<-") && is.name(x[[2]]) && is.call(x[[3]]) && x[[3]][[1]] == "Vectorize")
sapply(p[w], function(x) as.character(x[[2]]))

# how many for loops?
table(sapply(p, class)) # 7 for loops
# 5 of the calls are setwd, the rest are useless - not doing anything

w = sapply(p, class)
which(w == "for")
p[[6]]
length(p[[6]])
p[[6]][[1]]
p[[6]][[2]]
p[[6]][[3]]
p[[6]][[4]]



# ** NOTE
# When you're designing a function, you need to build tests before you run the code,
# and confirm that your functions are working as expected
# ** Run the code and get the results (which form your comparizon base)
# Compare:
# BEST - the last result (not really the best)
# All top level variables
# All intermediate values

# You need to build tests for each individual function
# Want to change a function and test it immediately

# Collect all calls, their arguments and their results, plus we need the
# global variables that were in effect at that time

# Can add a line at the beginning of a function
# args = list(x = x, base = base)
# then do the computation
# finally, get:
# functionNameInfoGet = list(args, answer, globals)

# There's a function in R that can be used to get the above: trace()
# trace("functionName")
# tells you if a function is called and how often

# trace("functionName", quote(ctr <<- ctr + 1))

ctr = 0
quote(ctr <<- ctr + 1)
ctr = 0
e = quote(ctr <<- ctr + 1)
eval(e)
trace(f2, quote(ctr <<- ctr + 1))

f2(10, 2)
ctr
replicate(3, f2(10, 4))
ctr
# can tell trace where to run within a function. So can tell it to do it in the end.
# then you get the result in the end
# Duncan created this package call "CallCounter"
# Collects no. of calls, arguments, and what it returns (you tell it what to collect)
# Also gives dimensions of objects (to see if it's vectorized)

# So once we collect the above arguments using testFunCalls()


# ** Vectorize eg
g = function(x) x + 1
f = Vectorize(g)
f(1:10)

# Look at f
f
# Vectorize broke our tracing mechanism in terms of function calls, and it
# also broke our code!!
# it's doing all scalar calls so it takes longer
# It's needed sometimes for functions like curve() coz there are some things
# you can't vectorize

trace(g)
f(1:10)

ls(environment(f))
environment(f)$FUN # that's g!!
b = environment(f)$FUN
environment(f)$FUN = b
trace(b) # not working right now

# ** Profiling code
Rprof("tmp.prof") # turn on profiler
replicate(1000, median(rnorm(1000)))
Rprof(NULL) # turn off profiler

p = summaryRprof("tmp.prof")
p$by.self # Duncan lives by this

# ** Post profiling
# Need to fix ifelse() and Vectorize()
# Only one functions calls it: finalinflowprep
# There's a nested ifelse() in the function above
# It took 15% of the time even though it was called only once

# Old function
SpringDeltaVcprep=function(RcstarWinter,month, RC,p){
  DeltaVc=ifelse((springcoeff[1]+springcoeff[2]*RC+springcoeff[3]*as.numeric(TaLookup(month,p))+springcoeff[4]*as.numeric(QLookup(month,p))) <0,0,
                 (springcoeff[1]+springcoeff[2]*RC+springcoeff[3]*as.numeric(TaLookup(month,p))+springcoeff[4]*as.numeric(QLookup(month,p))))
  AdjustedforRcWinter=mround(DeltaVc#-1.5*10^6
                             ,bin)
  return(AdjustedforRcWinter)
}


# Newer function
SpringDeltaVcprep=function(RcstarWinter,month, RC,p){
    tmp = (springcoeff[1]+springcoeff[2]*RC+springcoeff[3]*as.numeric(TaLookup(month,p))+springcoeff[4]*as.numeric(QLookup(month,p)))
  DeltaVc=ifelse( tmp < 0, 0, tmp)
  mround(DeltaVc, bin)
}


SpringDeltaVcprep=function(RcstarWinter,month, RC,p){
    tmp = (springcoeff[1]+springcoeff[2]*RC+springcoeff[3]*as.numeric(TaLookup(month,p))+springcoeff[4]*as.numeric(QLookup(month,p)))
    DeltaVc = pmax(tmp, 0)
    mround(DeltaVc, bin)
}


SpringDeltaVcprep=function(RcstarWinter,month, RC,p){
    tmp = (springcoeff[1]+springcoeff[2]*RC+springcoeff[3]*as.numeric(TaLookup(month,p))+springcoeff[4]*as.numeric(QLookup(month,p)))
    DeltaVc = tmp[ tmp < 0] = 0# now we'll have to check for NAs
    mround(DeltaVc, bin)
}

# ** Another function eg

finalinflowprep=function(month,Q){ #converts from cfs to taf
  monthlyQ=ifelse(month == "January"|| month=="March"|| month=="June"|| month=="July"||month=="August"||month=="October"||month=="December", 
                  Q*1.98*31, #Q*cfs to af* day number in month
                  ifelse(month=="February", 
                         Q*1.98*28,
                         #as.numeric(Lookupy[which(Lookupy[,1]==month),3])*1.98*28,
                         Q*1.98*30))
  return(monthlyQ)
}


MonthDays = c(January = 31)

finalinflowprep=function(month,Q)
    MonthDays[month] * Q*1.98

body(finalinflowprep)
body(finalinflowprep)[[2]][[2]]
body(finalinflowprep)[[2]][[2]][[2]]
class(body(finalinflowprep)[[2]][[2]][[2]])

body(finalinflowprep)[[2]][[2]][[2]] = MonthDays # It's inserted in the code directly then. Inlined the constant expression
# Constant propagation


# ** Another function

QLookupprep=function(month,p){ #this could be VC, Vc or Vcstates
  Qprep=Lookupy[which(Lookupy[,1]==month & Lookupy[,2]==p),3]
  as.numeric(Qprep)
}

# Lookupy is a dataframe with 4 columns
# `[` is expensive - doing a lot fo stuff with rownames, fixing the dataframe, etc

# Should use `[[` instead of [,1]
# It's pulling the vector and working with it, vs working with the dataframe
QLookupprep=function(month,p){ #this could be VC, Vc or Vcstates
  Qprep=Lookupy[[3]][ which(Lookupy[[1]]==month & Lookupy[[2]]==p)]
  as.numeric(Qprep)
}

# Now vectorize
QLookupprep=function(month,p){ #this could be VC, Vc or Vcstates
  Qprep=Lookupy[[3]][ Lookupy[[1]]==month & Lookupy[[2]]==p]
  as.numeric(Qprep)
}

# Then use matrix magic

