# Making R Code Efficient 
### Friday February 7, 9am - 12.
### Duncan Temple Lang


This is a hands-on case study in taking an existing R project/script and making it more efficient.
This means 
+ making it use less memory so it is actually feasible to run, and 
+ making it faster to run so that we can go to "larger" runs.

We start from an existing script and make it approximately 300 times faster
and uses significantly less memory.

This will involve 
+ becoming familiar with code we didn't write,
+ how to figure out what we need to change and in what order,
+ thinking about ways to do the same computations in R in different ways,
+ making changes and ensuring we get the same results – testing.

We think about
+ meta-programming to understand code
+ how to represent data in a "good" way
+ how to measure how long computations take
+ how to compare results across different implementations.

In addition to providing technical information, this workshop primarily focuses on how to reason
about changing computations to improve them.  While it is the thought process that is most interesting,
we also describe several R packages and functions that are generally useful.

This workshop is for people who are familiar with R. However,
one doesn't need to be an expert and we believe that all R users
will get something out of this workshop, but don't need to understand everything.

Also, we'll answer lots of questions about this case study and more general questions,
so it should be of value to anyone.


# Efficiency in R
 + Speed/run-time 
 + Memory usage

# Thought process
+ This workshop is all about the thought process
   + what to do next, 
   + how to decide where to spend time,
   + when to give up on an approach and shift to something else.
   
# Case Study
+ Given an R script
   + ask for the required data, but have to identify all the files.
   + separate the functions from the script - why?
       + To allow us to change and test them individually, separately from running the script.
	   + Need ways to verify any changes don't change the results.

# Version Control 
+ Want to make changes in a controlled and reversible manner
   + Use git (or other VC system)
   + Use `git commit` for files on a single branch so we can revert back to previous version, OR
   + Use branches to do experimental changes we can abandon or merge.

# Too Much Memory
+ Run the code in BIG.R Let it fail.
  + show the size of the object being allocated.
  + How can we understand the dependencies
  + createVariableGraph()
 
# Choosing Appropriate Data Structures 
+ allstages as a 4-D array or a list of lists of 2-D matrices.
   + if desired, write intermediate results to disk - `saveRDS()`
+ Later
    + 2-D Lookup tables rather than data frame
    + 2-column matrix indexing

# Understanding code
   + Find the function definitions.
      + Functions in separate files so can source.
      + Can find them programmatically without sourcing the file - CodeDepends::getFunctionDefs()
   + Find the data input commands
       + No hard-coded file paths
       + Make data available.
   + find the for() loops, if()	statements in the script.
   + don't have `rm(list = ls(all = TRUE))` or source() in separate environment 
   ```
   e = new.env()
   source("rm.R", e)
   ```
	   
# Compare two scripts
   + Do a diff of TOY.R and BIG.R
       + See ~/Book/CodeReview/compareScripts.xml
 
# Speed
 + So going to make TOY.R faster before we go about making BIG.R run at all.
 + Why - because we want to guarantee the same results before we change the data structures.
     + Why - so we can compare the results. If we change everything, we can't necessarily compare
       the results.
	 +  Especially the intermediate results.
 
# Measuring Run-time Code Characteristics
+ Profiling - [`Rprof()`](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/Rprof) and [`summaryRprof()`](https://www.rdocumentation.org/packages/utils/versions/3.6.2/topics/summaryRprof)
+ [CallCounter](https://github.com/duncantl/CallCounter) package
   + Determine how many times each function is called.
+ [`trace()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/trace) and collecting information about calls, e.g. length, class, typeof, dim of each argument.
  + Need to know the maximum number of elements in each call.
+ Time the individual top-level expressions
   + `timeEval()` in timeEval.R
+ Test the results are the same
   + At the end
   + At each top-level expression
      + Modify script/code to create temporary assignments/write to RDS files.
	  + Do it programmatically
	      +  `storeEval()` and storeEval2.R
+ cbind() variations
  + matrix
  + more specific context
  + C code.
+ inlining global constants
  + determining they are constants
  + modifying functions programmatically
  
# Avoid
+ [`ifelse()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/ifelse)
+ [`Vectorize()`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/Vectorize)
    + If you do use this, specify the arguments being looped over and the scalars.
    + Vectorize functions using vectorized operations, not loops over arguments
+ `||` when you mean `|` (and `&&` versus `&`)


# Meta-programming
+ variable timelines
+ Variable graphs


# Improvements
+ replace ifelse() with | and %in% and subsetting. 
+ replace Vectorize() with vectorized computations.
+ [,1] versus [[1]] in a data.frame. - Lookupy
+ accumulate - loop over fewer elements in a single loop rather than all elements in 2 nested
  loops. Matrix indices from `which(, arr.ind = TRUE)`
+ improve cbind()
+ subset inputs to call, not the results of the call.
+ Create matrices outside of loop when same size.  `x[] = 0` versus `x = matrix(...)` each time.
+ removing variables we create but don't use.
+ `benefit(x, y,z)[w]` versus `benefit(x[w], y[w], z)`
+ `MaxQw=max(Lookupy[Lookupy[,1]=="April" | Lookupy[,1]=="May" | Lookupy[,1]=="June" |
  Lookupy[,1]=="July" | Lookupy[,1]=="August" | Lookupy[,1]=="September" |
  Lookupy[,1]=="October",3])`
   versus
  `MaxQw=max(LookupyQ[c("April" , "May" , "June" , "July", "August", "September", "October"),])`

# Packages
 + [CodeDepends](https://github.com/duncantl/CodeDepends)
 + [CodeAnalysis](https://github.com/duncantl/CodeAnalysis)
 + [rstatic](https://github.com/nick-ulle/rstatic)
 + codetools
 + R's own base package.
