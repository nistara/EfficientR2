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
+ thinking about ways to do the same thing in R in different ways,
+ making changes and ensuring we get the same results - testing.

This require thinking about
+ meta-programming to understand code
+ how to represent data in a "good" way
+ how to meaure how long computations take
+ how to compare results across different implementations.

In addition to providing technical information, this workshop primarily focuses on how to reason
about changing computations to improve them.  While is the thought process that is most interesting,
we also describe several R packages that are generally useful.

This workshop is for people who are familiar with R. However,
one doesn't need to be an expert and we believe that all R users
will get something out of this workshop, but don't need to understand everything - yet.

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
   + as for the required data, but have to identify all the files.
   + separate the functions from the script - why?
       + To allow us to change and test them individually, separately from running the script.
	   + Need ways to verify any changes don't change the results.

# Version Control
+ Want to make changes in a controlled and reversible manner
   + Use git (or other VS system)
   + Use `git commit` for files on a single branch so we can revert back to previous version.
   + Use branches to do experimental changes we can abandon or merge.


## Too Much Memory
+ Run the code in BIG.R Let it fail.
  + show the size of the object being allocated.
  + How can we understand the dependencies
  + createVariableGraph()
 
# Choosing Appropriate Data Structures 
+ allstages as a 4-D array or a list of lists of 2-D matrices.
+ Later
    + 2-D Lookup tables rather than data frame
    + 2-column matrix indexing


## Understanding code
   + Find the function definitions.
      + Functions in separate files so can source.
      + Can find them programmatically without sourcing the file - CodeDepends::getFunctionDefs()
   + Find the data input commands
       + No hard-coded file paths
       + Make data available.
   + find the for() loops, if()	statements in the script.
	   
## Compare two scripts
   + Do a diff of TOY.R and BIG.R
       + See ~/Book/CodeReview/compareScripts.xml
 
## Speed
 + So going to make TOY.R faster before we go about making BIG.R run at all.
 + Why - because we want to guarantee the same results before we change the data structures.
     + Why - so we can compare the results. If we change everything, we can't necessarily compare
       the results.
	 +  Especially the intermediate results.
 
## Measuring Run-time Code Characteristics
+ Profiling - `Rprof()` and `summaryRprof()`
+ CallCounter package
   + Determine how many times each function is called.
+ `trace()` and collecting information about calls, e.g. length, class, typeof, dim of each argument.
  + Need to know the maximum number of elements in each call.
+ Time the individual top-level expressions
   + `timeEval()` in timeEval.R
+ Test the results are the same
   + At the end
   + At each top-level expression
      + Modify script/code to create temporary assignments/write to RDS files.
	  + Do it programmatically
	      +  `storeEval()` and storeEval2.R


## Avoid
+ `ifelse()`
+ `Vectorize()`
    + If you do use this, specify the arguments being looped over and the scalars.
    + Vectorize functions using vectorized operations, not loops over arguments
+ `||` when you mean `|` (and `&&` versus `&`)


# Meta-programming
+ variable timelines
+ Variable graphs



# Packages
 + [CodeDepends](https://github.com/duncantl/CodeDepends)
 + [CodeAnalysis](https://github.com/duncantl/CodeAnalysis)
 + [rstatic](https://github.com/nick-ulle/rstatic)
 + codetools
 + R's own base package.
