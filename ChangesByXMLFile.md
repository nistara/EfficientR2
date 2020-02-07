# Sequence of Changes

This (attempts to) details the changes to TOY2.R and functions.R at each step corresponding to our XML documents.

+ Topics and timings are not necessarily in the specified or appropriate order.
+ Some topics are not actually implemented.

1. pkgsTutorials.xml
   *  just information.
1. ExploringCode.xml
   1. create TOY2.R
   1. move function definitions in TOY2.R to functions.R, including Vectorize() calls.
   1. comment out/remove the setwd() and rm() calls.
   1. change require(zoo) and na.locf to zoo:na.locf
   1. add the 3 CSV data files
1. VariableDepends.xml
   + no changes.
1. Tests.xml
   + collect the arguments for all calls.
     + Note that these will be wrong for finalinflowprep as there is a bug in the original.
1. GlobalVariables.xml
   + comment out fs and fstarvalue in 2 loops, and define locally in accumulate and firststageaccumulate
1. VariableLifetimes.xml
   + no changes
1. profiling.xml
   + no changes
1. repeatedComputations.xml
   + remove all of the ans = ; return(ans) in functions.R
   + SpringDeltaVcprep and ifelse(a < 0, 0, a) computes a twice.
   + OutgoingVcprep, OutgoingVwprep and repeating computation of lakeseasonbin(month)
   + fallsolve repeating VC + deltaVc, summersolvprep VW + deltaVw
   + Notice in springsolveprep minor deviation from idiom which is probably a mistake in the code -
     `VW+deltaVw` and then `Vw+deltaVw`, i.e., VW and Vw in same ifelse() call   
   + fallsolveprep  VC + deltaVc in 5 places
   + choosesolveprep - seasonbin(month) repeated in numerous places.
   + Didn't update Jelly, Balls, ClrCk, RBDD.
1. Vectorize.xml
   + focus on finalinflowprep as a) first, b) called directly from script and not indirectly via
     other functions which may be scalar, c) already almost vectorized.
   + finalinflowprep - change || to |, leave ifelse() call.
   + finalinflow = finalinflowprep, no Vectorize() version.
   + Didn't do the same for fishtemp, but could have.
1. IfElse_Or_Lookup.xml
   + Use named-indexing into vector of days in months.
   + Fix error for June and May being confused with which has 31 days in original code.
1. MonthAsFactor.xml
   + change TOY2.R and climate = read.csv( stringsAsFactors = FALSE)
   + Don't actually change functions. (Could drop the as.character(month) later.)
1. Vectorize2.xml
   + Find functions that are Vectorized() but which are always called with all scalar values.
   + Remove 10 Vectorize() calls for TaLookup, QLookup, WinterDeltaVc, SpringDeltaVc, ReleaseTemp,
    benefit, mixedsolve, springsolve, summersolve, fallsolve as currently only called with all
    scalar values so can call the prep versions directly.
    ```
    TaLookup = TaLookuprep
    ```
   + Leave Vectorize() definitions for finalinflow, ColdDelta, OutgoingVc, OutgoingVw, choosesolve
   + Do not replace calls to ifelse() just vectorize the three arguments in each call.
   + Ultimately want these functions to be called with vectors, but take the easy win for now.
1. IfElse2.xml   (down to 293 seconds)
   + remove all ifelse() calls in all (21) functions. Replace with 
     switch() or ifelse() or 
   ```
   ans = rep(defaultValue, len); ans[cond] = val1; ans[cond2] = val2
   ```
    but this last one is more for vectorization which will come later. But we use it here also.
     So is this vectorized and scalar.
   + pmax() in √ springsolve, √ mixedsolve for AvailableVw. Note correction to springsolve VW + deltaVw
     rather than Vw + deltaVw
	   + springsolve not giving same answer, as we expected because the bug fix change.
       + NO: ¿or switch(). (Look at final functions to see where we did that. OutgoingVcprep, choosesolve. And
               OutgoingVc but at this point that is still Vectorize()
   + √Jelly, Balls, ClrCk, RBDD - `ans = comp; ans[ans < 30] = 0; ans` - so also vectorized.
   + √monthcounter - also vectorized.
   + √fishtemp - lookup table, create from 4 variables in TOY2.R. Move them to functions. Also vectorized
   + √lakeseasonbin, √seasonbin - lookup named vector (¿Do we inline the lookup vectors?),
   + √mixedsolve (use VC for length even though should    be RC? or R? as we discover later),
   + √summersolve, √ springsolve
       + summersolve: 331 of 730K calls not the same, but all results were NAs, just numeric not logical
       + springsolve: 8976 of 224K calls not the same, fixed logic inversion for one case of whether to call benefit.
   + √ finalinflow - Wasn't finalinflow already done in Vectorize? Yes. Gives different results since
       fixed May/June reversal of number of days.
   + √ choosesolveprep - use switch? But then not vectorized in month, but focus is on removing `ifelse()`.
	    + 9347 of 182K calls were not the same.
   + √ OutgoingVcprep - switch(), OutgoingVwprep - if-else-if-else and %in%.
   + √ ColdDeltaprep - simple if - else but not vectorized in month.
   + √ SpringDeltaVc - pmax() replaces ifelse().
   + √ ReleaseTemp - replace ifelse(  ifelse, ifelse) with if-else if-else in 3 places. !!! Come back
     and vectorize and subset and introduce error below.
   + √ benefit  
   + √ firststageaccumulate, √ accumulate
       + How far do we go with accumulate and firststageaccumulate? Just ifelse() replacements with if-else.
   + [Didn't do - related to vectorization later] Recreate error	 
      + in ReleaseTemp. `w = w & Tw != 0` instead of `!w & TW != 0`.
      + accumulate and `>` rather than `>=` (See Test2.xml#794)
1. SubsetDataFrame.xml  (gets about 70 seconds back - 30 with first 3 changes, additional 40 with the inlining)
   + Separate commits for each of these.
      + Change Lookupy[,col] to Lookupy[[col]] in QLookup 
      + remove which() 
	  + remove as.numeric() - not necessary since the lookup is on a numeric
      + Inline the Lookupy columns object into these functions. Do in TOY2.R
	  + Do the same for TaLookup. (1 commit for all changes to TaLookup)
	     +  (*No - later in Inlining.xml*) ¿inline other lookup vectors - fishtemp, ...?
1. Vectorize3.xml
   + Create vectorized functions of √ OutgoingVc, √ OutgoingVw, choosesolve
   + choosesolve becomes ¿
   ```
   choosesolve=
    function(month, VW, VC, RC, RW, R, V, RcstarWinter, K, DP,p)
        mapply(choosesolveprep, month, VW, VC, RC, RW, R, V, RcstarWinter, K, DP, p)
   ```
      ¿ But apparently we if(FALSE) it out and left the Vectorize() version?
   + √ QLookup - LookupyQ[cbind()]
   + √ TaLookup
   + Inline LookupyQ and LookupyTa
   + √ Vectorize ColdDelta with rep() for the winter case and vectorized version of SpringDeltaVc which
     already vectorized math operations.
   + √ summersolve - we already vectorized, but can do better later when we call benefit, as with
     other solve functions.
   + √ benefit - vectorized.
   + √ ReleaseTemp vectorized using a 2-way lookup table and cbind. 
   + Note that since choosesolve() is still not really vectorized but uses mapply(), all the calls
     to other functions from choosesolveprep are still called with scalar values.  We have been
     creating vectorized functions, but they are not being called with vectors.
   + <font color="red">We may want to not change the ColdDeltaprep here and leave the ColdDelta as a
   Vectorize()-created function to get the timings to match the text better. We can leave this
   definition but just if(FALSE) it or whatever without hopefully generating a lot of conflicts in
   the rebaseing of the following commits.
     </font>
1. TimeEval.xml
   + <font color="orange">When implementing these, look at parallelLoops.xml and try to make the code
     match the code shown there</font>
   + √ speed up one loop in TOY2.R - contains Rcaccumstar and Rwaccumstar 
      + chain definitions of `Rcaccumstar = Rwaccumstar = matrix(NA, ..)`
	  + move outside of loop, reinitialize with `Rcaccumstar[] = NA`
      + set Rcaccumstar to NA by default and only 
	  + compute the 
	  ```
	   if(!is.na(isaccum[v,r])) {
          w = Vcstates == Vcoutdirect[v,r] & Vwstates == Vwoutdirect[v,r]
          Rcaccumstar[v,r] = Rcstar[w , S+1]
          Rwaccumstar[v,r]= Rwstar[w,S+1]
       }
	  ```
	  and no need to set the values if isaccum[v,r] is NA since done this in the reset of
	  Rcaccumstar and Rwaccumstar.
	  + √ Replace nested loop with a single loop over the indices of the cells of isaccum which are
	       not NA identified with which(, TRUE).  Small % of the cells are actually not NA so non-trivial savings.
	       + Didn't vectorize `Vcstates   ==  ... & Vwstates == ` so do in loop, but significantly
             fewer iterations.
   + √ preallocate `rangeVc = numeric(pn)` before loop calling OutgoingVcprep. And for rangeVw. rangex
   + √ Row indices for setting Best[S,3] and Best[S,4]	identical so compute only once. Same in
        final 1:pn loop.
   + √ Change Best to a numeric matrix by assigning the 6th column as rownnames.
       + √ remove as.numeric() in 8 calls to `as.numeric(Best[,])`
1. Vectorize4.xml
    + <font color="red"><b>Tests2.xml may come before this - see first sentence in this
      document.</b></font>
	   +  *We have a quite different profiling table in the code repos than at the start of Vectorize4.xml*
    + firststageaccumulate
 	    + √ move call to choosesolve() out of loop with vector arguments. Leave loop over elements.
		+ √ Uncovers bug in ReleaseTemp, but actually in mixedsolve and we need to use the length of
          RC not VC to define the default value.
		+ (already introduced this in Vectorize3?) Change Vectorize() version of ColdDelta to
	    ```
		ColdDelta = ColdDeltaprep =
          function(month, RcstarWinter, RC, p)
          { 
            if(lakeseasonbin(month)=="winter")
               rep(WinterDeltaVc(month, p), length(RC))
            else
              SpringDeltaVc(RcstarWinter, month, RC, p) 
          }
		```
		+ change the loop in firststageaccumulate. Should have left this to Vectorize5 apparently,
          but seems like it is mentioned in the text as if we did it.
		   + operate on the subset of elements from choosesolve() result that are positive non-NA
		   + update these in a loop, then with vectorized assignment but have to use mapply() to get
             the indices.
1. Inlining.xml
   + **Note:** The benefit of inlining decreases with vectorization changes since we call each function fewer
     times (but with vectors rather than scalars) and so the number of times we lookup these global variables
	 is less and each time is constant, so the cumulative time goes down. It is still useful, but
	 not if the computations to inline the objects takes more time than the cumulative lookup time.
   + √ Do these conditionally in TOY2.R with `if(UseInline)` and define this at the top.
   + √ Inline LookupyQ in QLookup. ¿Didn't we already do this?
   + Directly in functions.R without waiting for variables to be created in TOY2.R first.
       + √ seasonbin - SeasonbinLookup and LakeseasonbinLookup
       + √ monthcounter - month.name
   + 2 ways to inline into SpringDeltaVc & WinterDeltaVc	   
       + Inline each element of springcoeff (and wintercoeff) into the expression within the body.
	      + SpringDeltaVc - springcoeff, bin. Drop term with springcoeff[2] since it is 0.
	      + WinterDeltaVc - not worth it as no gain in time.
	   + √ Add a parameter named springcoeff whose default value is the literal vector
        	 after these variables are defined in TOY2.R. It **seems** to make a little difference (.6
        	 of a second).
 		   + Do this for other functions to via code analysis of the global variables and which
              variables change in TOY2.R.
   + √ <font color="red">Add to text</font> Collect all of the computations that do the inlining in TOY2.R into a separate
     script and conditionally source that. This makes it all more modular and reusable w/o
     source'ing TOY2.R. But there are dependency issues.
1. cbind.xml
   + Inline `as.character(max(Vw))` and `as.character(max(Vc))` - conditionally on if(UseInline)
      + perhaps make Vw and maxVw parameter and inline those values when we have computed them
	    and rewrite the code to use maxVw.
   + [didnd't bother implementing in this repos. See LaurenOrig.] switch calls to cbind() to use matrix(c(a, b), ...) . We will revert this.
1. pmax.xml  (22 - 25 seconds)
   + √ remove pmax() calls in SpringDeltaVc, springsolve, mixedsolve, summersolve, fallsolve, to use, e.g.,
    ```
 	AvailableVc = VC + deltaVC
    AvailableVc[AvailableVc < 0] = 0
	```
1. Vectorize5.xml ( 3 seconds)
   + At this point, there are still a large number of scalar calls, e.g. 421,875 to choosesolve
     (from accumulate),   168750 to mixedsolve, but some vectorized calls to each of these. 
   + [Done in Vectorize 4] firststageaccumulate
   		+ loop over only the values for which we will change `fs`, i.e. `!is.na(tmp) & tmp > 0` values where tmp is the value from
          choosesolve. See the final version in functions.R (LaurenOrig). 
		+ Vectorize calls OutgoingVc and OutgoingVw before the loop.
		+ See Vectorize5.xml text for details.
   + √ accumulate
      + NOTE: Appears from text that we hadn't changed this earlier to vectorize choosesolve and loop over
        only the values that need updating in result.
      + Similar approach as firststageaccumulate.
	  + Use `which(, arr.ind = TRUE)` to get rows and columns to update and replace nested loop with
        single loop.
1. parallelLoops.xml - focus on TOY2.R script, not functions.
    + <font color="blue">Somewhat mis-named as we actually change some of the code in the loop,
      but not to be parallel</font>
	+ Primary loop taking time  
	   + [Seem to have done this specific step in Vectorize3.xml step above] Change a loop like accumulate in Vectorize5.xml above to only loop over the values that will
         be changed, not all values.
       + √ compute `apply(onestage, 1, max, na.rm = TRUE)` just once in call to `fstar[,S] = ...`
       + √ whichxstar[,S] computed once, not three times.
       + √ remove matrix() in call matrix(accumulate(....))
	   + √ Move to outside loop
     	 ```
	     Vwoutdirect=matrix(OutgoingVw(S,VwSpace,Rwspace,p),
                   nrow=length(Vstates),ncol=length(Rdecs))  
 	     ```
	     and update within loop via `Vwoutdirect[] = OutgoingVw(S,VwSpace,Rwspace,p)`
   	     and same for Vcoutdirect
       + [done in earlier step] Rcaccumstar and Rwaccumstar - similarly outside loop.
       + √ Then move these all the way outside of the outer loops.
	   + √ Same for currentB, accumB, Vcoutacc, Vwoutacc
	   + √ Replace `finalx = matrix(choosesolve(..)); allstages[,,i, S] = finalx` with `allstages[,,i,S]
            = choosesolve()`
	      dropping the matrix() and the intermediate variable. (In BiggerGrid.xml, we have to add back
          the matrix())
       + √ chain definition of fstar, etc. with one call to stagepolicy()
   	      ```
     	  fstar = whichxstar = xstar = Rcstar = Rwstar = stagepolicy(Vstates,NoofStages)
      	  ```
       + √ Fix repeated calls in `ifelse(apply(LastStage))...` with subsetting assignment.
       + √ Replace numerous ifelse() computations.
       +  Note about onestage summation and using apply() which is apparently slower in this case.	   
    + Second longest-running loop.
	   + √ Merge two loops for rangeVc and rangeVw into one that sets each element of each vector in
         the body. See 599 of parallelLoops.xml.
	   + √ allocate rangex, rangeVc and rangeVw in chain outside of loop.
	   + ‡ These may have slightly slowed this loop down!		 
    + 3rd loop - for(i in 1:pn) - Rcaccumstarfirst
	   + √ re-implement ifelse() statements
	   + √ check if any elements in isaccumfirst are not NA. Turns out no! So avoid that nested loop.
 	       + Get rid of the inner loop since from 1 to 1.
 	   + ¿ `rangex[i] = ...` ???? Read more of section in parallelLoops.xml
    + √ Change `VcSpace=apply(basics,2,function(x)Vcstates)` to  `VcSpace = matrix(Vcstates, length(Vstates), ncol(basics))`
1. Test2.xml - <font color="red">This maybe comes after TestsCheckpoints.xml</font>
    + Fix error in accumulate replacing `>` with `>=`, i.e., `w = !is.na(tmp) & tmp >= 0`
1. TestsCheckpoints.xml
    + Fix error in ReleaseTemp() - `!w & Tw != 0` See IfElse2.xml
1. cbind3.xml
    + bind2.c - C code for cbind2.
	+ cbind2 R function
	+ inline the external ptr to the address of R_cbind2 in the cbind2() R function
1. benefitCallSubset.xml
    + In the solve functions, change the call to `ans[w] = benefit(x, y, z)[w]` to
	   `ans[w] = benefit(x[w], y[w], z[w])`
	    + Eventhough subsetting numerous vectors, the number of true values in `w` is
		  so small that we save a lot of computation in benefit() and its sub-functions
		+ Get a speed of almost 2.
1. cbind4.xml
    + Didn't implement this in the timings, so no need for changes to functions.R or TOY2.R
1. cbind5.xml
    + √ Recognize that QLookup and TaLookup are called only with scalar values, so cbind2() is
      unnecessary and can use regular subsetting.
	+ √ Similarly, in 2 calls in ReleaseTemp, we have a vector and a scalar for subsetting and can
      use    `table[vec, scalar]` and don't have to replicate the scalar and use cbind2().  
    + Comment out the unused variables in TOY2.R
1. BiggerGrid.xml
   + Data structure changes in TOY2.R to allstages to create as list of lists of matrices and then
     just one list for the current stage. But not related to timings.
1. Lessons.xml
   + no changes
1. GitStrategy.xml
   + no changes

## Left over.

 + Start with BIG.R and TOY.R and CSV files. (b155cff)
 + Copy TOY.R to TOY2.R                      (ff082b9)
 + Move the functions from TOY2.R to functions_orig.R with the Vectorize() calls left in TOY2.R (c276f91)
 + Move the Vectorize calls from TOY2.R to functions.R  (b1889c8) (7564092)
 + Copy functions_orig.R to functions.R  (8df4ac2)
 + Comment out setwd() and rm() calls   (264c6184)
   + Timing: 586.48 seconds
 + Comment out fstar and fstarvalue in TOY2.R and make them local to accumulate() and
  firststageaccumulate() (GlobalVariables.xml)  (c8f0aa0b88)
   + Timing: 611.94 seconds. (and 624 seconds)
 + Initial vectorization of finalinflowprep. That function is 13-21 times faster than the Vectorize()
  version.  No timing on TOY2.R as finalinflow is only called once. o(b7a7964)
 + Improved version of finalinflow using lookup, corresponding to IfElse_Or_Lookup.xml. Again, no
  timing for TOY2.R for same reason. (76497cf)
