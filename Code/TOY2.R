if(!exists("UseInline"))
  UseInline = FALSE

###########
###table of contents
###########
#1. prep
####a. libraries
####b. system directories
#2. exogenous information
####a. var definitions
####b. user defined variables
####c. LookupTables
#######i. climate
#######ii. T|V
#######iii. x|R,V
#######iv. season characteristics|month|Stage 
#######v. N|x
#3. dp solver components
####a. approximations/discretizations
####b. state variables
####c. state transitions
####d. infeasible space 
####e. R|(V,y)
####f. keep R|min(x,f*) only if feasible except Oct else R|(V,y) (and R|min(N,N*{t+1}))
####g. calc end of period storage|(R,V,y)
####h. stores outcome variables (f*, x)
#4. dp solvers
####a. direct benefit function
####b. accumulative obj function
#5. dp model
####a. last stage
####b. intermediate stage
####c. first stage
#6. results
####a. best policy
####b. libraries
####c. 
###########

#profvis({

########
##prep
#######
####a. libraries
####b. set workspace
#rm(list=ls(all=TRUE)) #start with empty workspace
#setwd("C:/Users/leadams/Desktop/ShastaDam/DP")
getwd()
####c.global functions
#Un monthlist = month.name
#######
##exogenous information
#######

#### user defined variables
#####
NoofStages=17 #number of stages in program
bin=10^6 #storage discretization unit
#Un Tsep=52 #temperature at which layers were separated into "cold' and "warm" with file "ExploringInputDataEachVwVc"
#Un roundingdigits=6 #10^6
Vcinitial=bin #initial conditions
Vwinitial=0#need to check the month bin #initial conditions
Vinitial=Vcinitial+Vwinitial #total storage in reservoir
K=4*10^6 #reservoir storage capacity
DP=0.2*10^3 #water is no longer passing through the turbines
Max=4*10^6#bin#mround(K+4*10^6,bin)+2*10^6 #max reservor storage + max(Qin)=4*10^6
RcstarWinter=0 #need to put in model at later date
p=0.99
#releases are at end of period
######### approximations/discretizations
######
#state vars
Vc=seq(0,Max, bin) #sequence of Vc choices 
Vw=seq(0,Max, bin) #sequence of Vw choices
#Un V=Vc+Vw
L=seq(-0.5*bin, Max+0.5*bin, bin) #makes bins 
###state space

Voptions=t(Vdiscretizations(Vc,Vw)[c(1:2),])
colnames(Voptions)=c("Vc","Vw")
#apprximated states (vector)
Vcstates=Voptions[,1] #Vc is vector of discretized release possibliities
Vwstates=Voptions[,2]#Vw is vector of discretized release possibliities
Vstates=Vcstates+Vwstates

####action(decision) vars
#Rc=Vc #decision variables must equal bin discretizations to look forward and backward in time
nospillRc=Vc #seq(0,Max+4*10^6, bin)
nospillRw=Vw #Rc
#Un nospillR=nospillRc+nospillRw
#approximated set of action choices
nospillRdiscretizations=Vdiscretizations(nospillRc,nospillRw)
#Un nospillRoptions=t(Vdiscretizations(nospillRc,nospillRw)[c(1:2),])
#Un colnames(nospillRoptions)=c("nospillRc","nospillRw")
#rownames(Rdiscretizations)=c("Rw","Rc","R")
#decisoin space options (vector)
nospillRcdecs=nospillRdiscretizations[2,]#Vcstates
nospillRwdecs=nospillRdiscretizations[1,]#"Vwstates
#Un nospillRdecs=nospillRcdecs+nospillRwdecs

####d. non-release dependent LookupTables
#####
#######i. climate
###########

#joint probability of historical air temp and inflow
#synthetically calculatead per month based on pdf. each record of 0.01 occurs with 0.01 joint prob/month for each month
#setwd("C:/Users/leadams/Desktop/ShastaDam/Input_data")
climate=read.csv("climateinputs5probs.csv", stringsAsFactors = FALSE) #climate data generated with InflowCleaningv2 oct 19 2017
#setwd("C:/Users/leadams/Desktop/ShastaDam") #paper output goes to this file directory
#print(xtable(climate,type="latex"), file="yeartypes.tex") 
#probability of occurence
probs=c(0.01,0.1, 0.5, 0.9, 0.99) #probs

####change in air temp for mixed pool state transition
rawTa=climate[,5]
pn=length(climate[,5])/12 #number of probs in analysis calculated with inflowcleaningv2 in input_data
ifelse(length(probs)!=pn, print("PROBSERROR"),"")
deltaTa=matrix(0,nrow=length(rawTa),ncol=1) #creates matrix within which to store deltaTa per month for each historical climate record with prob p
for(i in (pn+1):length(rawTa)){ #calculates deltaTa
  deltaTa[i]=rawTa[i]-rawTa[i-pn]
}
for(i in pn:1){ #calculates deltaTa for the first observation 
  deltaTa[i]=rawTa[length(rawTa)-(i-1)]-rawTa[i] #digits=-3) 
}
deltaTa=mround(deltaTa,0.1) #uses mround function from global functions in this script

######monthly inflow


rawQ=climate[,6]
Q=mround(finalinflow(climate[,2],rawQ), bin)

###Markov Chains
ystates=probs

Lookupyprep=cbind(climate,Q,deltaTa)
Lookupy=Lookupyprep[,c(2,10,11,12)]


# Build the 2-way tables for column 3 and column 4 (Q and Ta)
# where rows are indexed by month and columns indexed by p corresponding
# to values in columns 1 and 2 respectively.
tmp.m = as.character(unique(Lookupy$month))
tmp.p = unique(Lookupy$p)
LookupyQ = matrix(0, length(tmp.m), length(tmp.p), dimnames = list(tmp.m, as.character(tmp.p)))
LookupyQ[ cbind(as.character(Lookupy$month), as.character(Lookupy$p)) ] = Lookupy$Q

LookupyTa = matrix(0, length(tmp.m), length(tmp.p), dimnames = list(tmp.m, as.character(tmp.p)))
LookupyTa[ cbind(as.character(Lookupy$month), as.character(Lookupy$p)) ] = Lookupy$deltaTa
rm(tmp.m, tmp.p)


#########
#######ii. T|V
##########
####get raw temp and volume data
#setwd("C:/Users/leadams/Desktop/ShastaDam/SDP")
all=read.csv("TwopoolAlldatav4.csv") #full reservoir storage and temp dataset generated from file exploringinputdataeachvwvc 2 Sept
#notes: all was created with 51F warm and cold pool stratification temperature
abbreviated=read.csv("GroupidExpectedVwVcTcTw.csv") #cleaned and analyzed/aggregated dataset for analysis 6 Sept

#setwd("C:/Users/leadams/Desktop/ShastaDam/DP") #run analysis on VwVc to TcTw from this directory
#Un groupingcorrelation=lm(all$Groupid~all$month) #checking to see if i can run the analysis by group id
#good correlation OK to group Vc and Vw by groupid
#Residual standard error: 22.91 on 217 degrees of freedom
#Multiple R-squared:  0.7583,	Adjusted R-squared:  0.7472 
#F-statistic: 68.09 on 10 and 217 DF,  p-value: < 2.2e-16

#order by storage volume size, first of Vw then Vc - hopefully this is relatively monotonic
OrderedVwVc=abbreviated[order(abbreviated$totalVw,-abbreviated$totalVc),] 
#Un d=OrderedVwVc 
observedVc=all[,100] #raw Vc
observedVw=all[,99] #raw Vw
ObservedTc=all[,98] #raw reservoir Tc
ObservedTw=all[,97] #raw reservoir Tw
#full raw table of temp to volume relationships w 51F warm and cool split
ObservedLookupTable=cbind(observedVc,ObservedTc,observedVw,ObservedTw) 

#######iii. season characteristics|month|Stage
###############
#gets month based on stage number

#gets climate season based on month


#####################
#####-finding max Qc and Qw to include spill
##and to estimate range of possible end of period storage before releases (e.g., Vc+Qc-Rc)

######### mixed season coldpool transition
######
#coefficients
NickelsSpringCoefficients=c(3.324,0, -0.372, 0.264) #intercept, fall bypass(removed rather than replaced by -0.855*Rc releases), spring air temp, spring volume
#km3, km3/km3, km3/changeC, km3/km3 1km=810.714 TAF
#current units needed AF, AF/AF, AF/F, AF/AF
NickelsWinterCoefficients=c(0.887, -0.264, 0.322,-0.253) #intercept, winter air temp, winter inflow, oct/nov reservoir temp
#km3, km3/changeC, km3/km3, km3/C
#AF, AF/F, AF/AF, AF/F
afconversion=810714 #AF
springcoeff=NickelsSpringCoefficients*c(afconversion,1, afconversion,1) 
wintercoeff=NickelsWinterCoefficients*c(afconversion,afconversion,1,afconversion) 

#####################
###including spill/transition states in model state and action spaces
#####################

#get discretizations with spill (remember to line up the index)
#get max available VC during fall overturn #calculations from fallsolve
MaxQcprep=max(ColdDelta("January",RcstarWinter, nospillRc, p),
              ColdDelta("February",RcstarWinter, nospillRc, p),
              ColdDelta("March",RcstarWinter, nospillRc, p),
              max(QLookup("November",p)+max(Vwstates)-min(nospillRw)+max(Vcstates)),
              ColdDelta("December",RcstarWinter, nospillRc, p))
MaxQw=max(Lookupy[Lookupy[,1]=="April" | Lookupy[,1]=="May" | Lookupy[,1]=="June" | Lookupy[,1]=="July" | Lookupy[,1]=="August" | Lookupy[,1]=="September" | Lookupy[,1]=="October",3])

##need max Q included in lookup temp table 
AvailVcprep=seq(0, MaxQcprep + max(Vc), bin)
AvailVw=seq(0,MaxQw+max(Vw),bin) #vc and Vw includes inflow and atm conditions

MaxQc=max(MaxQcprep,MaxQcprep+max(AvailVcprep)-max(AvailVw)) #fall overturn creates situations in which more than maxQ enters

#Un AvailVc=seq(0,MaxQc+max(Vc),bin)
#Un AvailVw=seq(0,MaxQw+max(Vw),bin) #vc and Vw includes inflow and atm conditions

#Voutoptions=Vdiscretizations(AvailVc,AvailVw)
#Voutoptions=t(Vdiscretizations(AvailVc,AvailVw)[c(1:2),])
#colnames(Voutoptions)=c("Vc","Vw")

########make V| T to include start and end of period storage possibilities (excluding releas options)
#aggregated into uniform bins of width "bin" for DP 


allbins=MakingBins(ObservedLookupTable,observedVc,observedVw,Vc,Vw,L)
preppeddata=allbins
preppeddata[is.na(preppeddata)]=0 #this includes observations for which there is some of one pool but not the other
#gets median cold pool temperature for each possible combination of cold and warm pool volumes
Tc=aggregate(as.numeric(preppeddata[,4])~Vcbin+Vwbin, 
             preppeddata, FUN=median, na.action=na.pass,na.rm=TRUE)
#gets median warm pool temperature for all observationseach possible combination of cold and warm pool volumes
Tw=aggregate(as.numeric(preppeddata[,6])~Vcbin+Vwbin, 
             preppeddata, FUN=median, na.action=na.pass,na.rm=TRUE)

LookupTableprep=merge(Tc,Tw,by=c("Vcbin","Vwbin"), all=TRUE)
colnames(LookupTableprep)=c("Vc","Vw","Tc","Tw")
LookupTableprep2=merge(Voptions,LookupTableprep, by=c("Vc","Vw"), all=TRUE) #get full range of DP states

#catches observations with volumes less than lowest discretization
LookupTableprep2[,3]=ifelse(LookupTableprep2[,3]==0,NA,LookupTableprep2[,3])
LookupTableprep2[,4]=ifelse(LookupTableprep2[,4]==0,NA,LookupTableprep2[,4])
LookupTable=LookupTableprep2


#either have all possible Vc and Vw functions in the lookup table or have a function that assumes
#anything above K etc as equivalent to a full Vc /Vw


#interpolate missing values either by carrying forward last observation or spline (polynomial interpolation)
LookupTableprep2[,3]=ifelse(LookupTableprep2[,3]>0, LookupTableprep2[,3],NA)
LookupTableprep2[,4]=ifelse(LookupTableprep2[,4]>0, LookupTableprep2[,4],NA)
#require(zoo)
#install.packages("zoo")
#library(zoo)

#the data was prepped to use this form of interpolation. all 0 are NA
#such that the data is carried forward always
LookupTableprep3 = zoo::na.locf(LookupTableprep2) #last observation carried forward 

#if vol is 0 then no release with 0 temp 
LookupTableprep3[,3]=ifelse(LookupTableprep3[,1]==0, 0, LookupTableprep3[,3])
LookupTableprep3[,4]=ifelse(LookupTableprep3[,2]==0, 0, LookupTableprep3[,4])
LookupTable=LookupTableprep3[order(LookupTableprep3[,1],LookupTableprep3[,2],LookupTableprep3[,3],LookupTableprep3[,4]),]



tmp.Vc = unique(LookupTable$Vc)
tmp.Vw = unique(LookupTable$Vc)
LookupTableTc = LookupTableTw = matrix(NA, length(tmp.Vc), length(tmp.Vw), dimnames = list(tmp.Vc, tmp.Vw))

LookupTableTc[cbind(as.character(LookupTable$Vc), as.character(LookupTable$Vw))] = LookupTable$Tc
LookupTableTw[cbind(as.character(LookupTable$Vc), as.character(LookupTable$Vw))] = LookupTable$Tw


greaterTc=LookupTable[LookupTable[,1]==max(Vc) & LookupTable[,2]==0,3]
greaterTw=LookupTable[LookupTable[,1]==0 & LookupTable[,2]==max(Vw),4]
##########get T from R and V


###plot RC and RW and TC and TW
 #plot(LookupVRT[,3],LookupVRT[,7],ylim=c(40,65)) #Rc v T
 #plot(LookupVRT[,4],LookupVRT[,7],ylim=c(40,65)) #Rw v T

#####add maximum spill potential option to spill 

#Rcprep=seq(0,MaxQcprep+max(nospillRc),bin) #releases include spill
Rc=seq(0,MaxQc+max(nospillRc),bin) #releases include spill
Rw=seq(0,MaxQw+max(nospillRw),bin)

#LookupVR=VRdiscretizations(AvailVc,AvailVw,Rc,Rw)
#LookupVRTraw=merge(LookupVR,LookupTable,by=c("Vc","Vw"), all=TRUE)

#######################
#######iii. x|(R,V)
##########################
#make table of Rc, Rw, Vc, Tc, Tw, month, x (of several policy choices,x, e.g. x from monthly thresholds and constant thresh)
#get x from T

##consistent monthly target of 56F


#Lookupx=cbind(ClrCk,Balls,Jelly)#,#Bend,
                    #RBDD)
#LookupVRTx=LookupVRTxsac

 #plot(LookupVRT$ReleaseT,ClrCk, ylim=c(45,65), xlim=c(45,65), type="l", ylab="",xlab="")
 #par(new=T)
 #plot(LookupVRT$ReleaseT,Balls,ylim=c(45,65),xlim=c(45,65),col="red", type="l", ylab="",xlab="")
 #par(new=T)
 #plot(LookupVRT$ReleaseT,Jelly,ylim=c(45,65),xlim=c(45,65),col="blue",type="l", ylab="",xlab="")
 #par(new=T)
 #plot(LookupVRT$ReleaseT,Bend,ylim=c(45,65), xlim=c(45,65),col="green",type="l", ylab="",xlab="")
 #par(new=T)
 #plot(LookupVRT$ReleaseT,RBDD,ylim=c(45,65), xlim=c(45,65),col="purple", type="l", ylab="",xlab="")

rivermiles=c(302,302,289,276,266#,#256,
             #243
             )#,39,37,6) #shasta dam river mile not provided but not important since mgmt starts with keswick
distance=302-rivermiles #distance from Shasta

#######
###fish temp reqs
#########


#LookupVRTx=cbind(LookupVRT,ConstantTempx)
 #plot(ReleaseT,ConstantTempx)
###############
#################
#######v. N|x
##############
#eventually add N to the x|(T,month,R,V) table

#######
##DP solver components
######
###########
####a. including spill in action space
#########
#approximated set of action choices
Rdiscretizations=Vdiscretizations(Rc,Rw)
#Un Roptions=t(Vdiscretizations(Rc,Rw)[c(1:2),])
#Un colnames(Roptions)=c("Rc","Rw")
#rownames(Rdiscretizations)=c("Rw","Rc","R")
#decision space options (vector)
Rcdecs=Rdiscretizations[2,]#Vcstates
Rwdecs=Rdiscretizations[1,]#"Vwstates
Rdecs=Rcdecs+Rwdecs
############
####b. state, action and outcome spaces
###########
#to create state and action spaces
basics=matrix(0, nrow=length(Vstates),ncol=length(Rdecs))
onestage=matrix(0, nrow=length(Vstates), ncol=length(Rdecs)) #*length(ystates)
allstages=array(0, dim=c(length(Vstates),length(Rdecs),length(ystates), NoofStages)) #length(ystates)
#dimnames(allstages)[[1]]=as.list(Vstates)
#dimnames(allstages)[[2]]=as.list(Rdecs)
#dimnames(allstages)[[3]]=as.list(ystates)
#Un stageslist=seq(1,NoofStages,1)
#dimnames(allstages)[[4]]=as.list(stageslist)

#state space (matrix)
#!VcSpace=apply(basics,2,function(x)Vcstates)
#!VwSpace=apply(basics,2,function(x)Vwstates)
#!VSpace=apply(basics,2,function(x)Vstates)

#choice space (matrix)
#!Rcspace=t(apply(basics,1,function(x)Rcdecs))
#!Rwspace=t(apply(basics,1,function(x)Rwdecs))
#!Rspace=t(apply(basics,1,function(x)Rdecs))

NR = length(Vstates)
NC = length(Rdecs)
VcSpace = matrix(Vcstates, NR, NC)
VwSpace = matrix(Vwstates,  NR, NC)
VSpace = matrix(Vstates,  NR, NC)

#choice space (matrix)
#! Changed
Rcspace=t(matrix(Rcdecs, NC, NR))
Rwspace=t(matrix(Rwdecs, NC, NR))
Rspace=t(matrix(Rdecs, NC, NR))

#outcome spaces (matrices)
##creates matrix for holding each stage calculation
#can include y

fstar = whichxstar = xstar = Rcstar = Rwstar = stagepolicy(Vstates,NoofStages)
#whichf    #could possible shorten code by saying whichxstar=f, for all below options
#!whichxstar=stagepolicy(Vstates,NoofStages)
#stores the position of the optimal f to find xstar
#!xstar=stagepolicy(Vstates,NoofStages) #river mile number
#!Rcstar=stagepolicy(Vstates,NoofStages)
#!Rwstar=stagepolicy(Vstates,NoofStages)
#############
####c. state transitions
#############
####d. end of period storage |(R,V,y)
############



##############

#######
##DP solvers
######
####a. direct benefit function
############



####each lake and climate season's benefits | end of month storage





if(UseInline) 
    source("inlineVariables.R")



#objective function/calculate current benefits



###############
####b. accumulative obj function
#######################

#######
##DP model
######

####a. last stage
#p=0.5
S=NoofStages
##############
  #for(S in 11:11){
month=monthcounter(S)

for(i in 1:pn){
  p=ystates[i]
  lastStage=matrix(choosesolve(month,VwSpace,VcSpace,Rcspace, Rwspace, Rspace,VSpace,RcstarWinter,K, DP,p)
                   ,nrow=length(Vcstates),ncol=length(Rdecs))
  #colnames(lastStage)=Rdecs
  #rownames(lastStage)=Vstates
  #  print(LastStage) checking results
  allstages[,,i,S]=lastStage
  #print(allstages[,,i,S])
}
LastStage=#ifelse(obj==11, 
  (allstages[,,1,S]+allstages[,,2,S]+allstages[,,3,S]+allstages[,,4,S]+allstages[,,5,S])/pn #,

#  print(LastStage) checking results
#store data
fstar[,S]=ifelse(apply(LastStage,1,max, na.rm=TRUE)<0, -9999, apply(LastStage,1,max, na.rm=TRUE))
whichxstar[,S]=apply(LastStage,1,which.max) 
Rcstar[,S]=ifelse(fstar[,S]<0, -9999, Rcdecs[whichxstar[,S]])
Rwstar[,S]=ifelse(fstar[,S]<0, -9999, Rwdecs[whichxstar[,S]])
#Un Rstar=Rcstar+Rwstar
 # }

#####################
####b. intermediate stage
#####################

Rcaccumstar= Rwaccumstar = matrix(NA,nrow=length(Vstates),ncol=length(Rdecs))
currentB = accumB = Vcoutdirect = Vwoutdirect = Vcoutacc = Vwoutacc = matrix(0,nrow=length(Vstates),ncol=length(Rdecs))

for(S in (NoofStages-1):2){
  month=monthcounter(S)
  for(i in 1:pn){
    p=ystates[i]
#!    currentB=matrix(choosesolve(month,VwSpace,VcSpace,Rcspace, Rwspace, Rspace,VSpace,RcstarWinter,K, DP,p),       nrow=length(Vstates),ncol=length(Rdecs))
    currentB[] = choosesolve(month,VwSpace,VcSpace,Rcspace, Rwspace, Rspace,VSpace,RcstarWinter,K, DP,p)

#!  accumB = matrix(accumulate(month,S,Vcstates,Vwstates,VcSpace,RcstarWinter,VwSpace,Rcspace,Rwspace,VSpace,Rspace,p), nrow=length(Vstates),ncol=length(Rdecs))

    accumB[] = accumulate(month, S, Vcstates,Vwstates,VcSpace,RcstarWinter,VwSpace,Rcspace,Rwspace,VSpace,Rspace,p)

    intstage = pmin(currentB, accumB)
  
  ##get Rcstar
  #if choose accumulate then need to pick R from accumulate not choose, unless R is infeasible from accumulate (in which case use R from choose)
  isaccum = ifelse(intstage < currentB, 1, NA)
  
  ##get Rc_t+1
  
  #1. get Vc and Vw out for this S
  #! Vcoutdirect=matrix(OutgoingVc(S,VcSpace, RcstarWinter,VwSpace,Rcspace,Rwspace,p),nrow=length(Vstates),ncol=length(Rdecs)) #gets end period storage VC
  #! Vwoutdirect=matrix(OutgoingVw(S,VwSpace,Rwspace,p),nrow=length(Vstates),ncol=length(Rdecs)) #gets end period storage Vw

    Vcoutdirect[] = OutgoingVc(S,VcSpace, RcstarWinter,VwSpace,Rcspace,Rwspace,p) #gets end period storage VC
    Vwoutdirect[] = OutgoingVw(S,VwSpace,Rwspace,p) #gets end period storage Vw
    
  #2. get Rc_t+1 for those combinations for which accB < direct B
      #rule out infeasible outs with directR


    Rcaccumstar[] = Rwaccumstar[] = NA
    #lookup Rc_t+1 for each x based on its start of period storage ==current period end of storage (Vc and Vw out), only for feasible R_t+1 options

    idx = which(!is.na(isaccum), TRUE)
    vc = Vcoutdirect[idx]
    vw = Vwoutdirect[idx]

    for(tmp1 in seq_len(nrow(idx))) {
        tmp = idx[tmp1, , drop  = FALSE]
        w = Vcstates == vc[tmp1] & Vwstates == vw[tmp1]
        Rcaccumstar[tmp[1], tmp[2]] = Rcstar[w, S+1]
        Rwaccumstar[tmp[1], tmp[2]] = Rwstar[w, S+1]        
    }    
    
  #colnames(Rcaccumstar)=Rcdecs
  #colnames(Rwaccumstar)=Rwdecs
  
  #3. use R direct t if R acc t+1 is infeasible
  #!Vcoutacc=matrix(OutgoingVc(S, VcSpace, 0, VwSpace, Rcaccumstar, Rwspace, p), nrow=length(Vstates), ncol=length(Rdecs))
  #!Vwoutacc=matrix(OutgoingVw(S,VwSpace,Rwaccumstar,p),nrow=length(Vstates),ncol=length(Rdecs))

   Vcoutacc[] = OutgoingVc(S, VcSpace, 0, VwSpace, Rcaccumstar, Rwspace, p)
   Vwoutacc[] = OutgoingVw(S,VwSpace,Rwaccumstar,p)

    Voutacc = Vcoutacc + Vwoutacc
    Voutacc[ Vcoutacc < 0 | Vwoutacc < 0 ] = NA
#!    Voutacc = ifelse(Vcoutacc<0, NA,
#!                    ifelse(Vwoutacc<0, NA,
#!                           Vcoutacc+Vwoutacc))
  
    #4. is R acc infeasible
    

    feasibleRcac = feasibleRwac = rep(NA, length(Voutacc))
    w = !( Vcoutacc < Rcspace | Vwoutacc < Rwspace | (Voutacc > K | Voutacc <DP))
    w = w & !is.na(w) # Handle 
    feasibleRcac[w] = Rcaccumstar[w]
    feasibleRwac[w] = Rwaccumstar[w]    
    
#! feasibleRcac=ifelse(Vcoutacc<Rcspace, NA,
#!                    ifelse(Vwoutacc<Rwspace,NA, 
#!                           ifelse(Voutacc > K | Voutacc <DP, NA, 
#!                                     Rcaccumstar)))    
#! feasibleRwac=ifelse(Vcoutacc<Rcspace, NA,
#!                     ifelse(Vwoutacc<Rwspace,NA, 
#!                            ifelse(Voutacc > K | Voutacc <DP, NA, 
#!                                   Rwaccumstar)))
   feasibleR=feasibleRcac*feasibleRwac
  
  #5. get final R
    #!  finalRc = ifelse(is.na(feasibleRcac),Rcspace, feasibleRcac)
    #! finalRw=ifelse(is.na(feasibleRwac),Rwspace, feasibleRwac)
    finalRc = feasibleRcac
    w = is.na(feasibleRcac)
    finalRc[w] = Rcspace[w]

    finalRw = feasibleRwac
    w = is.na(feasibleRwac)
    finalRw[w] = Rwspace[w]    
  
    finalR=finalRc+finalRw
  
  ###get final x with final R
  #finalx=matrix(choosesolve(month,VwSpace,VcSpace,finalRc,finalRw, finalR, VSpace, 0, K, DP, p), nrow=length(Vstates),ncol=length(Rdecs))
  
  allstages[,,i,S] = choosesolve(month, VwSpace, VcSpace, finalRc, finalRw, finalR, VSpace, 0, K, DP, p) # finalx
  #print(head(allstages[,,i,S]))[,1:20]
  }
  onestage=(allstages[,,1,S]+allstages[,,2,S]+allstages[,,3,S]+allstages[,,4,S]+allstages[,,5,S])/pn 
  #Could use
  #     onestage = apply(allstages[,,,S], c(1, 2), sum)/pn
  # or add at the end in the inner loop
  #   onestage = onestage + allstages[,,i,S]
  # But curiously the apply() is not faster. Slower by .2 of a second.
  
  #store data
  tmp = apply(onestage, 1, max, na.rm = TRUE)
  tmp[tmp < 0] = -9999
  fstar[,S] = tmp   # ifelse(apply(onestage,1,max, na.rm=TRUE)<0, -9999, apply(onestage,1,max, na.rm=TRUE))
  idx = whichxstar[,S] = apply(onestage, 1, which.max)  
  Rcstar[,S] = Rcdecs[idx]
  Rwstar[,S] = Rwdecs[idx]
}
########################
####c. first stage
#######################
S=1
month=monthcounter(S)
firststageholding=array(0, dim=c(1, #length(Vstates)==1
                                 length(Rdecs),
                                 pn)) #length(ystates)==pn
#1)) 
#dimnames(firststageholding)[[1]]=as.list(Vinitial)
#dimnames(firststageholding)[[2]]=as.list(Rdecs)
#dimnames(firststageholding)[[3]]=as.list(ystates)

for(i in 1:pn){
 p=ystates[i]
#!!!  fs=vector(length=length(Rdecs))
#!!!  fstarvalue=vector(length=length(Rdecs))
currentfirstB=matrix(choosesolve(month,Vwinitial, Vcinitial,Rcdecs, Rwdecs, Rdecs, Vinitial,RcstarWinter,K, DP,p),nrow=length(1),ncol=length(Rdecs))
accumfirstB=matrix(firststageaccumulate(month, S, Vcstates,Vwstates,Vcinitial,RcstarWinter,Vwinitial,Rcdecs,Rwdecs,Rdecs, Vinitial,p),
                   nrow=length(1), ncol=length(Rdecs))
firststageprep=matrix(pmin(currentfirstB,accumfirstB),nrow=length(1),ncol=length(Rdecs))
firststage=firststageprep

##get Rcstar
#if choose accumulate then need to pick R from accumulate not choose, unless R is infeasible from accumulate (in which case use R from choose)
isaccumfirst=ifelse(firststage<currentfirstB, 1, NA)

##get Rc_t+1

#1. get Vc and Vw out for this S
Vcoutdirectfirst = matrix(OutgoingVc(S,Vcinitial,RcstarWinter,Vwinitial,Rcdecs,Rwdecs,p),nrow=length(1),ncol=length(Rdecs)) #gets end period storage VC
Vwoutdirectfirst = matrix(OutgoingVw(S,Vwinitial,Rwdecs,p),nrow=length(1),ncol=length(Rdecs)) #gets end period storage Vw

#2. get Rc_t+1 for those combinations for which accB < direct B
#rule out infeasible outs with directR
Rcaccumstarfirst = Rwaccumstarfirst = matrix(NA,nrow=length(1),ncol=length(Rdecs))

#lookup Rc_t+1 for each x based on its start of period storage ==current period end of storage (Vc and Vw out), only for feasible R_t+1 options
 if(any(!is.na(isaccumfirst[1,]))) {
   for (r in 1:length(Rdecs)){
     if(!is.na(isaccumfirst[1,r])) {
         Rcaccumstarfirst[1,r]= Rcstar[which(Vcstates==Vcoutdirectfirst[1,r] & Vwstates==Vwoutdirectfirst[1,r]),S+1]
         Rwaccumstarfirst[1,r]= Rwstar[which(Vcstates==Vcoutdirectfirst[1,r] & Vwstates==Vwoutdirectfirst[1,r]),S+1]
     }
   }
  }


#3. use R direct t if R acc t+1 is infeasible
 Vcoutaccfirst=matrix(OutgoingVc(S, Vcinitial,0,Vwinitial,Rcaccumstarfirst,Rwdecs,p),nrow=length(1),ncol=length(Rdecs))
 Vwoutaccfirst=matrix(OutgoingVw(S,Vwinitial,Rwaccumstarfirst,p),nrow=length(1),ncol=length(Rdecs))
 Voutaccfirst= Vcoutaccfirst + Vwoutaccfirst
 Voutaccfirst[ Vcoutaccfirst<0 | Vwoutaccfirst<0] = NA

 #4. is R acc infeasible

  feasibleRcacfirst=Rcaccumstarfirst
 w = ( Vcoutaccfirst < Rcdecs| Vwoutaccfirst < Rwdecs | Voutaccfirst > K | Voutaccfirst <DP)
 w = is.na(w) | w
 feasibleRcacfirst[w] = NA
 feasibleRwacfirst = Rwaccumstarfirst
 feasibleRwacfirst[w] = NA

 feasibleRfirst=feasibleRcacfirst*feasibleRwacfirst

 #5. get final R
 finalRcfirst = feasibleRcacfirst
 finalRcfirst[w] = Rcdecs[w]
 finalRwfirst = feasibleRwacfirst
 finalRwfirst[w] = Rwdecs[w] 
 finalRfirst=finalRcfirst+finalRwfirst

    ###get final x with final R
     # finalxfirst=matrix(,nrow=length(1),ncol=length(Rdecs))
     #allstages[,,i,S]=firststageprep
  firststageholding[,,i] = choosesolve(month,Vwinitial,Vcinitial,finalRcfirst,finalRwfirst, finalRfirst, Vinitial, 0, K, DP, p)
} # end if for(i in 1:pn)


firststage=(
  firststageholding[,,1]
  + firststageholding[,,2]+ firststageholding[,,3]+ firststageholding[,,4]+ firststageholding[,,5])/pn 
fstarone=ifelse(max(firststage)<=0, -9999, max(firststage))
whichxstarone=which.max(firststage) 
Rcstarone=ifelse(fstarone<0, -9999, Rcdecs[whichxstarone])
Rwstarone=ifelse(fstarone<0, -9999, Rwdecs[whichxstarone])
#Rwstar[,1]=Rwstarone
#Rcstar[,1]=Rcstarone
xstarone=matrix(choosesolve(month,Vwinitial, Vcinitial,finalRcfirst, finalRwfirst, finalR, Vinitial,RcstarWinter,K, DP,p))[whichxstarone]

#6. results
####a. best policy
Best=matrix(0,nrow=NoofStages,ncol = 5)
colnames(Best)=c("Vc","Vw","Rc","Rw","x") # ,"month")
rownames(Best) = monthcounter(1:NoofStages)
Best[1,]=c(Vcinitial, Vwinitial,Rcstarone,Rwstarone,xstarone) # ,monthcounter(1))
rangeVc = rangeVw = rangex = numeric(pn)  
for(S in 2:NoofStages){
  ###for Vc 
  month=monthcounter(S-1)
  for(i in 1:pn){
    p=ystates[i]
    rangeVc[i]=OutgoingVcprep((S-1), Best[(S-1),1], RcstarWinter, Best[(S-1),2], Best[(S-1),3], Best[(S-1),4],p)
    rangeVw[i]=OutgoingVwprep((S-1), Best[(S-1),2], Best[(S-1),4], p)    
  }
  Best[S,1] = mround(sum(rangeVc)/pn,bin)
  ###for Vw
  Best[S,2] = mround(sum(rangeVw)/pn,bin)
  
  #month=monthcounter(S)
  #plocation=3 #p=0.5 is the third value
  idx = which(Vcstates== Best[S,1] & Vwstates == Best[S,2])
  Best[S,3] = Rcstar[idx, S]
  Best[S,4] = Rwstar[idx, S]
  month = monthcounter(S)

  for(i in 1:pn){
    p = ystates[i]
    # idx = which(Vcstates== Best[S,1] & Vwstates == Best[S,2])
    rangex[i]=matrix(choosesolve(month,VwSpace,VcSpace,Rcspace, Rwspace, Rspace,VSpace,RcstarWinter,K, DP,p)
                       ,nrow=length(Vcstates),ncol=length(Rdecs))[idx, whichxstar[idx, S]]

#    print(rangex[i])
  }
  Best[S,5] = sum(rangex)/pn
  # Best[S,6]=monthcounter(S)
  # print(S)
}
print(Best)
#write.csv(Best,"VariableTempRBDDToyxp01.csv")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#})

#file traken from 
#BASELINEBaselinePersistencev2cleanedDP_realdatav5wInflowsMONTHLYDPJaysdisplayUpdatedBenFxYearTypesv2
#jan 4 2018
##rewriting old code
#updating to run off several loopkup tables to reduce number of computations
#includes persistence functions
#could also re-write to run off different include files so that basic file is updated when choosing policies

#januar 29 2018
#rewriting lookup tables for V|R and improving efficiency of persistence constraint
#from DP_Code7VariableTempToy_persistencev4_EVv2
