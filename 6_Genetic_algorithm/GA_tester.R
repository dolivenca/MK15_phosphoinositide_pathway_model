#########################################################
###									###
###				MK15					###
###									###
### 	Genetic algoritm to search parameter space	###
###									###
#########################################################

rm(list=ls())	#clean memory
graphics.off() 	#close plot windows

### Basic model - basal

library(deSolve)	# library for solve diferential equations
library(units)	# library to do combinations
library(rootSolve) # to calculate Jacobians

#citation("deSolve")

# Insert the number of time units that the simulation will run
tmax=7000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

a=1
k=50*1.66*10^-18

# Inicial states for variables
# insert, in the parentesis, the inicial value of the dependent variables (ex: X=2) 
state=c(PI=319618.6,PI3P=92.42918,PI4P=10098.730,PI5P=110.29159,PI35P2=20.52992,PI45P2=9480.915,PI34P2=37.58107,PI345P3=763.1520) 

# insert, in the parentecis, the parameters and independent variables and their values (ex: const=5, rate=2). 
# If no parameters and independent variables just put a zero.

parameters=c(
		pi_3KI = 75 ,			# this and PTEN control the level of PI(345)P3
		pi_3KII = 5 ,      		# creates PI34P2 from PI4P. In vitro can phosphorilate PI to PI4P so effectly that the original name of this enzyme was PI_4K_type_I 
		pi_3KII_III = 90 ,		# Phosphorilate PI into PI3P.
		pi_4K = 360 ,   			# type II. Phosphorilate PI into PI4P. Very similar to PI_3Ks in structure. Type I was after all, PI_3KII 
		pi_Kfyve = 20 ,			# type III PIP_5K. Really diferent from the other PIP_5K.
		pip_5KI = 350 ,			# Phosphorilate PI4P into PI45P2. Activated by PA.
		pip_5KII = 750 ,    		# Phosphorilate PI5P into PI45P2 
		pi_4K_pip_5KI = 15 ,		#this is here because PI45P2 should not be inflenced by PI4P depletion

		SYNJ = 2 ,   			# phosphatase that transforms PI45P2 directly into PI
		SYNJ_SAC1 = 120 ,
		SYNJ_SAC1_SAC3 = 121 ,
		SYNJ_SAC1_MTMR = 200 ,
		SYNJ_TMEM55 = 50 ,
		SIOSS = 35 ,
		SIOSS_SHIP2 = 36 ,
		MTMR = 80 ,
		INPP4 = 50 ,
		PTEN = .5 ,				# (phosphatase and tensin homologue deleted on chromosome 10) PI_3P family. This and pip_3K control the level of PI(345)P3 

		gamma1 = a*15000, 
		gamma16 = a*300,
		gamma18 = a*10,

		gamma2 = a*1.091975e+15, f2_PI = 1.138439e-01 , 		
		gamma3 = a*3.067331e+12, f3_PI3P = 9.993362e-01,
		gamma4 = a*5.100640e+14, f4_PI = 2.864618e-01,  
		gamma5 = a*6.650532e+12, f5_PI4P = 9.093929e-01,
		gamma6 = a*2.870688e+14, f6_PI = 6.388895e-02, 
		gamma7 = a*3.067331e+12, f7_PI5P = 9.993362e-01,
		gamma8 = a*1.33952e+14, f8_PI3P = 9.986181e-01,
		gamma9 = a*3.060265e+12, f9_PI35P2 = 9.999336e-01,
		gamma10 = a*8.493641e+15, f10_PI4P = 4.596175e-02,
		gamma11 = a*6.650532e+12, f11_PI45P2 = 9.093929e-01,
		gamma12 = a*2.946989e+13, f12_PI5P = 8.784401e-01,
		gamma13 = a*4.133048e+12, f13_PI45P2 = 6.487218e-01,   	
		gamma14 = a*1.888500e+14, f14_PI45P2 = 3.063327e-01,
		gamma15 = a*4.814076e+14, f15_PI345P3 = 9.813754e-01,
		gamma17 = a*4.15923e+14 , f17_PI35P2 = 9.998672e-01, 
		gamma19 = a*1.333061e+12, f19_PI34P2 = 9.981983e-01,
		gamma21 = a*1.678636e+11, f21_PI345P3 = 9.981984e-01,
		gamma_e = a*0.045,
		gamma30 = a*2.666243e+14, f30_PI = 2.864618e-01,
		gamma31 = a* 1.045965e+13, f31_PI4P = 5.009150e-01,
		gamma32 = a*7.597450e+12, f32_PI34P2 = 9.999769e-01,
		gamma33 = a*6.650532e+12, f33_PI45P2 = 9.093929e-01						
		)  
Pars=as.data.frame(t(parameters)) #create a data.frame p with the same information in vector parameters
attach(Pars) # put all colunms of p as variables in memory


equations=function(t,state,parameters) 	# function containing the diferential equations
	{with(as.list(c(state,parameters)),
		{
		#rate of change (velocities of the variables concentration changes), the actual diferencial equations
		# insert the diferentical equations (ex: dX=const*X^rate or dX=5*X^2)
		
		V1 = gamma1
		V16 = gamma16
		V18 = gamma18

		V2 = gamma2*PI^f2_PI*pi_3KII_III* 142.185 * 1.66*10^-18
		V3 = gamma3*PI3P^f3_PI3P* (SYNJ_SAC1_MTMR) * 98.220 * 1.66*10^-18

		V4 = gamma4*PI^f4_PI* pi_4K * 109.711 * 1.66*10^-18
		V5 = gamma5*PI4P^f5_PI4P* (SYNJ_SAC1) * 135.202 * 1.66*10^-18

		V6 = gamma6*PI^f6_PI* pi_Kfyve * 237.136 *1.66*10^-18
		V7 = gamma7*PI5P^f7_PI5P* (SYNJ_SAC1) * 135.202 * 1.66*10^-18

		V8 = gamma8*PI3P^f8_PI3P* pi_Kfyve * 237.136 * 1.66*10^-18
		V9 = gamma9*PI35P2^f9_PI35P2* (SYNJ_SAC1_SAC3) * 127.311 * 1.66*10^-18

		V10 = gamma10*PI4P^f10_PI4P*pip_5KI * 65.643 * 1.66*10^-18
		V11 = gamma11*PI45P2^f11_PI45P2* (SIOSS) * 108.156 * 1.66*10^-18 

		V12 = gamma12*PI5P^f12_PI5P*pip_5KII * 46.968 * 1.66*10^-18
		V13 = gamma13*PI45P2^f13_PI45P2* (SYNJ_TMEM55) * 99.048 * 1.66*10^-18

		V14 = gamma14*PI45P2^f14_PI45P2*pi_3KI * 123.245 * 1.66*10^-18
		V15 = gamma15*PI345P3^f15_PI345P3* PTEN * 47.166 * 1.66*10^-18

		V17 = gamma17*PI35P2^f17_PI35P2* MTMR * 83.226 * 1.66*10^-18

		V19 = gamma19*PI34P2^f19_PI34P2* INPP4 * 107.347 * 1.66*10^-18

		V21 = gamma21*PI345P3^f21_PI345P3* (SIOSS_SHIP2) * 111.538 * 1.66*10^-18

		V22 = gamma_e*PI45P2
		V23 = gamma_e*PI
		V24 = gamma_e*PI4P
		V25 = gamma_e*PI345P3
		V26 = gamma_e*PI3P
		V27 = gamma_e*PI35P2
		V28 = gamma_e*PI5P
		V29 = gamma_e*PI34P2

		V30 = gamma30*PI^f30_PI* pi_4K_pip_5KI * 248.608 * 1.66*10^-18

		V31 = gamma31*PI4P^f31_PI4P* pi_3KII * 147.99 * 1.66*10^-18
		V32 = gamma32*PI34P2^f32_PI34P2* PTEN * 47.166 * 1.66*10^-18 

		V33 = gamma33*PI45P2^f33_PI45P2* (SYNJ) * 185.272 * 1.66*10^-18

		dPI = V1 + V3 + V5 + V7 + V33 - V2 - V4 - V6 - V23 - V30 
		dPI3P = V18 + V2 + V9 + V19 - V3 - V8 - V26
		dPI4P = V16 + V4 + V11 + V32 - V5 - V10 - V24 - V31
		dPI5P = V6 + V17 + V13 - V7 - V12 - V28
		dPI35P2 = V8 - V9	- V17 - V27
		dPI45P2 = V10 + V12 + V15 + V30 - V11 - V13 - V14 - V22 - V33
		dPI34P2 = V21 + V31 - V19 - V29 - V32
		dPI345P3 = V14 - V15 - V21 - V25		
	
		# return the rate of change 
		# insert, in the parentecis, the variables that you want to store the values (ex:dX)
		list(dy=c(dPI, dPI3P, dPI4P, dPI5P,dPI35P2,dPI45P2,dPI34P2,dPI345P3),
				count=c(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V21,V22,V23,V24,V25,V26,V27,V28,V29,V30,V31,V32,V33))
		})
	}





##################################### Perturb function ####################################

Perturb = function(ind)
	{

	perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations

	# Mega perturbation
	### PI45P2 sensitive to pi_4K, pi_4K_plus_pip_5KI.
	# Finaly remove the # form the front of all the next code lines
	# On the first column put the time of the perturbation
	ptime=c(0,
		1000,1500,
		2000,2005,2500,2505,
		3000,3005,3500,3505,
		4000,4500,
		5000,5500,6000,6500,
		max(t));

	# On the second column put the variables to be altered. (ex: 'X')
	pvar=c('novar',
		'gamma1','gamma1',
		'pi_4K','pi_4K_pip_5KI','pi_4K','pi_4K_pip_5KI',
		'pip_5KI','pi_4K_pip_5KI','pip_5KI','pi_4K_pip_5KI',
		'MTMR','MTMR',
		'pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve',
		'novar')
	
	# cbind(1:length(parameters),parameters)
	# On the third the new value for the variable.
	pval=as.numeric(c(NA,
				ind[19]*.5,ind[19],
				ind[4]*0.1,ind[8]*0.1,ind[4],ind[8],
				ind[6]*0.5,ind[8]*.5,ind[6],ind[8],
				ind[16]*.65,ind[16],
				ind[5]*.1,ind[5],ind[5]*.001,ind[5],
				NA))

	perturbations=data.frame(time=ptime,var=pvar,val=pval)
	perturbations
	}

# Perturb(parameters)





################################## Cruncher function - ode integrator that deals with perturbations ###############################

Cruncher = function (state,t,equations,parameters,perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA')))
	{

	# running the ode solver
	out_test=NA #cleanning out_MC vector

	# initializating the matrix out1 and out, necessary for the cycle
	out1=rbind(rep(0,1,length(state)+1),c(t=0,state))
	out_test=matrix()

	if (perturbations[1,1]=='NA') # if no perturbations an ordinary ODE is run
		{out_test <- ode(y = state, times = t, func = equations, parms=parameters)}

	if (perturbations[1,1]!='NA') # if there are perturbations, then go to the for cycle
		{
		if (tmax<max(perturbations[1:(nrow(perturbations)-1),1])) {cat("Time error: Last perturbation later than tmax.Please alter tmax or perturbation time.", "\n")}
		for (i in 1:(nrow(perturbations)-1)) # the number of perturbations is nrow(perturbations)-1 because the first line on the perturbations matrix is just to permit the inicial values to be used
			{
			t=seq(max(0,tail(out1[,1],1)+1),perturbations[i+1,1]-1,tstep)  	# set the time intervals #?#
	
			state=out1[nrow(out1),(2:(length(state)+1))]			# retrive the last values of the variables
			for (j in 1:length(state)) 						# check for variables perturbations
				{
				if (names(state)[j]==perturbations[i,2])			# if a perturbation is in a variable this will alter the variable value in the state vector
					{
					state[j]=perturbations[i,3]
					}
				}
			for (k in 1:length(parameters)) 						# check for variables perturbations
				{
				if (names(parameters)[k]==perturbations[i,2])			# if a perturbation is in a variable this will alter the variable value in the state vector
					{
					parameters[k]=perturbations[i,3]
					}
				}

			out1 <- ode(y = state, times = t, func = equations,parms=parameters)	# good and old ode
			if (is.na(out_test[1,1])) {out_test=out1} else {out_test=rbind(out_test,out1)} 			# store the values of the last interval	
			}
		names(out_test)[1]='time'
		}
	out_test
	}

# out_test = Cruncher(state,t,equations,parameters)						# no perturbations
# perturbations1 = Perturb(parameters)
# out_test = Cruncher(state,t,equations,parameters,perturbations1)			# with perturbations
# edit(out_test)





############### score of the inicial solution to normalize #################
# creating the base_score to normalize the socres of the diferent solutions. With this normalization each check will be 1. The total will be the number of checks.

perturbations = Perturb(parameters)
out_test = Cruncher(state,t,equations,parameters,perturbations)

base_score = c(
	abs(state[3]-state[6]),			# similar levels of PI4P and PI45P2
	abs(state[2]-state[4]),			# similar levels of PI3P and PI5P
	abs(state[4]-5*state[5]),		# PI5P levels are 5 fold of PI35P2 levels
	abs(state[3]-100*state[4]),	 	# PI4P or PI45P2 is 100 times more than PI5P or PI3P
	abs(state[1]-300000),			# Ss value of PI close to 300000
	abs(state[6]-10000),			# Ss value of PI45P2 close to 10000
	abs(state[2]-100),			# Ss value of PI3P close to 100
	abs(state[5]-state[7]),			# similar levels of PI35P2 and PI34P2

	abs(out_test[1490,2]/state[1]-out_test[1490,7]/state[6]),	# PI45P2 will decrease in the same percentage as PI

	abs(0.5-out_test[2490,4]/state[3]),	# PI4P drops to .5 after PI4K knockout
	abs(0.5-out_test[2490,7]/state[6]),	# PI45P2 drops to .5 after PI4K knockout
	abs(0.5-out_test[3490,7]/state[6]),	# PI45P2 drops to .5 after PIP_5KI knockdown

	abs(0.2-out_test[4490,5]/state[4]),	# PI5P drops to .2 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
	abs(1.5-out_test[4490,6]/state[5]),	# PI35P2 raizes to 1.5 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
	abs(0.5-out_test[5490,5]/state[4]),	# PI5P should drop to 50% if pi_Kfive is reduced to 10%
	abs(0.5-out_test[5490,6]/state[5]),	# PI35P2 should drop to 50% if pi_Kfive is reduced to 10%
	abs(0.15-out_test[6490,5]/state[4]),# PI5P should drop to .15 if pi_Kfyve is knockout
	abs(0.001-out_test[6490,6]/state[5]),# PI35P2 should drop to .001 if pi_Kfyve is knockout
	abs(0.8-out_test[6490,7]/state[6]),	# PI45P2 should drop to .8 if pi_Kfyve is knockout 
	abs(5-out_test[6490,3]/state[2]),	# PI3P should increase 5 fold if pi_Kfyve is knockout
	sum(parameters[1:18])			# Occam's razor
	)
names(base_score)=NULL
(base_score=as.numeric(base_score))	# make pval a numeric vector, the usual problem with data formats





########################### Mutation function ##########################
# parameters - parameter set to be mutated
# pop_size - number of mutatnts to be produced
# stand_dev - standard deviation of the normal that will cause the mutations
# mut_ins - if we will mutate the inffluxes
# mut_outs - if we will mutate the outfluxes
# num_par_mut - number of parameters to be mutated. Between 1 and 18.

Mutation = function (parameters,pop_size=1,stand_dev=1,mut_ins=TRUE,mut_outs=TRUE,num_par_mut=18)			# min of the variablility allowed to the parameters	# max of the variablility allowed to the parameters
	{
	population=matrix(rep(NA,length(parameters)*pop_size),ncol=length(parameters))						# new vector to store the population to be tested
	colnames(population)=names(parameters)						# put the names of the parameters in the columns of the population matrix
	population=as.data.frame(population)	

	for (ii in 1:pop_size)
		{
		# mutating parameters
		morula=parameters 								# creating the vector were the mutations will occur
		
		# INs
		if (mut_ins==TRUE) {		
			morula[19]= parameters [19] *max(0.001,rnorm(1,1,stand_dev))	      # gamma1
			morula[20]= parameters [20] *max(0.001,rnorm(1,1,stand_dev))		# gamma16
			morula[21]= parameters [21] *max(0.001,rnorm(1,1,stand_dev))		# gamma18
			}			 

		# OUTs
		if (mut_outs==TRUE) 
			{
			morula[56]=parameters[56]* max(0.001,rnorm(1,1,stand_dev))	# out gamma_e
			}	

		# others
		for (i in sample(1:18,num_par_mut))
			{
			morula[i]=parameters[i]* max(0.001,rnorm(1,1,stand_dev))
			}
		population[ii,]=morula
		}
	as.data.frame(population)
	}
# (population = Mutation (parameters, pop_size=10))





#########################   Viability function   ###############################

Viability = function (embryo)
	{
	out_test=NA #cleanning out_MC vector
	
	# Insert the number of time units that the simulation will run
	tmax=3000
	# Insert the step of the simulation
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time

	# initializating the matrix out1 and out, necessary for the cycle
	out1=rbind(rep(0,1,length(state)+1),c(t=0,state))
	out=matrix()
	
	out_v <- ode(y = state, times = t, func = equations, parms=embryo)
	# edit(out_test)

	Sstate_out<<-out_v[tmax,2:ncol(out_v)]
	Sstate<-out_v[tmax,2:ncol(out_v)]	
	slopes=equations(out_v[tmax,],state=out_v[tmax,2:ncol(out_v)],embryo)

	# calculation of eigenvalues
	Sstate_V=unlist(Sstate[1:9])
	if (sum(is.na(Sstate_V))==0)
		{
		class(Sstate_V)
		j=jacobian.full(Sstate_V,func=equations,parms=embryo)
		if (exists('j')) { eigensum = sum(Re(eigen(j)$values[2:9])<=0) } else { eigensum = 0 }
		} else {eigensum = 0}
	
	# testing parameter sets to see if they are fit
	if 	(
		exists('out_v')								# if the solution exists
		&
		!is.na(out_v[nrow(out_v),ncol(out_v)])				# if the solution is not made of NAs	
		&
		!is.nan(out_v[nrow(out_v),ncol(out_v)])
		&	
		sum(abs(slopes$dy)<=1E-5)==8						# if the variation of the variables are small (steady state)
		&
		(tail(out_v,1)[2]>200000)&(tail(out_v,1)[2]<400000)		# if 200000 < PI < 400000
		&
		(tail(out_v,1)[3]>.1)&(tail(out_v,1)[3]<6000)			# if .1 < PI3P < 6000
		&
		(tail(out_v,1)[4]>5000)&(tail(out_v,1)[4]<20000)		# if 5000 < PI4P < 20000
		&
		(tail(out_v,1)[5]>.1)&(tail(out_v,1)[5]<200)			# if .1 < PI5P < 200
 		&
		(tail(out_v,1)[6]>.1)&(tail(out_v,1)[6]<200)			# if .1 < PI35P2 < 200
		&
		(tail(out_v,1)[7]>5000)&(tail(out_v,1)[7]<20000)		# if 5000 < PI45P2 < 20000
		&
		(tail(out_v,1)[8]>.1)&(tail(out_v,1)[8]<200)			# if .1 < PI34P2 < 200
		&
		(tail(out_v,1)[9]>.1)&(tail(out_v,1)[9]<5000)			# if .1 < PI345P3 < 5000
		&
		(sum(tail(out_v,1)[2:9])>220000)&(sum(tail(out_v,1)[2:9])<450000) # If the sum of all PIs is between 220000 and 450000
		&
		eigensum==8			   	             		     # if real part of eigenvalues is negative (Steady state stable)
		)			
		{TRUE} else {FALSE}	
	}		
# Viability(parameters)





###########################   Score   ##################################

# ind = parameters
# Ss_ind = state
# bscore = base_score

Score = function (ind,Ss_ind,bscore=base_score,Occam=FALSE,only_score=FALSE)
	{
	# Insert the number of time units that the simulation will run
	tmax=7000
	# Insert the step of the simulation
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time

	names(Ss_ind)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3') # put names on steady states, if they do not have them.

	perturbations2 = Perturb(ind)
	out_score = Cruncher(Ss_ind,t,equations,ind,perturbations2)
	# edit(out_score)

	score = c(
			abs(Ss_ind[3]-Ss_ind[6])/bscore[1],				# similar levels of PI4P and PI45P2
			abs(Ss_ind[2]-Ss_ind[4])/bscore[2],				# similar levels of PI3P and PI5P
			abs(Ss_ind[4]-5*Ss_ind[5])/bscore[3],			# PI5P levels are 5 fold of PI35P2 levels
			abs(Ss_ind[3]-100*Ss_ind[4])/bscore[4],	 		# PI4P or PI45P2 is 100 times more than PI5P or PI3P
			abs(Ss_ind[1]-300000)/bscore[5],				# Ss value of PI close to 300000
			abs(Ss_ind[6]-10000)/bscore[6],				# Ss value of PI45P2 close to 10000
			abs(Ss_ind[2]-100)/bscore[7],					# Ss value of PI3P close to 100
			abs(Ss_ind[5]-Ss_ind[7])/bscore[8], 			# similar levels of PI35P2 and PI34P2

			abs(out_score[1490,2]/Ss_ind[1]-out_score[1490,7]/Ss_ind[6])/bscore[9],		# PI45P2 will decrease in the same percentage as PI

			abs(0.5-out_score[2490,4]/Ss_ind[3])/bscore[10],	# PI4P drops to .5 after PI4K knockout
			abs(0.5-out_score[2490,7]/Ss_ind[6])/bscore[11],	# PI45P2 drops to .5 after PI4K knockout
			abs(0.5-out_score[3490,7]/Ss_ind[6])/bscore[12],	# PI45P2 drops to .5 after PIP_5KI knockdown

			abs(0.2-out_score[4490,5]/Ss_ind[4])/bscore[13],	# PI5P drops to .2 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
			abs(1.5-out_score[4490,6]/Ss_ind[5])/bscore[14],	# PI35P2 raizes to 1.5 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)
			abs(0.5-out_score[5490,5]/Ss_ind[4])/bscore[15],	# PI5P should drop to 50% if pi_Kfive is reduced to 10%
			abs(0.5-out_score[5490,6]/Ss_ind[5])/bscore[16],	# PI35P2 should drop to 50% if pi_Kfive is reduced to 10%
			abs(0.15-out_score[6490,5]/Ss_ind[4])/bscore[17],	# PI5P should drop to .15 if pi_Kfyve is knockout
			abs(0.001-out_score[6490,6]/Ss_ind[5])/bscore[18],	# PI35P2 should drop to .001 if pi_Kfyve is knockout
			abs(0.8-out_score[6490,7]/Ss_ind[6])/bscore[19],	# PI45P2 should drop to .8 if pi_Kfyve is knockout 
			abs(5-out_score[6490,3]/Ss_ind[2])/bscore[20],		# PI3P should increase 5 fold if pi_Kfyve is knockout
			sum(ind[1:18])/base_score[21]					# Occam's razor
			)
	names(score)=NULL
	score=as.numeric(score)

	score_names=c(
	'Similar levels of PI4P and PI45P2',
	'Similar levels of PI3P and PI5P',
	'PI5P levels are 5 fold of PI35P2 levels',
	'PI4P or PI45P2 is 100 times more than PI5P or PI3P',
	'Ss value of PI close to 300000',
	'Ss value of PI45P2 close to 10000',
	'Ss value of PI3P close to 100',
	'Similar levels of PI35P2 and PI34P2',
	'PI45P2 will decrease in the same percentage as PI',
	'PI4P drops to .5 after PI4K knockout',
	'PI45P2 drops to .5 after PI4K knockout',
	'PI45P2 drops to .5 after PIP_5KI knockdown',
	'PI5P drops to .2 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)',
	'PI35P2 raizes to 1.5 after MTMR2 knockdown (MTMR estimated to be reduced to 65%)',
	'PI5P should drop to 50% if pi_Kfive is reduced to 10%',
	'PI35P2 should drop to 50% if pi_Kfive is reduced to 10%',
	'PI5P should drop to .15 if pi_Kfyve is knockout',
	'PI35P2 should drop to .001 if pi_Kfyve is knockout',
	'PI45P2 should drop to .8 if pi_Kfyve is knockout',	
	'PI3P should increase 5 fold if pi_Kfyve is knockout',
	"Occam's razor",
	"Score"
	)

	if (Occam==TRUE)												# if we want to calculate score with Occam's razor (less proteins is best)
		{
		score=c(score,sum(score))									# put the total score in the end of the score vector
		if (only_score==TRUE) {tail(score,1)} else {data.frame(score_names,score)}	# if I just want the total score put only_score=TRUE in the function parameters, else we will see the complete score vector
		} else {												# if we want to calculate score without Occam's razor 
		score[length(score)]=sum(score[1:(length(score)-1)])					# put the total score in the place of Occom's 
		score_names = score_names[-(length(score_names)-1)] 
		if (only_score==TRUE) {tail(score,1)} else {data.frame(score_names,score)}	# if I just want the total score put only_score=TRUE in the function parameters, else we will see the complete score vector
		}

	}

# Score(parameters,state,bs=base_score,Occam=FALSE,only_score=FALSE)





##########################   The Founding fathers   #####################

pop_size=10

population=matrix(rep(NA,length(parameters)*pop_size),nrow=pop_size)
population[1,]=as.numeric(parameters)
colnames(population)=names(parameters)
population=as.data.frame(population)
Ss_population=matrix(rep(NA,8*pop_size),nrow=pop_size)
Ss_population[1,]=state
colnames(Ss_population)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
	
par_trash=matrix(rep(NA,length(parameters)),nrow=1)				# matrix to store the parameters sets that are NOT OK
colnames(par_trash)=names(parameters)
par_trash=as.data.frame(par_trash)
counter2=2

while (counter2<=pop_size)
	{	
	Born = FALSE
	
	while (Born==FALSE)
		{
		embryo = Mutation (parameters, pop_size=1,stand_dev=.5)	
		if(Viability(embryo)==TRUE)			
			{
			population[counter2,] = embryo
			Ss_population[counter2,] = Sstate_out[1:8]
			Born = TRUE
			counter2=counter2+1
			writeLines(paste("Born? ",Born,'   ',round(100*(counter2-1)/(nrow(population)),2),'%'))
			} else {
				par_trash=rbind(par_trash,embryo)	# saving the trash
				writeLines(paste("Born? ",Born,'   ',round(100*(counter2-1)/(nrow(population)),2),'%'))
				}	
			flush.console()		

		}
	}
		
generation=0
scores=NA


# load("genetic_algoritm_MK15_temp_problemsolver2.RData")


# to create the best score variables. The initial values are not the best, just values to initiate the variables.
best=population[1,]
Ss_best=Ss_population[1,]
names(Ss_best)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
best_score=scores[1]



for (generation in 1:200) {

###########################   Parameter plot   ##########################

population
data1 = population[,c(1:18,19,20,21,56)]
data1[,19]=data1[,19]/100

windows() # quartz()
boxplot(data1,range=0,ylim=c(0,1500),col=c(rep('white',18),'blue',rep('white',3))
	 ,main='Population parameters evolution',cex.axis=.6,yaxt="n",las=2)
axis(2,cex.axis=1)
legend("topleft",
	c(paste('GEN',generation),paste('Best score = ',min(round(scores,1),na.rm=TRUE))),
	text.col=c('black','blue'),bty='n',cex=.5)	# legend to the plot, not used

my_solution= parameters[c(1:18,19,20,21,56)]
my_solution[19]=my_solution[19]/100
points(my_solution,col='orange')





###########################   Create offspring   ########################

pop_size=190		# desired population size. must be greater than 100

offspring=matrix(rep(NA,length(parameters)*pop_size),nrow=pop_size)
offspring=data.frame(offspring)
offspring[1,]=parameters
colnames(offspring)=names(parameters)

Ss_offspring=matrix(rep(NA,8*pop_size),nrow=pop_size)
Ss_offspring=data.frame(Ss_offspring)
Ss_offspring[1,]=state
colnames(Ss_offspring)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')

par_trash=matrix(rep(NA,length(parameters)),nrow=1)				# matrix to store the parameters sets that are NOT OK
colnames(par_trash)=names(parameters)
par_trash=as.data.frame(par_trash)

counter1=1

while (counter1<=choose(nrow(population),2))
	{
	#cbind(population[combn(nrow(population),2)[1,counter1],],
	#	population[combn(nrow(population),2)[2,counter1],],
	#	(population[combn(nrow(population),2)[1,counter1],]+population[combn(nrow(population),2)[2,counter1],])/2
	#	)
	embryo=(population[combn(nrow(population),2)[1,counter1],]+population[combn(nrow(population),2)[2,counter1],])/2

	Born = FALSE
	
	while (Born==FALSE)
		{	
		if(Viability(embryo)==TRUE)			
			{
			offspring[counter1,] = embryo
			Ss_offspring[counter1,] = Sstate_out[1:8]
			Born = TRUE
			counter1=counter1+1
			writeLines(paste("Born? ",Born,'   ',round(100*counter1/(choose(nrow(population),2)),2),'%'))
			} else {
				par_trash=rbind(par_trash,embryo)	# saving the trash
				embryo = Mutation(parameters=population[sample(1:10,1),],pop_size=1,stand_dev=0.3)
				writeLines(paste("Born? ",Born,'   ',round(100*counter1/(choose(nrow(population),2)),2),'%'))
				}	
			flush.console()		

		}
	}

population=rbind(population,offspring)
Ss_population=rbind(Ss_population,Ss_offspring)



writeLines("Mutantes")

counter3=counter1+10+1

while (counter3<=(counter1+10+1+10))
	{
	Born = FALSE
	while (Born==FALSE)
		{
		mutant=Mutation(parameters=population[sample(1:(counter1+10),1),],pop_size=1,stand_dev=0.3)	
		if(Viability(mutant)==TRUE)			
			{
			population[counter3,] = mutant
			Ss_population[counter3,] = Sstate_out[1:8]
			Born = TRUE
			counter3=counter3+1
			writeLines(paste("Born? ",Born,' || ',(counter3-(counter1+10+1)),' of 10.'))
			} else {
				writeLines(paste("Born? ",Born,' || ',(counter3-(counter1+10+1)),' of 10.'))
				}	
			flush.console()		
		}
	}



writeLines("Mutantes & Normals")

counter4=counter3+1

while (counter4<=100)
	{
	Born = FALSE
	while (Born==FALSE)
		{
		half_breed=(population[sample(1:55,1),]+population[sample(56:65,1),])/2	
		if(Viability(half_breed)==TRUE)			
			{
			population[counter4,] = half_breed
			Ss_population[counter4,] = Sstate_out[1:8]
			Born = TRUE
			counter4=counter4+1
			writeLines(paste("Born? ",Born,' || ',(counter4-(counter3+1)),' of 35.'))
			} else {
				writeLines(paste("Born? ",Born,' || ',(counter4-(counter3+1)),' of 35.'))
				}	
			flush.console()		
		}
	}


writeLines("Minor Mutations")

counter5=counter4+1

while (counter5<=(pop_size+10))
	{
	Born = FALSE
	while (Born==FALSE)
		{
		minor_mut=Mutation(parameters=population[sample(1:100,1),],pop_size=1,stand_dev=3,mut_ins=FALSE,mut_outs=FALSE,num_par_mut=sample(1:5,1))	
		if(Viability(minor_mut)==TRUE)			
			{
			population[counter5,] = minor_mut
			Ss_population[counter5,] = Sstate_out[1:8]
			Born = TRUE
			counter5=counter5+1
			writeLines(paste("Born? ",Born,' || ',(counter5-(counter4+1)),' of ',(pop_size+10-100)))
			} else {
				writeLines(paste("Born? ",Born,' || ',(counter5-(counter4+1)),' of ',(pop_size+10-100)))
				}	
			flush.console()		
		}
	}



scores=rep('NA',dim(population)[1])
for (iii in 1:dim(population)[1])
	{
 	Ss=as.numeric(Ss_population[iii,])
 	names(Ss)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3') 
	scores[iii]=Score(ind=population[iii,],Ss_ind=Ss,Occam=FALSE,only_score=TRUE)
	writeLines(paste(scores[iii],'   ',100*iii/(dim(population)[1]),'%'))
	flush.console()
	}
scores=as.numeric(scores)
scores[is.na(scores)]=1000

ordered=sort(scores,index.return = TRUE,)$ix

### save the best population in the GA that kills the parents ###
if (best_score>scores[ordered[1]])
	{
	best=population[ordered[1],]
	colnames(best) = names(parameters)
	Ss_best=Ss_population[ordered[1],]
	colnames(Ss_best)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
	best_score=scores[ordered[1]]
	}


new_pop=data.frame(matrix(rep(NA,length(parameters)*10),nrow=10))
new_pop[1,] = population[ordered[1],]
colnames(new_pop) = names(parameters)
Ss_new_pop = data.frame(matrix(rep(NA,8*10),nrow=10))
Ss_new_pop[1,]=Ss_population[ordered[1],]
colnames(Ss_new_pop)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
new_scores=scores[ordered[1]]
counter6 = 2
counter7 = 1
while (counter7<=9)
	{
	temp1=as.numeric(new_pop[counter7,])
	temp2=as.numeric(population[ordered[counter6],])
	if (identical(temp1,temp2)) 
		{
		counter6=counter6+1
		} else 
		{
		counter7=counter7+1
		new_pop[counter7,] = population[ordered[counter6],]
		Ss_new_pop[counter7,] = Ss_population[ordered[counter6],]
		new_scores=c(new_scores,scores[ordered[counter6]])
		counter6=counter6+1	
		}		
	}
population=new_pop
Ss_population=Ss_new_pop
scores=new_scores


### script to verify steady states ###
new_Ss=data.frame(matrix(rep(NA,10*8),nrow=10))
for (i in 1:10)
	{
	tmax=7000
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time

	p=population[i,]
 	Ss=as.numeric(Ss_population[i,])
 	names(Ss)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
	out1=Cruncher(Ss,t,equations,p)
	new_Ss[i,]=tail(out1,1)[2:9]
	}
colnames(new_Ss)=names(state)
Ss_population=new_Ss



save(population,Ss_population,scores,file = "genetic_algoritm_MK15_temp.RData")	
# load( "genetic_algoritm_MK15_temp_12_28.RData")
}







# save(population,Ss_population,scores,file = "genetic_algoritm_MK15_G30.RData")	
# load( "genetic_algoritm_MK15_G30.RData")











######################################### Test_me function #####################################

# ind = parameters
# Ss_ind = state

Test_me = function (ind,Ss_ind)
	{

	# Insert the number of time units that the simulation will run
	tmax=7000
	# Insert the step of the simulation
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time

	perturbations3 = Perturb(ind)
	out_test_me = Cruncher(Ss_ind,t,equations,ind,perturbations3)
	
	# edit(out_test_me)

	varnames=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
	varcolors=c('gray','darkgoldenrod1','coral','coral1','coral2','orange','orange2','darkorange3','gold','deepskyblue1','deepskyblue3')
	
	# normalized graphs
	stst=out_test_me[500,]					# getting steady state at time 500
	out_n=cbind(time=out_test_me[1],out_test_me[,2:ncol(out_test_me)])
	tail(out_n)
	out_norm=out_test_me[,1]				# creating matrix out_norm to store the normalized information
	for (i in 2:ncol(out_test_me))
		{
		newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
		out_norm=cbind(out_norm,newc)			# storing normalized information
		}
	colnames(out_norm)=colnames(out_test_me)
	head(out_norm)

	windows()
	par(mfrow=c(3,3))
	for (i in 2:9) 
		{
		plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
		text(tmax-2,out_norm[tmax-2,i],varnames[i],col=varcolors[i])
		}
	par(mfrow=c(1,1))

	# check for PI4P and PI45P2 up's and down's 
	check_pert = c('PI 50%',NA,'pi4K 10%','pi4K 10%',NA,'pip_5KI')   	
	check_names = c('PI45P2',NA,'PI4P','PI45P2',NA,'PI45P2')
	check_targets = c(.5,NA,.5,.5,NA,.5)
	check_darts = c(out_norm[1490,7],NA,out_norm[2490,4],out_norm[2490,7],NA,out_norm[3490,7])
	names(check_darts)=NULL
	check_big = data.frame(check_pert,check_names,check_targets,check_darts) 

	# check for small PIs
	check_pert = c('MTMR 10%','MTMR 10%',NA,'pi_Kfyve 10%','pi_Kfyve 10%',NA,'pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%')   	
	check_names = c('PI5P','PI35P2',NA,'PI5P','PI35P2',NA,'PI5P','PI35P2','PI45P2','PI3P')
	check_targets = c(.2,1.5,NA,.5,.5,NA,.15,.001,.80,5)
	check_darts = c(out_norm[4490,5],out_norm[4490,6],NA,out_norm[5490,5],out_norm[5490,6],NA,out_norm[6490,5],out_norm[6490,6],out_norm[6490,7],out_norm[6490,3])
	names(check_darts)=NULL
	check_small = data.frame(check_pert,check_names,check_targets,check_darts) 

	test_me_out = list (check_big,check_small,Ss_ind)
	test_me_out
	}

# Test_me(parameters,state)

# Score(parameters,state,bs=base_score,Occam=FALSE,only_score=FALSE)




### commands to test solutions ###
 p=population[1,]
 Ss=as.numeric(Ss_population[1,])
 names(Ss)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3') 
 Score(ind=p,Ss_ind=Ss,bscore=base_score,only_score=FALSE)
 Test_me(p,Ss)






### Script top add elements to the founding fathers ###

# load("GA_MK15_foundind_fathers_no_reduc.RData")
# pop=population
# Ss=Ss_population

# load( "genetic_algoritm_MK15_temp1.RData")
# pop[1,]=population[1,]
# Ss[1,]=Ss_population[1,] 
# population=pop
# Ss_population=Ss

scores2=rep('NA',dim(population)[1])
for (iii in 1:dim(population)[1])
	{
	pop=as.numeric(population[iii,])
	names(pop)=names(parameters)
 	Ss=as.numeric(Ss_population[iii,])
 	names(Ss)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3') 
	scores2[iii]=Score(ind=pop,Ss_ind=Ss,Occam=FALSE,only_score=TRUE)
	writeLines(paste(scores2[iii],'   ',100*iii/(dim(population)[1]),'%'))
	flush.console()
	}
scores2=as.numeric(scores2)

# save(population,Ss_population,scores,file = "GA_MK15_foundind_fathers_no_reduc.RData")	






### script to verify steady states ###
new_Ss=data.frame(matrix(rep(NA,10*8),nrow=10))
for (i in 1:10)
	{
	tmax=7000
	tstep=1
	t=seq(0,tmax+1,tstep) 	# time

	p=population[i,]
 	Ss=as.numeric(Ss_population[i,])
 	names(Ss)=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
	out1=Cruncher(Ss,t,equations,p)
	new_Ss[i,]=tail(out1,1)[2:9]
	}
Ss_population[10,]==new_Ss[10,]

#colnames(new_Ss)=names(state)
#Ss_population=new_Ss





