#############################################
###							###
###	Code for Monte Carlo parameter 	###
###   search and scores for the 		###
###   found parameter sets			###
###							### 
###			MK15				###
### 							###
### 	 Basolateral configuration		###
###							###
#############################################


rm(list=ls())	#clean memory
graphics.off() 	#close plot windows


############################# Basic model - this is basal ####################################

# run this code lines if you do not have the packges installed
# install.packages('deSolve')

library(deSolve)	#library for solve diferential equations
# citation("deSolve")

tmax=10000			# Insert the number of time units that the simulation will run
tstep=1			# Insert the step of the simulation
t=seq(0,tmax+1,tstep) 	# time

# Inicial states for variables
# insert, in the parentesis, the inicial value of the dependent variables (ex: X=2) 
state=c(PI=304372.0,
	PI3P=100.0054,
	PI4P=10040.13,
	PI5P=100.4866,
	PI35P2=20.21146,
	PI45P2=10000.97,
	PI34P2=20.61771,
	PI345P3=639.5332
	) 

# insert, in the parentecis, the parameters and independent variables and their values (ex: const=5, rate=2). 
# If no parameters and independent variables just put a zero.
parameters=c(
		pi_3KI = 74.71967,		
		pi_3KII = 4.997565,      	
		pi_3KII_III = 94.62836,		
		pi_4K = 152.4052,   		
		pi_Kfyve = 21.68522,		
		pip_5KI = 159.7755,		
		pip_5KII = 884.2539,	    	
		pi_4K_pip_5KI = 100.478,	

		SYNJ = 0.009059019,   		
		SYNJ_SAC1 = 39.78906,
		SYNJ_SAC1_SAC3 = 809.0106,
		SYNJ_SAC1_MTMR = 157.1438,
		SYNJ_TMEM55 = 63.7418,
		SIOSS = 35.34966,
		SIOSS_SHIP2 = 41.92146,
		MTMR = 87.86368,
		INPP4 = 213.8983,
		PTEN = 0.6764485,			 

		gamma1 = 16588.74, 
		gamma16 = 297.2636,
		gamma18 = 12.37069,
		gamma2 = 1.091975e+15, f2_PI = 1.138439e-01 , 		
		gamma3 = 3.067331e+12, f3_PI3P = 9.993362e-01,
		gamma4 = 5.100640e+14, f4_PI = 2.864618e-01,  
		gamma5 = 6.650532e+12, f5_PI4P = 9.093929e-01,
		gamma6 = 2.870688e+14, f6_PI = 6.388895e-02, 
		gamma7 = 3.067331e+12, f7_PI5P = 9.993362e-01,
		gamma8 = 1.33952e+14, f8_PI3P = 9.986181e-01,
		gamma9 = 3.060265e+12, f9_PI35P2 = 9.999336e-01,
		gamma10 = 8.493641e+15, f10_PI4P = 4.596175e-02,
		gamma11 = 6.650532e+12, f11_PI45P2 = 9.093929e-01,
		gamma12 = 2.946989e+13, f12_PI5P = 8.784401e-01,
		gamma13 = 4.133048e+12, f13_PI45P2 = 6.487218e-01,   	
		gamma14 = 1.888500e+14, f14_PI45P2 = 3.063327e-01,
		gamma15 = 4.814076e+14, f15_PI345P3 = 9.813754e-01,
		gamma17 = 4.15923e+14 , f17_PI35P2 = 9.998672e-01, 
		gamma19 = 1.333061e+12, f19_PI34P2 = 9.981983e-01,
		gamma21 = 1.678636e+11, f21_PI345P3 = 9.981984e-01,
		gamma_e = 0.051948,
		gamma30 = 2.666243e+14, f30_PI = 2.864618e-01,
		gamma31 =  1.045965e+13, f31_PI4P = 5.009150e-01,
		gamma32 = 7.597450e+12, f32_PI34P2 = 9.999769e-01,
		gamma33 = 6.650532e+12, f33_PI45P2 = 9.093929e-01						
		)  
Pars=as.data.frame(t(parameters)) #create a data.frame with the same information in vector parameters
attach(Pars) # put all colunms of Pars as variables in memory


equations=function(t,state,parameters) 	# function containing the diferential equations
	{with(as.list(c(state,parameters)),
		{

		# Fluxes		

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



		# insert the diferentical equations (ex: dX=const*X^rate or dX=5*X^2)

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





perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations
out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out)







#######################################
### 	Monte Carlo parameter seach 	###
### 						###
###	Do no run. The sets found can ###
###	be loaded in the end of the 	###
###	code block				###
#######################################
#
#
# run basic model please
#
#cbind(1:64,parameters)									# numerating the parameters
#
#library(rootSolve) 									# library for the jacobian function
#par_vault=matrix(parameters,nrow=1)							# new vector to store the parameters that are OK
#colnames(par_vault)=names(parameters)						# put the names of the parameters in the columns of the par_vault matrix
#par_trash=matrix(rep(NA,length(parameters)),nrow=1)				# matrix to store the parameters sets that are NOT OK
#Sstate_vault=matrix(tail(out[,2:9],1),nrow=1)					# new vector to store the steady states
#colnames(Sstate_vault)=colnames(out[,2:9])					# put the names of the dep. variables in the columns of the Sstate_vault matrix 
#SS_OK=vector()										# vector to store if the parameter set conplies with all the tests
#par_golden=matrix(parameters,nrow=1)						# matrix to store the parameter sets that pass the tests
#colnames(par_golden)=names(parameters)						# put names of the parameters to the par_golden matrix
#Sstate_golden=matrix(tail(out[,2:9],1),nrow=1)					# matrix to store the steady states for the parameter sets that pass the tests
#colnames(Sstate_golden)=colnames(out[,2:9])					# put names of dep. variables in the Sstate_golden matrix
#
#sample_size = 12500									# sample size 
#count1 = 0
#
#while (dim(par_vault)[1]<sample_size)	  					# cycling iterations
#	{
#	count = count + 1
#	writeLines(paste(dim(par_vault)[1]*100/sample_size,"%")) 		# percentage of iterations done so far  
#	ok=FALSE										# initializing variable ok. ok=TRUE if the parameter set is acceptable and FALSE otherwise
#
#	unif_min=0										# min of the variablility allowed to the parameters
#	unif_max=2										# max of the variablility allowed to the parameters
#
#	# mutating parameters
#	offspring=parameters 								# creating the vector were the mutations will occur
#	
#	# INs
#	offspring[19]= parameters [19] *runif(1,unif_min,unif_max)	      # gamma1
#	offspring[20]= parameters [20] *runif(1,unif_min,unif_max)		# gamma16
#	offspring[21]= parameters [21] *runif(1,unif_min,unif_max)		# gamma18
#
#	# OUTs
#	offspring[56]=parameters[56]* runif(1,unif_min,unif_max)		# out gamma_e
#
#	# others
#	offspring[29]=parameters[29]* runif(1,unif_min,unif_max)		# f5_PI4P
#	offspring[35]=parameters[35]* runif(1,unif_min,unif_max)		# f8_PI3P
#	offspring[43]=parameters[43]* runif(1,unif_min,unif_max)		# f12_PI5P
#	offspring[47]=parameters[47]* runif(1,unif_min,unif_max)		# f14_PI45P2
#	offspring[51]=parameters[51]* runif(1,unif_min,unif_max)		# f17_PI35P2
#	offspring[60]=parameters[60]* runif(1,unif_min,unif_max)		# f31_PI4P
#	offspring[41]=parameters[41]* runif(1,unif_min,unif_max)		# f11_PI45P2
#
#
#
#	out_MC=NA #cleanning out_MC vector
#
#	out_MC <- ode(y = state, times = t, func = equations,parms=offspring)	# good and old ode
#
#	Sstate=tail(out_MC,1)[,2:ncol(out_MC)]	
#	slopes=equations(tail(out_MC,1)[1],state=tail(out_MC,1)[,2:ncol(out_MC)],offspring)
#
#	# calculation of eigenvalues
#	Sstate_V=unlist(Sstate[1:9])
#	if (sum(is.na(Sstate_V))==0)
#		{
#		class(Sstate_V)
#		j=jacobian.full(Sstate_V,func=equations,parms=parameters)
#		if (exists('j')) { eigensum = sum(Re(eigen(j)$values[2:9])<=0) } else { eigensum = 0 }
#		} else {eigensum = 0}
#
#	# testing parameter sets to see if they are fit
#	if 	(
#		exists('out_MC')								# if the solution exists
#		&
#		!is.na(out_MC[nrow(out_MC),ncol(out_MC)])				# if the solution is not made of NAs
#		&
#		!is.nan(out_MC[nrow(out_MC),ncol(out_MC)])
#		&	
#		sum(abs(slopes$dy)<=1E-5)==8						# if the variation of the variables are small (steady state)
#		)
#		{	
#		ok=TRUE									# ok=TRUE if the parameter set is acceptable and FALSE otherwise
#		par_vault=rbind(par_vault,offspring)				# storing the parameter set if is fit
#		Sstate_vault=rbind(Sstate_vault,tail(out_MC[,2:9],1))		# storing the steady state if is fit
#		if 	(
#		(tail(out_MC,1)[2]>200000)&(tail(out_MC,1)[2]<400000)		# if 200000 < PI < 400000
#		&
#		(tail(out_MC,1)[3]>.1)&(tail(out_MC,1)[3]<6000)			# if .1 < PI3P < 6000
#		&
#		(tail(out_MC,1)[4]>5000)&(tail(out_MC,1)[4]<20000)		# if 5000 < PI4P < 20000
#		&
#		(tail(out_MC,1)[5]>.1)&(tail(out_MC,1)[5]<200)			# if .1 < PI5P < 200
#		&
#		(tail(out_MC,1)[6]>.1)&(tail(out_MC,1)[6]<200)			# if .1 < PI35P2 < 200
#		&
#		(tail(out_MC,1)[7]>5000)&(tail(out_MC,1)[7]<20000)		# if 5000 < PI45P2 < 20000
#		&
#		(tail(out_MC,1)[8]>.1)&(tail(out_MC,1)[8]<200)			# if .1 < PI34P2 < 200
#		&
#		(tail(out_MC,1)[9]>.1)&(tail(out_MC,1)[9]<5000)			# if .1 < PI345P3 < 5000
#		&
#		(abs(tail(out_MC,1)[4]-tail(out_MC,1)[7])<2000)			# if the diference between PI4P and PI45P2 is less than 2000
#		&
#		(abs(tail(out_MC,1)[3]-tail(out_MC,1)[5])<50)			# if the diference between PI3P and PI5P is less than 50
#		&
#		(sum(tail(out_MC,1)[2:9])>220000)&(sum(tail(out_MC,1)[2:9])<450000) # If the sum of all PIs is between 220000 and 450000
#		&
#		eigensum==8			   	             		     # if real part of eigenvalues is negative (Steady state stable)
#		)			
#			{
#			SS_OK=c(SS_OK,1)
#			par_golden=rbind(par_golden,offspring)				# storing the parameter set if is fit
#			Sstate_golden=rbind(Sstate_golden,tail(out_MC[,2:9],1))	# storing the steady state if is fit
#			} else {SS_OK=c(SS_OK,0)}
#
#		} else {par_trash=rbind(par_trash,offspring)}			# saving the trash
#	print(ok)										# sending ok variable to consule
#	writeLines(paste("Perfect solution? ",tail(SS_OK,1)))
#	writeLines(paste("Lists of parameters stored",nrow(par_vault)))   # sending the quantity of parameter sets stored so far
#	flush.console()									# send everything to the consule
#	}

#if (is.na(par_trash[1,1])) {par_trash=par_trash[2:nrow(par_trash),]}


#tail(par_vault,1)


# BIG_output=cbind(par_vault,Sstate_vault)
# write.table(BIG_output,file="BIG_sens_analysis_for_MK11") 		# write to disk Monte Carlo output (par_vault + Sstate_vault)
# MC2_output=read.table(file="BIG_sens_analysis_for_MK11") 			# read from disk Monte Carlo output (par_vault + Sstate_vault)

# save(par_vault,Sstate_vault,par_golden,Sstate_golden,par_trash,SS_OK,file = "parameter_space_explorer_MK15_1.RData")	
# load( "parameter_space_explorer_MK15_1.RData")

### loading the found parameter sets ###
load(file = "parameter_space_explorer_MK15_join3.RData")
par_vault = joint[[1]]
Sstate_vault = joint[[2]]
par_golden = joint[[3]]
Sstate_golden = joint[[4]]




###########################################################################
### More tests to select the parameter sets with the best fit (golden).	###
################################# more tests ##############################

par_golden[1,]
dim(par_golden)
Sstate_golden[1,]
dim(Sstate_golden)




sum(
	(Sstate_golden[,1]>200000)&(Sstate_golden[,1]<400000)		# if 200000 < PI < 400000
	&
	(Sstate_golden[,2]>.1)&(Sstate_golden[,2]<6000)			# if .1 < PI3P < 6000
	&
	(Sstate_golden[,3]>5000)&(Sstate_golden[,3]<20000)		# if 5000 < PI4P < 20000
	&
	(Sstate_golden[,4]>.1)&(Sstate_golden[,4]<200)			# if .1 < PI5P < 200
	&
	(Sstate_golden[,5]>.1)&(Sstate_golden[,5]<200)			# if .1 < PI35P2 < 200
	&
	(Sstate_golden[,6]>5000)&(Sstate_golden[,6]<20000)		# if 5000 < PI45P2 < 20000
	&
	(Sstate_golden[,7]>.1)&(Sstate_golden[,7]<200)			# if .1 < PI34P2 < 200
	&
	(Sstate_golden[,8]>.1)&(Sstate_golden[,8]<5000)			# if .1 < PI345P3 < 5000
	&


	(abs(Sstate_golden[,3]-Sstate_golden[,6])<2000)			# if the diference between PI4P and PI45P2 is less than 2000
	&
	(abs(Sstate_golden[,2]-Sstate_golden[,4])<50)			# if the diference between PI3P and PI5P is less than 50
	&
	(Sstate_golden[,6]/100 >= Sstate_golden[,2])			# if PI45P2 is 100 larger than PI3P
	&
	(Sstate_golden[,4]/4 >= Sstate_golden[,5])			# if PI5P is 5 times (4 times) larger than PI35P2
	&
	(Sstate_golden[,4] >= Sstate_golden[,6]*0.005) & (Sstate_golden[,4] <= Sstate_golden[,6]*0.02)		# PI5P between 0.5 and 2% of PI45P2
	&


	(((par_golden[,19]*100)/Sstate_golden[,1])>.1)&(((par_golden[,19]*100)/Sstate_golden[,1])<25)		# PI input between 0.1 and 30% of the PI pool
	&
	(((par_golden[,20]*100)/Sstate_golden[,3])>.1)&(((par_golden[,20]*100)/Sstate_golden[,3])<25)		# PI4P input between 0.1 and 30% of the PI4P pool	
	&
	(((par_golden[,21]*100)/Sstate_golden[,2])>.1)&(((par_golden[,21]*100)/Sstate_golden[,2])<25)		# PI3P input between 0.1 and 30% of the PI3P pool
	&


	(((par_golden[,56]*Sstate_golden[,1]*100)/Sstate_golden[,1])>=0) & (((par_golden[,56]*Sstate_golden[,1]*100)/Sstate_golden[,1])<=7)		# PI output between 0 and 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,2]*100)/Sstate_golden[,2])<7)		# PI3P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,3]*100)/Sstate_golden[,3])<7)		# PI4P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,4]*100)/Sstate_golden[,4])<7)		# PI5P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,5]*100)/Sstate_golden[,5])<7)		# PI35P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,6]*100)/Sstate_golden[,6])<7)		# PI45P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,7]*100)/Sstate_golden[,7])<7)		# PI34P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,8]*100)/Sstate_golden[,8])<7)		# PI345P3 output less than 7% of the respective pool
    )



par_golden_index=which(
	(Sstate_golden[,1]>200000)&(Sstate_golden[,1]<400000)		# if 200000 < PI < 400000
	&
	(Sstate_golden[,2]>.1)&(Sstate_golden[,2]<6000)			# if .1 < PI3P < 6000
	&
	(Sstate_golden[,3]>5000)&(Sstate_golden[,3]<20000)		# if 5000 < PI4P < 20000
	&
	(Sstate_golden[,4]>.1)&(Sstate_golden[,4]<200)			# if .1 < PI5P < 200
	&
	(Sstate_golden[,5]>.1)&(Sstate_golden[,5]<200)			# if .1 < PI35P2 < 200
	&
	(Sstate_golden[,6]>5000)&(Sstate_golden[,6]<20000)		# if 5000 < PI45P2 < 20000
	&
	(Sstate_golden[,7]>.1)&(Sstate_golden[,7]<200)			# if .1 < PI34P2 < 200
	&
	(Sstate_golden[,8]>.1)&(Sstate_golden[,8]<5000)			# if .1 < PI345P3 < 5000
	&


	(abs(Sstate_golden[,3]-Sstate_golden[,6])<2000)			# if the diference between PI4P and PI45P2 is less than 2000
	&
	(abs(Sstate_golden[,2]-Sstate_golden[,4])<50)			# if the diference between PI3P and PI5P is less than 50
	&
	(Sstate_golden[,6]/100 >= Sstate_golden[,2])			# if PI45P2 is 100 larger than PI3P
	&
	(Sstate_golden[,4]/4 >= Sstate_golden[,5])			# if PI5P is 5 times (4 times) larger than PI35P2
	&
	(Sstate_golden[,4] >= Sstate_golden[,6]*0.005) & (Sstate_golden[,4] <= Sstate_golden[,6]*0.02)		# PI5P between 0.5 and 2% of PI45P2
	&


	(((par_golden[,19]*100)/Sstate_golden[,1])>.1)&(((par_golden[,19]*100)/Sstate_golden[,1])<25)		# PI input between 0.1 and 30% of the PI pool
	&
	(((par_golden[,20]*100)/Sstate_golden[,3])>.1)&(((par_golden[,20]*100)/Sstate_golden[,3])<25)		# PI4P input between 0.1 and 30% of the PI4P pool	
	&
	(((par_golden[,21]*100)/Sstate_golden[,2])>.1)&(((par_golden[,21]*100)/Sstate_golden[,2])<25)		# PI3P input between 0.1 and 30% of the PI3P pool
	&


	(((par_golden[,56]*Sstate_golden[,1]*100)/Sstate_golden[,1])>=0) & (((par_golden[,56]*Sstate_golden[,1]*100)/Sstate_golden[,1])<=7)		# PI output between 0 and 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,2]*100)/Sstate_golden[,2])<7)		# PI3P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,3]*100)/Sstate_golden[,3])<7)		# PI4P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,4]*100)/Sstate_golden[,4])<7)		# PI5P output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,5]*100)/Sstate_golden[,5])<7)		# PI35P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,6]*100)/Sstate_golden[,6])<7)		# PI45P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,7]*100)/Sstate_golden[,7])<7)		# PI34P2 output less than 7% of the respective pool
	&
	(((par_golden[,56]*Sstate_golden[,8]*100)/Sstate_golden[,8])<7)		# PI345P3 output less than 7% of the respective pool
    )

par_golden_index

(par_super_golden=par_golden[par_golden_index,])
(Sstate_super_golden=Sstate_golden[par_golden_index,])

super_golden = list(par_super_golden,Sstate_super_golden)	# list with the best fit parameter sets

# save(super_golden,file = "parameter_space_explorer_MK15_super_golden.RData")
# load(file = "parameter_space_explorer_MK15_super_golden.RData")
# par_super_golden = super_golden[[1]]
# Sstate_super_golden = super_golden[[2]]





				 ################
##########################	Score	    #############################
				 ################


library(deSolve)		# library for solve diferential equations

# load(file = "parameter_space_explorer_MK15_super_golden.RData")
# par_super_golden = super_golden[[1]]
# Sstate_super_golden = super_golden[[2]]



# Insert the number of time units that the simulation will run
tmax=7000
# Insert the step of the simulation
tstep=1
t=seq(0,tmax+1,tstep) 	# time

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

		gamma1 = 15000, 
		gamma16 = 300,
		gamma18 = 10,

		gamma2 = 1.091975e+15, f2_PI = 1.138439e-01 , 		
		gamma3 = 3.067331e+12, f3_PI3P = 9.993362e-01,
		gamma4 = 5.100640e+14, f4_PI = 2.864618e-01,  
		gamma5 = 6.650532e+12, f5_PI4P = 9.093929e-01,
		gamma6 = 2.870688e+14, f6_PI = 6.388895e-02, 
		gamma7 = 3.067331e+12, f7_PI5P = 9.993362e-01,
		gamma8 = 1.33952e+14, f8_PI3P = 9.986181e-01,
		gamma9 = 3.060265e+12, f9_PI35P2 = 9.999336e-01,
		gamma10 = 8.493641e+15, f10_PI4P = 4.596175e-02,
		gamma11 = 6.650532e+12, f11_PI45P2 = 9.093929e-01,
		gamma12 = 2.946989e+13, f12_PI5P = 8.784401e-01,
		gamma13 = 4.133048e+12, f13_PI45P2 = 6.487218e-01,   	
		gamma14 = 1.888500e+14, f14_PI45P2 = 3.063327e-01,
		gamma15 = 4.814076e+14, f15_PI345P3 = 9.813754e-01,
		gamma17 = 4.15923e+14 , f17_PI35P2 = 9.998672e-01, 
		gamma19 = 1.333061e+12, f19_PI34P2 = 9.981983e-01,
		gamma21 = 1.678636e+11, f21_PI345P3 = 9.981984e-01,
		gamma_e = 0.045,
		gamma30 = 2.666243e+14, f30_PI = 2.864618e-01,
		gamma31 = 1.045965e+13, f31_PI4P = 5.009150e-01,
		gamma32 = 7.597450e+12, f32_PI34P2 = 9.999769e-01,
		gamma33 = 6.650532e+12, f33_PI45P2 = 9.093929e-01						
		)  

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

# Score(parameters,state,bs=base_score,Occam=FALSE,only_score=TRUE)


scores=vector()
for (i in 1:dim(par_super_golden)[1])
	{
	 scores=c(scores,Score(ind=par_super_golden[i,],Ss_ind=Sstate_super_golden[i,],bs=base_score,Occam=FALSE,only_score=TRUE))
	}

scores

# save(joint,super_golden,scores,file="parameter_space_explorer_MK15_join_supergolden_score.Rdata")
# load(file='parameter_space_explorer_MK15_join_supergolden_score.Rdata')





############
# Figure 5 #
############
# load(file='parameter_space_explorer_MK15_join_supergolden_score.Rdata')
boxplot(scores[2:length(scores)],horizontal=TRUE,ylim=c(0,110),range=0,
	 staplewex = 1,xlab="Scores",frame=FALSE,col='darkgrey')
points(20,1,type='p',col='cyan',pch="|",cex=3)
points(scores[1],1,type='p',col='blue',pch=17,cex=3)
points(0,1,type='p',col='black',pch="|",cex=3)
text(x=c(0,scores[1],fivenum(scores[2:length(scores)])), labels =c(0,round(scores[1],0),round(fivenum(scores[2:length(scores)]),0)), y=1.25, col=c('black','blue',rep('black',5)),cex=1.5)



