# Monte-Carlo simulation MK2 for pip model MK11
# Global sensitivites for MK11
# three steps : finding list of parameters, finding sensitivities for these parameters, counting the times each parameter has a significant sensitivity.
# should run a basic model first



# rm(list=ls())	#clean memory
#graphics.off() 	#close plot windows
#
#
#library(beepr)
#library(rootSolve) 									# library for the jacobian function
#par_vault=matrix(parameters,nrow=1)							# new vector to store the parameters that are OK
#colnames(par_vault)=names(parameters)
#par_trash=matrix(rep(NA,length(parameters)),nrow=1)
#Sstate_vault=matrix(tail(out[,2:9],1),nrow=1)					# new vector to store the steady states
#colnames(Sstate_vault)=colnames(out[,2:9])
#
#
#sample_size = 20									# sample size 
#
#
#while (dim(par_vault)[1]<sample_size)	  					# cycling iterations
#	{
#	writeLines(paste(dim(par_vault)[1]*100/sample_size,"%")) 		# percentage of iterations done so far  
#	ok=FALSE										# initializing variable ok. ok=TRUE if the parameter set is acceptable and FALSE otherwise
#	offspring=parameters * runif(length(parameters),.5,1.5)					# mutating parameters
#
#	out_MC=NA #cleanning out_MC vector
#
#	out_MC <- ode(y = state, times = t, func = equations,parms=offspring)	# good and old ode
#
#	Sstate=tail(out_MC,1)[,2:ncol(out_MC)]	
#	slopes=equations(tail(out_MC,1)[1],state=tail(out_MC,1)[,2:ncol(out_MC)],offspring)
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
#		}
#		else {par_trash=rbind(par_trash,offspring)}			# saving the trash
#	print(ok)										# sending ok variable to consule
#	writeLines(paste("Lists of parameters stored",nrow(par_vault)))   # sending the quantity of parameter sets stored so far
#	flush.console()									# send everything to the consule
#	}
#
#beep(1)
#
#
#if (is.na(par_trash[1,1])) {par_trash=par_trash[2:nrow(par_trash),]}
#
#
#tail(par_vault,1)
#tail(Sstate_vault,1)
#dim(par_vault)
#
#save and load 
# as table 
#BIG_output=cbind(par_vault,Sstate_vault)
# write.table(BIG_output,file="BIG_sens_analysis_for_MK15") 		# write to disk Monte Carlo output (par_vault + Sstate_vault)
# MC2_output=read.table(file="BIG_sens_analysis_for_MK15") 			# read from disk Monte Carlo output (par_vault + Sstate_vault)
#
# save and load
# as variable
# save(BIG_output,par_trash,file = "BIG_sens_analysis_for_MK15_50p_extra.RData")	
# load( "BIG_sens_analysis_for_MK15_50p_extra.RData")





############################################# calculate sensitivities of the diferent solutions #########################################
#
#
#library(beepr)
#library(deSolve)	#library for solve diferential equations
#
# load("6model_components.RData")	# read components for the model in disk
#
# load("BIG_sens_analysis_for_MK11.RData")	# read from disk
#ls()					# see what variables are in memory
#
#dim(BIG_output)
#par=BIG_output[,1:length(parameters)]		# read parameters from disk. Can be MC1 or MC2.
#SstateFromPar=BIG_output[,(length(parameters)+1):(length(parameters)+8)]	# read steady state from disk
#
#MSA_raw_list=list()
#MSA_simple_list=list()
#
#
#
#for (i in 1:nrow(par))
#	{
#	param=par[i,]
#	Ss_old=as.numeric(SstateFromPar[i,])
#	names(Ss_old)=colnames(SstateFromPar)
#
#	Sens=matrix(rep(NA,length(par[1,])*length(Ss_old)),ncol=8)	# creating response matrix
#	colnames(Sens)=names(Ss_old)							# naming the columns after the dependent variables
#	rownames(Sens)=names(par[1,])						# naming the rows after the parameters
#
#	for (ii in 1:length(param))						# cycling the parameters
#		{
#		Pars=param								# creating discardable parameter vector to be changed
#		Pars[ii]=param[ii]*(1.01)					# changing one parameter in the discardable parameter vector in 1%
#		out_new <- ode(y = Ss_old, times = t, func = equations, parms=Pars)		#runing the ode solver
#		Ss_new=tail(out_new,1)[,2:9]					# new steady state
#		tsens=vector()							# creating sensitivity vector for this turn 
#		for (iii in 1:8)							# cycling all dependent variables
#			{
#			tsens=c(tsens,(Ss_new[iii]-Ss_old[iii])*100/Ss_old[iii])	#calculating sensitivities
#			}
#		Sens[ii,]=tsens							# puting this parameter sensitivities in sentitivity matrix
#		}	
#	MSA_raw_list[[length(MSA_raw_list)+1]] <- Sens			# add a matrix of sensitivities to the list MSA_raw_list
#
#	sum(abs(Sens)>=1)								# amount of significant sensibilities (greater then 1%)
#	Sens_show=abs(Sens)>=1							# mark significant sensibilities 
#	dim(Sens_show)								# dimentions of sensitivity matrix
#	for(iv in 1:length(parameters))					# cycling the rows of the sensitivity matrix
#		{	
#		for(v in 1:8)								# cycling the columns of the sensitivity matrix
#			{
#			if (Sens_show[iv,v]==1 & !is.na(Sens_show[iv,v])==TRUE){Sens_show[iv,v]=Sens[iv,v]}	# substitute the TRUEs for the significant sensibilities
#			}
#		}
#
#	MSA_simple_list[[length(MSA_simple_list)+1]] <- Sens_show	# add a matrix of sensitivities greater than 1 to the list MSA_simple_list
#
#	writeLines(paste(i*100/nrow(par),"%"))		# percentage of iterations done so far  
#	flush.console()						# send everything to the consule
#	}
#beep(1)
#
# third try for storage of parameters and steady states
# save(BIG_output,MSA_raw_list,MSA_simple_list,file = "MSA_for_MK15_50p_extra.RData")		# write to disk the sensitivities and the parameters
# load("MSA_for_MK15_50p.RData")										# read from disk





######################### Massive Sensitivity analysis counts #####################################
########################### Used to clean the data from NAs #######################################
#
#
#
# load("MSA_for_MK15_100p.RData")
#
#length(MSA_raw_list)
# analysis of MSA
#BIG_scrap=matrix(rep(NA,72),nrow=1)
#
#
#
#Check=FALSE
#while (Check==FALSE)
#{
# use to clean
# MSA_raw_list[[vi]] = NULL; BIG_scrap=rbind(BIG_scrap,BIG_output[vi,]); BIG_output=BIG_output[-vi,]			# cleaning the data from NAs
# use to clean
#	
#nSol=length(MSA_raw_list)								# number of solutions to be analysed
#analysis_MSA=matrix(rep(NA,dim(MSA_raw_list[[1]])[1]*dim(MSA_raw_list[[1]])[2]),ncol=dim(MSA_raw_list[[1]])[2])						# initializing analysis matrix
#colnames(analysis_MSA)=colnames(MSA_raw_list[[1]])				# names for colunms
#rownames(analysis_MSA)=rownames(MSA_raw_list[[1]])				# names for rows
#
#for(v in 1:dim(MSA_raw_list[[1]])[1])						# cycling the parameters
#	{
#	for(iv in 1:dim(MSA_raw_list[[1]])[2])					# cycling the dependent variables
#		{
#		sumvar1=0
#		for (vi in 1:nSol)							# cycling the solutions
#			{ 
#			if (abs(MSA_raw_list[[vi]][v,iv])>=1) {sumvar1=sumvar1+1}
#			}
#		analysis_MSA[v,iv]=sumvar1	
#		}
#	}
#analysis_MSA
#
#
#MSA_raw_list[[vi]] 
#vi
#Check=!is.na(analysis_MSA[64,8])
#}


#length(MSA_raw_list)
#dim(BIG_scrap)


# third try for storage of parameters and steady states
# save(BIG_output,BIG_scrap,MSA_raw_list,MSA_simple_list,file = "MSA_for_MK15_100p_clean_cut.RData")			# write to disk
# load("MSA_for_MK15_100p_clean_cut.RData")												# read from disk


#dim(BIG_output)
#length(MSA_raw_list)
#length(MSA_simple_list)
# MSA_raw_list=c(MSA_raw_list,MSA_raw_list_e)
# MSA_simple_list=c(MSA_simple_list,MSA_simple_list_e)








































##############
# Figure S1a #
######################## condencing the information ###################


### Max min finder ###

	load("MSA_for_MK15_100p_clean_cut.RData")
#	load("MSA_for_MK15_50p_clean.RData")	
#	load("MSA_for_MK15_20p_clean.RData")	
#	load("MSA_for_MK15_10p_clean.RData")
#	load("MSA_for_MK15_5p_clean.RData")		
	
MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_max=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_max)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_max)=rownames(MSA_raw_for_ext[[1]])
sens_min=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_min)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_min)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_max[i,ii]=max(x_temp)
		sens_min[i,ii]=min(x_temp)
		}
	}
sens_max
sens_min

length(sens_max[19:64,])*2
sum(sens_max[19:64,]>=1)+sum(sens_min[19:64,]<=-1)

# only for 100% - number of parameters with no high sens
# sum(sens_max[19:64,]<=1 & sens_min[19:64,]>=-1)

###################### plot of the evolution of significant extremes #####################

# I calculated all max and mins in any uncertainty level and count the extremes that will be greater then 1 in absolut value.
# this means that each entry is 
x=c(29,39,49,136,318,736)


plot(x*100/max(x), xaxt = "n", xlab='Uncertainty',
	cex=1.5,
	main='Evolution of extreme sensitivities',
	ylab='% of extremes with abs value greater than one',
	col=c(rep('blue',5),'orange'),
	cex.axis=1.5,
	cex.lab=1.5,
	,pch=c(1,1,1,1,1,2))
axis(1, at=1:6, labels=c('5%','10%','20%','50%','100%','MAX'),cex.lab=1.5)





##############
# Figure S1b #
###################### plot of the evolution of standard deviations ################ 

	load("MSA_for_MK15_100p_clean_cut.RData")

MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_SD=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_SD)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_SD)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_SD[i,ii]=sqrt(var(x_temp))
		}
	}
sens_SD[19:64,]
sens_SD_vector_100=as.vector(sens_SD[19:64,])



	load("MSA_for_MK15_5p_clean.RData")		
	
MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_SD=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_SD)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_SD)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_SD[i,ii]=sqrt(var(x_temp))
		}
	}
sens_SD[19:64,]
sens_SD_vector_5=as.vector(sens_SD[19:64,])



	load("MSA_for_MK15_10p_clean.RData")	
	
MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_SD=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_SD)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_SD)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_SD[i,ii]=sqrt(var(x_temp))
		}
	}
sens_SD[19:64,]
sens_SD_vector_10=as.vector(sens_SD[19:64,])


	
	load("MSA_for_MK15_20p_clean.RData")		
	
MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_SD=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_SD)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_SD)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_SD[i,ii]=sqrt(var(x_temp))
		}
	}
sens_SD[19:64,]
sens_SD_vector_20=as.vector(sens_SD[19:64,])



	load("MSA_for_MK15_50p_clean.RData")		
	
MSA_raw_for_ext=MSA_raw_list

nSol=length(MSA_raw_for_ext)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_ext[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_ext[[1]])[2]							# number of variables

MSA_raw_for_ext[[1]]

sens_SD=matrix(rep(NA,nPar*nVar),ncol=nVar)
colnames(sens_SD)=colnames(MSA_raw_for_ext[[1]])
rownames(sens_SD)=rownames(MSA_raw_for_ext[[1]])
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		x_temp=vector()
		for (iii in 1:nSol)
			{
			x_temp=c(x_temp,MSA_raw_for_ext[[iii]][i,ii])
			}
		sens_SD[i,ii]=sqrt(var(x_temp))
		}
	}
sens_SD[19:64,]
sens_SD_vector_50=as.vector(sens_SD[19:64,])



boxplot(cbind(sens_SD_vector_5,sens_SD_vector_10,sens_SD_vector_20,sens_SD_vector_50,sens_SD_vector_100),
	range=1.5,
	ylim=c(0,2),
	xaxt = "n",
	main='Evolution of standard deviations of sensitivities',
	xlab='Uncertainty',
	ylab='Sensitivities Standard Deviation',
	cex.axis=1.5,
	cex.lab=1.5
	)
axis(1, at=1:5, labels=c('5%','10%','20%','50%','100%'),cex.axis=1.5)

max(sens_SD_vector_100)






##############
# Figure S1c #
################################# plot showing sens distribution alterations #########################

load("MSA_for_MK15_100p_clean_cut.RData")
MSA_raw_for_100=MSA_raw_list	

load("MSA_for_MK15_50p_clean.RData")
MSA_raw_for_50=MSA_raw_list
	
load("MSA_for_MK15_20p_clean.RData")	
MSA_raw_for_20=MSA_raw_list

load("MSA_for_MK15_10p_clean.RData")
MSA_raw_for_10=MSA_raw_list

load("MSA_for_MK15_5p_clean.RData")		
MSA_raw_for_5=MSA_raw_list

nSol=length(MSA_raw_for_100)								# number of solutions to be analysed
nPar=dim(MSA_raw_for_100[[1]])[1]							# number of parameters
nVar=dim(MSA_raw_for_100[[1]])[2]							# number of variables

MSA_raw_for_100[[1]]

sens_positive=matrix(rep(NA,5),ncol=5)
for (i in 1:nPar)
	{
	for (ii in 1:nVar)
		{
		temp_5=vector()		
		temp_10=vector()
		temp_20=vector()
		temp_50=vector()
		temp_100=vector()
		for (iii in 1:nSol)
			{
			temp_5=c(temp_5,MSA_raw_for_5[[iii]][i,ii])
			temp_10=c(temp_10,MSA_raw_for_10[[iii]][i,ii])
			temp_20=c(temp_20,MSA_raw_for_20[[iii]][i,ii])
			temp_50=c(temp_50,MSA_raw_for_50[[iii]][i,ii])
			temp_100=c(temp_100,MSA_raw_for_100[[iii]][i,ii])
			}
		sens_positive=rbind(sens_positive,c(sum(temp_5>0),sum(temp_10>0),sum(temp_20>0),sum(temp_50>0),sum(temp_100>0)))
		}
	}

sens_positive

sens_positive[sens_positive[,1]>(5000*.4)&sens_positive[,1]<(5000*.6),]
which(sens_positive[,1]>(5000*.4)&sens_positive[,1]<(5000*.6))

sens_positive[sens_positive[,4]>(5000*.4)&sens_positive[,4]<(5000*.6),]
which(sens_positive[,4]>(5000*.4)&sens_positive[,4]<(5000*.6))

MC=matrix(sens_positive[146:dim(sens_positive)[1],4],ncol=8)
rownames(MC)=rownames(MSA_raw_for_100[[1]])[19:length(rownames(MSA_raw_for_100[[1]]))]
colnames(MC)=colnames(MSA_raw_for_100[[1]])
MC
MC>(5000*.4)&MC<(5000*.6)



plot(1:5,sens_positive[146,]*100/5000,
	type='l',
	col='blue',
	ylim=c(0,100),
	xaxt = "n",
	main='Evolution of the percentage of positive sensitivities',
	xlab='Uncertainty',
	ylab='% of positive observations',
	cex.axis=1.5,
	cex.lab=1.5
	)
for (i in 146:nrow(sens_positive))
	{
	points(1:5,sens_positive[i,]*100/5000,type='l',col='blue')
	}
axis(1, at=1:5, labels=c('5%','10%','20%','50%','100%'))









