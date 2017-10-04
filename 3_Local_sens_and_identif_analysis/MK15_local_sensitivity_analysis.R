### Sensitivity analysis MK1 ###

# run this code lines if you do not have the packges installed
# install.packages('deSolve')

library(deSolve)	#library for solve diferential equations
# citation("deSolve")



										

#################################################################################
###												    	###
### run a model first 										### 
### you should have in memory the t, state, parameters, equations variables	###				
### Run the model to get a time course storede in the variable out.		###													
### Please run the basolateral setup.							###
###													###
#################################################################################

ls()	# see vars in memory
t
state
parameters
equations
out



Ss_old=tail(out,1)[,2:9]							# storing steady state to be studied 
Sens=matrix(rep(NA,length(parameters)*length(Ss_old)),ncol=8)	# creating matrix to store the sensitivities
colnames(Sens)=names(Ss_old)							# naming the columns after the dependent variables
rownames(Sens)=names(parameters)						# naming the rows after the parameters
Sens

t=seq(0,5000+1,1) 	# time

for (i in 1:length(parameters))				# cycling the parameters
	{
	Pars=parameters						# creating discardable parameter vector to be changed
	Pars[i]=parameters[i]*(1.01)				# changing one parameter in the discardable parameter vector in 1%
	out_new <- ode(y = state, times = t, func = equations,parms=Pars)		#runing the ode solver
	Ss_new=tail(out_new,1)[,2:9]				# new steady state
	tsens=vector()						# creating sensitivity vector for this turn 
	for (ii in 1:8)						# cycling all dependent variables
		{
		tsens=c(tsens,(Ss_new[ii]-Ss_old[ii])*100/Ss_old[ii])	#calculating sensitivities
		}
	Sens[i,]=tsens						# puting this parameter sensitivities in sentitivity matrix
	}
Sens									# show all sensitivity matrix
sum(abs(Sens))							# sum of all sensitivities absolute value
sum(abs(Sens)>=1)							# amount of significant sensibilities (greater then 1%)

Sens_show=abs(Sens)>=1						# mark significant sensibilities 
for(iii in 1:length(parameters))				# cycling the rows of the sensitivity matrix
	{
	for(iv in 1:8)						# cycling the columns of the sensitivity matrix
		{
		if (Sens_show[iii,iv]==1){Sens_show[iii,iv]=Sens[iii,iv]}	# substitute the TRUEs for the significant sensibilities
		}
	}
Sens_show								# show sensitivities matrix with just the significant sensibilities
sum(abs(Sens_show[Sens_show!=0]))				# sum of the significant sensitivities absolute value

dim(Sens_show)							# dimentions of sensitivity matrix
dim(Sens_show)[1]*dim(Sens_show)[2]				# total sensitivities accessed 
sum(abs(Sens)>=1)/(dim(Sens_show)[1]*dim(Sens_show)[2]) # percentage of high sensitivities





Pars-parameters
Ss_new-Ss_old



# save and read from disk	
# as a table
# write.table(Sens_show,file="6sensMK1")			# write to disk
# sens_MK1=read.table(file="6sensMK1")			# read from disk


SA_raw=Sens
SA_simple=Sens_show
# save and read from disk
# as a variable in memory
# save(SA_raw,SA_simple,file = "Local_sens_analysis_for_MK15_1.RData")	
# load( "Local_sens_analysis_for_MK15_1.RData")



























################################ parameters 10% more ############################


### run a model first ### 

parameters_110=parameters*1.1

out_ref=ode(y = state, times = t, func = equations,parms=parameters_110)
Ss_old_110=tail(out_ref,1)[,2:9]								# storing steady state to be studied 
Sens_110=matrix(rep(NA,length(parameters_110)*length(Ss_old_110)),ncol=8)	# creating response matrix
colnames(Sens_110)=names(Ss_old_110)							# naming the columns after the dependent variables
rownames(Sens_110)=names(parameters_110)							# naming the rows after the parameters
Sens_110

t=seq(0,5000+1,1) 	# time

for (i in 1:length(parameters_110))				# cycling the parameters
	{
	Pars=parameters_110					# creating discardable parameter vector to be changed
	Pars[i]=parameters_110[i]*(1.01)			# changing one parameter in the discardable parameter vector in 1%
	out_new <- ode(y = state, times = t, func = equations,parms=Pars)		#runing the ode solver
	Ss_new=tail(out_new,1)[,2:9]				# new steady state
	tsens=vector()						# creating sensitivity vector for this turn 
	for (ii in 1:8)						# cycling all dependent variables
		{
		tsens=c(tsens,(Ss_new[ii]-Ss_old_110[ii])*100/Ss_old_110[ii])	#calculating sensitivities
		}
	Sens_110[i,]=tsens					# puting this parameter sensitivities in sentitivity matrix
	}
Sens_110					 			# show all sensitivity matrix
sum(abs(Sens_110))						# sum of all sensitivities absolute value
sum(abs(Sens_110)>=1)						# amount of significant sensibilities (greater then 1%)

Sens_show_110=abs(Sens_110)>=1				# mark significant sensibilities 
for(iii in 1:length(parameters_110))			# cycling the rows of the sensitivity matrix
	{
	for(iv in 1:8)						# cycling the columns of the sensitivity matrix
		{
		if (Sens_show_110[iii,iv]==1){Sens_show_110[iii,iv]=Sens_110[iii,iv]}	# substitute the TRUEs for the significant sensibilities
		}
	}
Sens_show_110							# show sensitivities matrix with just the significant sensibilities
sum(abs(Sens_show_110[Sens_show_110!=0]))			# sum of the significant sensitivities absolute value

dim(Sens_show_110)						# dimentions of sensitivity matrix
dim(Sens_show_110)[1]*dim(Sens_show_110)[2]				# total sensitivities accessed 
sum(abs(Sens_110)>=1)/(dim(Sens_show_110)[1]*dim(Sens_show_110)[2]) # percentage of high sensitivities








################################ parameters 10% less ############################


### run a model first ### 

parameters_90=parameters*.9

out_ref=ode(y = state, times = t, func = equations,parms=parameters_90)
Ss_old_90=tail(out_ref,1)[,2:9]								# storing steady state to be studied 
Sens_90=matrix(rep(NA,length(parameters_90)*length(Ss_old_90)),ncol=8)		# creating response matrix
colnames(Sens_90)=names(Ss_old_90)								# naming the columns after the dependent variables
rownames(Sens_90)=names(parameters_90)							# naming the rows after the parameters
Sens_90

t=seq(0,5000+1,1) 	# time

for (i in 1:length(parameters_90))				# cycling the parameters
	{
	Pars=parameters_90					# creating discardable parameter vector to be changed
	Pars[i]=parameters_90[i]*(1.01)			# changing one parameter in the discardable parameter vector in 1%
	out_new <- ode(y = state, times = t, func = equations,parms=Pars)		#runing the ode solver
	Ss_new=tail(out_new,1)[,2:9]				# new steady state
	tsens=vector()						# creating sensitivity vector for this turn 
	for (ii in 1:8)						# cycling all dependent variables
		{
		tsens=c(tsens,(Ss_new[ii]-Ss_old_90[ii])*100/Ss_old_90[ii])	#calculating sensitivities
		}
	Sens_90[i,]=tsens						# puting this parameter sensitivities in sentitivity matrix
	}
Sens_90					 			# show all sensitivity matrix
sum(abs(Sens_90))							# sum of all sensitivities absolute value
sum(abs(Sens_90)>=1)						# amount of significant sensibilities (greater then 1%)

Sens_show_90=abs(Sens_90)>=1					# mark significant sensibilities 
for(iii in 1:length(parameters_90))				# cycling the rows of the sensitivity matrix
	{
	for(iv in 1:8)						# cycling the columns of the sensitivity matrix
		{
		if (Sens_show_90[iii,iv]==1){Sens_show_90[iii,iv]=Sens_90[iii,iv]}	# substitute the TRUEs for the significant sensibilities
		}
	}
Sens_show_90							# show sensitivities matrix with just the significant sensibilities
sum(abs(Sens_show_90[Sens_show_90!=0]))			# sum of the significant sensitivities absolute value

dim(Sens_show_90)							# dimentions of sensitivity matrix
dim(Sens_show_90)[1]*dim(Sens_show_90)[2]			# total sensitivities accessed 
sum(abs(Sens_90)>=1)/(dim(Sens_show_90)[1]*dim(Sens_show_90)[2]) # percentage of high sensitivities








######################## analysis of the three matrixes #########################

Sens_3=matrix(rep(NA,length(parameters)*length(Ss_old)),ncol=8)			# creating response matrix
colnames(Sens_3)=names(Ss_old)								# naming the columns after the dependent variables
rownames(Sens_3)=names(parameters)								# naming the rows after the parameters
Sens_3


for (i in 1:length(parameters))								# cycling lines
	{
	for (ii in 1:length(Ss_old))								# cycling columns
		{
		tempvar1=0										# initianing local variable tempvar1
		if (abs(Sens[i,ii])>=1) {tempvar1=100}					# if the sensitivity is significant in Sens mark 1
		if (abs(Sens_110[i,ii])>=1) {tempvar1=tempvar1+20}			# if the sensitivity is significant in Sens_110 mark 2
		if (abs(Sens_90[i,ii])>=1) {tempvar1=tempvar1+3}			# if the sensitivity is significant in Sens_90 mark 3
		Sens_3[i,ii]=tempvar1								# save the information in the matrix Sens_3
		}
	}

Sens_3

