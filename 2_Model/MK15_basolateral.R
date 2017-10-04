#############################################
###							###
###	Phosphoinositide pathway Model	###
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

str(out)
varnames=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
varcolors=c('gray','darkgoldenrod1','coral','coral1','coral2','orange','orange2','darkorange3','gold','deepskyblue1','deepskyblue3')

# naming fluxes in out
flux_names=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V21','V22','V23','V24','V25','V26','V27','V28','V29','V30','V31','V32','V33')
colnames(out)[10:dim(out)[2]]=flux_names

par(mfrow=c(3,1))
# time plot
plot(out[,1],out[,2],type='l',col=1,ylim=c(0,max(out[,2:9])),main='Time course',xlab='Time (m)',ylab='Molecule count')
for (i in 2:9) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,15000),xlab='Time (m)',ylab='Molecule count')
for (i in 2:9) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}

plot(out[,1],out[,2],type='l',col=1,ylim=c(0,1500),xlab='Time (m)',ylab='Molecule count')
for (i in 2:9) 
	{
	points(out[,1],out[,i],type='l',col=varcolors[i])
	text(tmax,out[tmax-1,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))
           
# list (primitive)
(head(out,10))[,1:9]
(tail(out,10))[,1:9]

sum(out[nrow(out),2:9])

# edit(out)
# View(out)

# turnover = concentration / flux		high turnover = small concentration / big flux
(PIturnover=tail(out,1)[2]/(tail(out,1)[11]+tail(out,1)[13]+tail(out,1)[15]+tail(out,1)[31]))
(PI3Pturnover=tail(out,1)[3]/(tail(out,1)[12]+tail(out,1)[17]+tail(out,1)[34]))
(PI4Pturnover=tail(out,1)[4]/(tail(out,1)[14]+tail(out,1)[19]+tail(out,1)[32]+tail(out,1)[38])) 
(PI5Pturnover=tail(out,1)[5]/(tail(out,1)[16]+tail(out,1)[21]+tail(out,1)[36])) 
(PI35P2turnover=tail(out,1)[6]/(tail(out,1)[17]))
(PI45P2turnover=tail(out,1)[7]/(tail(out,1)[20]+tail(out,1)[22]+tail(out,1)[23]+tail(out,1)[30]))
(PI34P2turnover=tail(out,1)[8]/(tail(out,1)[29]+tail(out,1)[38]))
(PI345P3turnover=tail(out,1)[9]/(tail(out,1)[23]))
turnover=c(PIturnover,PI3Pturnover,PI4Pturnover,PI5Pturnover,PI35P2turnover,PI45P2turnover,PI34P2turnover,PI345P3turnover) 
names(turnover)=names((tail(out,1))[,2:9])
turnover

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()
par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)





############################################## MAP ################################################

# run this code lines if you do not have the packges installed
# install.packages('igraph')

library(igraph) 	   #library for graphs
# citation("igraph")


# readkey function. To stop r until enter is pressed. Not used any more.
#readkey <- function()
#{
#    cat ("Press [enter] to continue")
#    line <- readline()
#}

##### use to create the matrix #####
# names of the colunms and rows of the matrix
#node_names=c('PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3','PI_aux','PI3P_aux','PI4P_aux','PI5P_aux','PI35P2_aux','PI45P2_aux','PI34P2_aux','PI345P3_aux','Time','Time_aux')
#m=matrix(0,nrow=18,ncol=18)							# creating the matrix
#dimnames(m)=list(node_names,node_names) 					# giving the names to the matrix cols and rows
#m
#m=edit(m)
#write.table(m, file="map_matrix_MK6_perfect.txt", row.names=TRUE, col.names=TRUE) #save the matrix in a file
(m=as.matrix(read.table("map_matrix_MK6_perfect.txt")))  #read the map_matrix. Must be on correct directory 

##### use to create the node coordenats matrix #####
#vertex_coords=matrix(0,nrow=16,ncol=2)
#vertex_coords=edit(vertex_coords)
#write.table(vertex_coords, file="vertex_coords_MK6_perfect.txt", row.names=TRUE, col.names=TRUE) #save the node coordenates matrix in a file
(vertex_coords=as.matrix(read.table("vertex_coords_MK6_perfect.txt")))  #read the node coordenates matrix. Must be on correct directory 
cw=600	# canvas width
ch=800	# canvas leight


g=graph.adjacency(m,mode="directed",weighted=NULL,diag=FALSE) 	# creating igraph graphs from adjacency matrices

g$layout = vertex_coords							# coordenats for the nodes

# settings for the vertexs
#V(g)$shape='rectangle'
V(g)$color="white"
V(g)$label=paste(sep="\n",V(g)$name,degree(g))				
V(g)$label.color=c(rep("blue",8),rep('white',8),'grey','white')


# settings for the edges
E(g)$color=c('red','red','red','red','darkgray','blue','red','darkgray','blue','red','red','darkgray','blue','red','darkgray','blue','blue','darkgray','blue','blue','blue','red','darkgray','blue','blue','darkgray','blue','blue','darkgray','darkgray','darkgray','darkgray','darkgray')
E(g)$label=c('V2','V4','V6','V30','V23','V3','V8','V26','V5','V10','V31','V24','V7','V12','V28','V9','V17','V27','V33','V11','V13','V14','V22','V19','V32','V29','V15','V21','V25','V1','V18','V16','ref')
E(g)$label.color='black'
E(g)$arrow.size=.5
E(g)$curved=c(.5,.5,.5,-3,1.5,.5,.5,1.5,.5,.5,.5,1.5,.5,.5,0,.5,0,0,3.5,.5,.5,.5,0,0,0.5,0,.5,0,0,1.5,1.5,1.5,0)
cbind(E(g)$label,E(g)$color,E(g)$curved)

list.edge.attributes(g)

# put the normalized fluxes in the graph order
namesV=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V21','V22','V23','V24','V25','V26','V27','V28','V29','V30','V31','V32','V33')   # names of the fluxes in out_norm order 
ordered_norm_fluxes=out_norm[,10:ncol(out_norm)]	# creating matrix to store the graph ordered normalized fluxes
for (i in 1:length(E(g)$label))				# cycle the names of the graph edges
	{
	for (ii in 1:ncol(out_norm))				# cycle the names of the out_norm colunms
		{
 		if (namesV[ii]==E(g)$label[i]) {ordered_norm_fluxes[,i]=out_norm[,ii]} 	# if the names are equal store the colunm of outnorm in the matrix
		}
	}

# put the actual fluxes in the graph order
namesV=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3','V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V21','V22','V23','V24','V25','V26','V27','V28','V29','V30','V31','V32','V33')   # names of the fluxes in out_norm order 
ordered_fluxes=out[,10:ncol(out)]				# creating matrix to store the graph ordered actual fluxes
for (i in 1:length(E(g)$label))				# cycle the names of the graph edges
	{
	for (ii in 1:ncol(out))					# cycle the names of the out colunms
		{
 		if (namesV[ii]==E(g)$label[i]) {ordered_fluxes[,i]=out[,ii]} 	# if the names are equal store the colunm of out in the matrix
		}
	}



# simple graph
i=4900
windows()
plot.igraph(g,
	vertex.label=paste(sep="\n",V(g)$name,c(round(out[i,2:9],1),rep(1,8),out[i,1],1)),
	vertex.size=c(rep(30,8),rep(1,8),35,1),
	vertex.label.cex=c(1,1,1,1,1,1,1,1,.1,.1,.1,.1,.1,.1,.1,.1,1,.1),
	edge.label=paste(sep="\n",E(g)$label,c(round(ordered_fluxes[i,1:ncol(ordered_norm_fluxes)],3))),
	edge.label.cex=.7,
	edge.arrow.size=cbind(ordered_norm_fluxes[i,1:ncol(ordered_norm_fluxes)],1),
	edge.width=c(5*ordered_norm_fluxes[i,1:ncol(ordered_fluxes)],5*1)
	)
















#################################   Perturbations   ###################################
perturbations=cbind(c('NA','NA'),c('NA','NA'),c('NA','NA'))	# initiation of perturbations matrix, if unaltered it will do no perturbations
# Define perturbations. 
# The first row of the matrix is allways 0,novar,NA for the cycle use the inicial values at the start
# The last row of the matrix is always max(t),novar,NA for the cycle to work.
# The values in the middle rows are the perturbations to the system of diferential equations.



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



##############
# Figure 3 a #
##############
### PI45P2 sensitive to PI, pi_4K, pi_4K_plus_pip_5KI.
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,2000,3000,3005,4000,4005,5000,5005,6000,6005,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','gamma1','gamma1','pi_4K','pi_4K_pip_5KI','pi_4K','pi_4K_pip_5KI','pip_5KI','pi_4K_pip_5KI','pip_5KI','pi_4K_pip_5KI','novar')
# On the third the new value for the variable.
pval=c(NA,Pars$gamma1*.5,Pars$gamma1,Pars$pi_4K*0.05,Pars$pi_4K_pip_5KI*0.05,Pars$pi_4K,Pars$pi_4K_pip_5KI,Pars$pip_5KI*.5,Pars$pi_4K_pip_5KI*.5,Pars$pip_5KI,Pars$pi_4K_pip_5KI,NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out_test)

str(out)
varnames=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
varcolors=c('gray','darkgoldenrod1','coral','coral1','coral2','orange','orange2','darkorange3','gold','deepskyblue1','deepskyblue3')

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()
par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

### bar plots for P1 ###

# PI45P2 and PI4P plots
# for perturbation 3
# Figure 3

# edit(out_norm)
sample1_PI45P2=out_norm[c(500,1500,3500,5500),7]
sample1_PI4P=out_norm[c(500,1500,3500,5500),4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

windows()
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.2),
	ylab="Percentage of steady state levels",
	cex.axis=1,cex.lab=1.5
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 1.3)
legend("top", c("PI(4,5)P2","PI(4)P"), cex=1.5, bty="n",fill=c("black","grey"))
points(x=c(4.5,5.5,7.5,8.5,10.5,11.5),y=c(rep(.5,5),1),col='blue', pch = "_", cex = 7)






##############
# figure 3 b #
##############
### pi_4K, pi_4K_plus_pip_5KI, and PI5P can mantain PI45P2
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,1000,1500,2000,2500,3000,4000,5000,6000,7000,8000,9000,9500,max(t));
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','gamma16','pi_4K_pip_5KI','gamma12','pi_4K','pi_4K','pi_4K','pi_4K_pip_5KI','pi_4K_pip_5KI','gamma16','gamma16','gamma12','gamma12','novar')
# On the third the new value for the variable.
pval=c(NA,0,0,0,Pars$pi_4K*0.2,Pars$pi_4K,Pars$pi_4K*0.2,Pars$pi_4K_pip_5KI,0,Pars$gamma16,0,Pars$gamma12,0,NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out_test)

str(out)
varnames=c('Time','PI','PI3P','PI4P','PI5P','PI35P2','PI45P2','PI34P2','PI345P3')
varcolors=c('gray','darkgoldenrod1','coral','coral1','coral2','orange','orange2','darkorange3','gold','deepskyblue1','deepskyblue3')

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()
par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

### bar plots for P1 ###

# PI45P2 and PI4P plots 
# for perturbation 1 
# (Figure 2)

# edit(out_norm)
time_sample1=c(0,2500,3000,5000,7000,9000)
sample1_PI45P2=out_norm[time_sample1+100,7]
sample1_PI4P=out_norm[time_sample1+100,4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

windows()
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.05),
	ylab="Percentage of steady state levels",
	cex.axis=1,cex.lab=1.5	
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 1.3)
legend("topright", c("PI(4,5)P2","PI(4)P"), cex=1, bty="n",fill=c("black","grey"))


################
# (Figure 3 d) #
################
# edit(out_norm)
time_sample1=c(0,1000,2500,7000)
sample1_PI45P2=out_norm[time_sample1+100,7]
sample1_PI4P=out_norm[time_sample1+100,4]
joint=rbind(sample1_PI45P2,sample1_PI4P)

windows()
bp=barplot(joint,beside = TRUE,ylim=c(0,max(joint)+.05),
	ylab="Percentage of steady state levels",
	cex.axis=2,cex.lab=1.5	
	)
bp
text(bp,joint, labels = format(round(joint,2), 4),pos = 3, cex = 2)
legend("topright", c("PI(4,5)P2","PI(4)P"), cex=1, bty="n",fill=c("black","grey"))






################
# (Figure 3 c) #
################
### knock down of MTMR (myotubularyns) to 70% should make PI5P 20% and PI35P2 150% 
### knock out of SYNJ_TMEM55 should decrease 15%
### Playing with PIKfyve. 50% of PIkfyve should put PI5P and PI35P2 between 35% and 39% of original steady state levels
# Finaly remove the # form the front of all the next code lines
# On the first column put the time of the perturbation
ptime=c(0,501,1000,1500,2000,2500,3000,3500,4000,5000,6000,max(t))
# On the second column put the variables to be altered. (ex: 'X')
pvar=c('novar','MTMR','MTMR','SYNJ_TMEM55','SYNJ_TMEM55','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','pi_Kfyve','novar')
# On the third the new value for the variable.
pval=c(NA,MTMR*.65,MTMR,SYNJ_TMEM55*0,SYNJ_TMEM55,pi_Kfyve*.1,pi_Kfyve,pi_Kfyve*.001,pi_Kfyve,pi_Kfyve*.5,pi_Kfyve,NA)
perturbations=data.frame(time=ptime,var=pvar,val=pval)
perturbations

out = Cruncher(state,t,equations,parameters,perturbations)						# no perturbations			# with perturbations
# edit(out)

# normalized graphs
stst=out[500,]						# getting steady state at time 500
out_n=cbind(time=out[1],out[,2:ncol(out)])
tail(out_n)
out_norm=out[,1]						# creating matrix out_norm to store the normalized information
for (i in 2:ncol(out))
	{
	newc=c(out_n[,i]/as.numeric(stst[i]))	# normalizing to steady state got at time 500
	out_norm=cbind(out_norm,newc)			# storing normalized information
	}
colnames(out_norm)=colnames(out)
head(out_norm)
# edit(out_norm)

windows()

par(mfrow=c(3,3))
for (i in 2:9) 
	{
	plot(out_norm[,1],out_norm[,i],type='l',col=varcolors[i],xlab='Time',ylab='Fold change',main=varnames[i])
	text(tmax,out_norm[tmax,i],varnames[i],col=varcolors[i])
	}
par(mfrow=c(1,1))

# edit(out_norm)

#       check_pert check_names check_targets check_darts
#1        MTMR 10%        PI5P         0.200 0.300001244
#2        MTMR 10%      PI35P2         1.500 1.744820658
#3            <NA>        <NA>            NA          NA
#4  SYNJ_TMEM55 0%        PI5P         0.850 0.708696492
#5            <NA>        <NA>            NA          NA
#6    pi_Kfyve 10%        PI5P         0.500 0.463738581
#7    pi_Kfyve 10%      PI35P2         0.500 0.377789832
#8            <NA>        <NA>            NA          NA
#9    pi_Kfyve .1%        PI5P         0.150 0.207333684
#10   pi_Kfyve .1%      PI35P2         0.001 0.005400449
#11   pi_Kfyve .1%      PI45P2         0.800 0.901334972
#12   pi_Kfyve .1%        PI3P         5.000 5.467992807

# check for small PIs
check_pert = c('MTMR 10%','MTMR 10%',NA,'SYNJ_TMEM55 0%',NA,'pi_Kfyve 10%','pi_Kfyve 10%',NA,'pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%','pi_Kfyve .1%',NA,'pi_Kfyve 50%','pi_Kfyve 50%')   	
check_names = c('PI5P','PI35P2',NA,'PI5P',NA,'PI5P','PI35P2',NA,'PI5P','PI35P2','PI45P2','PI3P',NA,'PI5P','PI35P2')
check_targets = c(.2,1.5,NA,.85,NA,.5,.5,NA,.15,.001,.80,5,NA,35,39)
check_darts = c(out_norm[900,5],out_norm[900,6],NA,out_norm[1900,5],NA,out_norm[2900,5],out_norm[2900,6],NA,out_norm[3900,5],out_norm[3900,6],out_norm[3900,7],out_norm[3900,3],NA,out_norm[4900,5],out_norm[4900,6])
names(check_darts)=NULL
check = data.frame(check_pert,check_names,check_targets,check_darts) 
check

perturbations

data=c(out_norm[900,5],out_norm[900,6],NA,out_norm[1900,5],NA,out_norm[2900,5],out_norm[2900,6],NA,out_norm[3900,5],out_norm[3900,6],out_norm[3900,7],out_norm[3900,3])	
windows()
bp=barplot(data,
		col=c('gray92','gray','blue','gray92','blue','gray92','gray','blue','gray92','gray','dimgray','black'),
		ylab="Percentage of steady state levels",ylim=c(0,max(data,na.rm = TRUE)+2),
		cex.lab=1.5
		)
bp
text(bp,data, labels = format(round(data,2), 4),pos = 3, cex = 1.5)
legend("topleft", c("PI(5)P","PI(3,5)P2","PI(4,5)P2","PI(3)"), cex=1, bty="n",fill=c('gray92','gray','dimgray','black'))
points(x=bp,
	y=c(.2,1.5,-1,.85,-1,.5,.5,-1,.15,.001,.80,5),
	col='blue', pch = "_", cex = 5)




