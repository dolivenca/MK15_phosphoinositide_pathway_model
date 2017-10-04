

################################# the search for identifiabilities ##########################
#################################		 Voit's way 		 ###########################

rm(list=ls())	#clean memory
graphics.off() 	#close plot windows


# load local sensitivities
load("Local_sens_analysis_for_MK15_1.RData")					# load data
# Z is parm == columns, dep. var. == rows

dim(SA_raw)
Z=t(SA_raw[19:64,])						# taking proteins out. With the proteins the result is the same.
Z

# Changing the names in the matrix	       f5_PI4P f12_PI5P  f8_PI3P  f14_PI45P2 f17_PI35P2  f11_PI45P2 f31_PI4P    gamma_e
 colnames(Z)[c(11,25,17,29,33,23,42,38)] = c('f4->0','f5->45','f3->35','f45->345', 'f35->5',   'f45->4',  'f4->34',   '?i->')
# if all goes well this are the names of the identifiable parameters

euclidean_norm = function (X) {sqrt(sum(X^2))}				# creating a function to calculate de euclidean mean






e_n=vector()					# vector to store the euclidean norms
for (i in 1:dim(Z)[2])				# cycling the columns (parameters) of Z
	{
	e_n[i]=euclidean_norm(Z[,i])		# calculation the euclidean norm of each column (parameter)
	}
max_e_n=which(e_n==max(e_n))			# storing the position of the maximum euclidean norm
pos_v=max_e_n
max_e_n_v=max(e_n)

X_l=Z[,max_e_n]					# X_l vector containning the column with the maximum euclidean mean that we are working with
X_l_v=X_l						# X_l_v matrix containing all maximum euclidean norm vectors
X_l_v_names=colnames(Z)[max_e_n]		# names of the parameters that presented maximum euclidean nerm vectors


Z_l_hat=(X_l%*%(t(X_l)%*%X_l)^(-1))%*% t(X_l) %*% Z   # calculating Z_l_hat, projection of the X_l vector on the other columns of matrix Z or R_l
dim(Z_l_hat)

R_l = Z - Z_l_hat					# calculating resedues matrix, the matrix we get after taking the projection of X_l from Z or R_l-1 
#R_l[,max_e_n]=rep(0,8)

trashold=10^-200					# trashold that if a matrix presents a maximum below it, we consider all parameters represented in it as structurally non-identifiable. 
control=TRUE					# initiating while cycle control
while (control==TRUE)				# while the trashold is not crossed we repit the cycle
	{
	e_n=vector()				# vector to store the euclidean norms
	for (i in 1:dim(R_l)[2])		# cycling the columns (parameters) of Z
		{
		e_n[i]=euclidean_norm(R_l[,i])# calculation the euclidean norm of each column (parameter)
		}
	max_e_n=which(e_n==max(e_n))		# storing the position of the maximum euclidean norm
	pos_v=c(pos_v,max_e_n)
	max_e_n_v=c(max_e_n_v,max(e_n))	

	if(max(e_n)>trashold)			# if the trashold is not crossed
		{
		X_l = R_l[,max_e_n]		# X_l vector containning the column with the maximum euclidean mean that we are working with
		X_l_v=cbind(X_l_v,R_l[,max_e_n])	# X_l_v matrix containing all maximum euclidean norm vectors
		X_l_v_names=c(X_l_v_names,colnames(Z)[max_e_n])		# names of the parameters that presented maximum euclidean nerm vectors
			
		Z_l_hat=(X_l%*%(t(X_l)%*%X_l)^(-1))%*% t(X_l) %*% R_l		# calculating Z_l_hat, projection of the X_l vector on the other columns of matrix Z or R_l
		
		R_l =  R_l - Z_l_hat		# calculating resedues matrix, the matrix we get after taking the projection of X_l from Z or R_l-1 
		for (ii in pos_v)
			{
			R_l[,ii]=rep(0,8)	
			}
		} else {control=FALSE}		# if trashold crossed, stop while cycle
	}

colnames(X_l_v)=X_l_v_names			# put the names of the parameters in the column vectors used
X_l_v							# the identifiable parameters 
dim(X_l_v)

#           f5_PI4P  f8_PI3P  f12_PI5P f14_PI45P2 f17_PI35P2 f31_PI4P  f11_PI45P2  gamma_e
# sens_names=c('f5->45','f3->35','f45->345','f4->0','f35->5', 'f4->34','f45->4','?i->',rep('',39))


##############
# Figure S1d #
##############
plot(1:length(max_e_n_v),max_e_n_v,col='blue',cex=1.5,cex.axis=1.5,cex.lab=1.5,
	main='Parameter identifiability',
	xlab='Parameters',
	ylab='Euclidean norm of the parameter sensitivities')
abline(h=.3,col='cyan')
text(1:length(max_e_n_v),max_e_n_v, labels=X_l_v_names, cex=1.2,pos=c(4,4,3,4,1,4,4,4,rep(1,39)))
#sens_names

#X_l_v_names[1:length(max_e_n_v)]
eigen(Z)


