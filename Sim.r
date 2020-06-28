###############################################################
###   Simulation des données (cas d'inventaires réguliers)  ###
###############################################################

#Choix du nombre des variables explicatives
L<-3
#Choix du nombre de groupe 
K<-3
#Définition des paramètre du melange (pi_k et theta_k de chaque groupe)

	#Choix des parametres des composantes de la loi melange
	#dans une matrice K x (L+1)
theta<-NULL
theta<- t(replicate(K, rnorm(L+1)))
	#Choix des paramètres pi_k avec sum(Pi)=1
Pi<- diff(c(0,sort(runif(K-1)),1))
	
#Fonction logistique pour donner les probabilite par ligne
logist<-function(theta,x,groupe){
	1/(1+exp(-(theta[groupe,]%*%c(1,x))))
}

simul<-function(ligne){
	# Choix de la loi par tirage multinomiale Pi 
	# Realisation de la loi de Z
	
	V <- rmultinom(1,1,prob=Pi) 
	Z <- which(V==1)
	
	# Nombre d'arbre observé de la ligne
	n <- sample(10:30, 1)
	
	
	#Realisation de la loi de Y|Z
	#Choix de l'espece dans son groupe les espece de 
	#10Z-9 à 10Z sont dans le groupe Z: 1 à 10 pour le groupe 1
	
	Esp<- sample((10*Z-9):(10*Z), 1)
	
	#Valeurs de surfaces terrieres
	
	DvariableExpl<- rnorm(L)
	
	#Probabilité générer pour l'espèce choisi
	
	p<-logist(theta,DvariableExpl,Z)
	
	#Réalisation binomiale: nombre de grandi et nombre d'arbre
	ligne <- c(sum(rbinom(n,1,p)),n,Esp,DvariableExpl)
}

#Definition du tableau de donnees
N<-100
Tableau <- NULL
Tableau <- matrix(nrow=N,ncol=L+3)

#Simulation
Tableau<-t(apply(Tableau,1,simul))


##pour les courbes logistiques
#M<-Tableau[which(Tableau[,3]%in%21:30),]
#M<-t(apply(M,1,function(x) c(x, logist(EM_theta,x[-c(1:3)],2))))
#M<-t(apply(M,1,function(x) c(x, logist(theta,x[-c(c(1:3),7)],1))))
#png("courbe3.png")
#par(mfcol=c(2,1))
#plot(logit(M[,7]),M[,7], xlab="logit(p)", ylab="p", main="Courbe logistique estimée du groupe3")
#plot(logit(M[,8]),M[,8], xlab="logit(p)", ylab="p", main="Courbe logistique réelle du groupe3")
#dev.off()