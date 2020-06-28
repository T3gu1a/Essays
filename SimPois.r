###############################################################
###   Simulation des donn�es (cas d'inventaires r�guliers)  ###
###############################################################

#Choix du nombre des variables explicatives
L<-3
#Choix du nombre de groupe 
K<-3
#D�finition des param�tre du melange (pi_k et theta_k de chaque groupe)

	#Choix des parametres des composantes de la loi melange
	#dans une matrice K x (L+1)
theta<-NULL
theta<- t(replicate(K, rnorm(L+1)))
	#Choix des param�tres pi_k avec sum(Pi)=1
Pi<- diff(c(0,sort(runif(K-1)),1))
	
#Fonction logistique pour donner les probabilite par ligne
poiss<-function(theta,x,groupe){
	exp((theta[groupe,]%*%c(1,x)))
}
simul_P<-function(ligne){
	# Choix de la loi par tirage multinomiale Pi 
	# Realisation de la loi de Z
	
	V <- rmultinom(1,1,prob=Pi) 
	Z <- which(V==1)
	
	#Realisation de la loi de Y|Z
	#Choix de l'espece dans son groupe les espece de 
	#10Z-9 � 10Z sont dans le groupe Z: 1 � 10 pour le groupe 1
	
	Esp<- sample((10*Z-9):(10*Z), 1)
	
	#Valeurs de surfaces terrieres
	
	DvariableExpl<- rnorm(L)
	
	#Probabilit� g�n�rer pour l'esp�ce choisi
	
	lamda<-poiss(theta,DvariableExpl,Z)
	
	#R�alisation binomiale: nombre de grandi et nombre d'arbre
	ligne <- c(rpois(1,lamda),Esp,DvariableExpl)
}

#Definition du tableau de donnees
N<-100
Tableau3 <- NULL
Tableau3 <- matrix(nrow=N,ncol=L+3)

#Simulation
Tableau3<-t(apply(Tableau3,1,simul_P))
