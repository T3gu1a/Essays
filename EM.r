 #####################################################################
 ###    Estimation d'un melange de regression logistique par EM    ###
 #####################################################################
 
 
 #fonction utilisé dans l'optimisation des theta pour le calcul des réels qui 
 #pondèrent les vecteurs lignes de données dans l'expression du gradient
 fvecteur<-function(x,u){
			(x[1]-x[2]/(1+exp(-x[-c(1,2)]%*%u)))*x[-c(1,2)]
		}
 #fonction a optimiser pour trouver les theta_k
 foptimtheta<-function(u,k){
			norme(colSums(probaAposteriori[,k]*t(sapply(listDparEspece1vaExpl,
			     function(M) colSums(t(apply(M,1,function(x) fvecteur(x,u))))))))
	}
 #calcul des scalaires par espece
 fMatriceScalaire<-function(x){
		t(apply(x,1,function(y) c(y[1],y[2],EM_theta%*%c(1,y[4:(L+3)]))))
		}
 #calcul des vraisemblances par espece
 fVraisParEspece<-function(x){
		 EM_Pi*apply(t(apply(x,1, function(y) dbinom(y[1],y[2],prob=1/(1+exp(-y[3:(K+2)]))))),2,prod)
	}

 #norme
 norme<- function(x) sqrt(t(x)%*%x)

 
 #Initialisation des parametres
	EM_theta<- t(replicate(K, rnorm(L+1)))
	EM_Pi<- diff(c(0,sort(runif(K-1)),1))
 
 #Séparation des donnees par espece
 listDparEspece <- lapply(split(Tableau, Tableau[,3]),function(x) matrix(x,ncol=L+3))
 
 	#liste des matrices des produits scalaires initiaux par espece
	listMatriceScalaire<-lapply(listDparEspece, fMatriceScalaire)
		 
	#liste des vraisemblances initiales par espece
	listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
	#Initialisation de la fonction a maximiser
	NQ<- 0

	repeat{
		 
		 #l'ancienne valeur est mise a jour
		 AQ<-NQ
		 
		 ###############
		 ### Etape E ###
		 #########################################################################

			 
			 #calcul des probabilites a posteriori
			 
			 probaAposteriori<-t(sapply(listVraisParEspece, function(x) x/sum(x)))
		 
		 ###############
		 ### Etape M ###  
		 #########################################################################
		 
			 #mise a jour des proportions pi_k
			 EM_Pi<-colMeans(probaAposteriori)
			 
			 #mise a jour des theta_k
			 
			 #Mise en forme des vecteurs de donnees des especes avec 1 comme premiere composante
			 listDparEspece1vaExpl<-lapply(listDparEspece, function(x) t(apply(x,1, function(y) c(y[c(1:2)],1,y[4:(L+3)]))))
			 
			 #optimisation des theta_k
			 EM_theta<-cbind(1:K, EM_theta)
			 EM_theta<-t(apply(EM_theta,1,function(A) optim(A[-1], foptimtheta, k=A[1])$par))
			 
	    
		
		### calcul de la nouvelle valeur de la fonction a maximiser : Q(theta|theta(m)) 
						#liste des matrices de produits scalaires par espece
						listMatriceScalaire<-lapply(listDparEspece, fMatriceScalaire)
		 
						#liste des vraisemblances par espece
						listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
						
						NQ<-sum(rowSums(probaAposteriori*log(t(sapply(listVraisParEspece, function(x) x)))))
		### Condition d'arret
		
			if(abs(NQ-AQ)<0.0000001) break
	}
 
