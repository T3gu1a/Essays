 

 EM_MortaliteCroissance<-function(Tableau, K){
 
 #####################################################################
 ###    Estimation d'un melange de regression logistique par EM    ###
 #####################################################################
 
 #nombre de variable explicative
 L<-ncol(Tableau)-3
 
 #fonction utilisee dans l'optimisation des theta pour le calcul des réels qui 
 #pondèrent les vecteurs lignes de données dans l'expression du gradient
 fvecteur<-function(x,u){
			 v<- x[-c(1,2)]
			(-x[1]+x[2]/(1+exp(-v%*%u)))*v
		}
 #fonction utilisee dans l'optimisation des theta pour le calcul des réels qui 
 #pondèrent les vecteurs lignes de données dans l'expression de la hessienne		
 fvecteur2<- function(x,u){
			 v <- x[-c(1,2)]		
		     ((x[2]/(2+2*cosh(v%*%u)))* v) %*% t(v)
	}
 #calcul du nabla en theta_k (u)
 nablaJ<-function(u,k){
			colSums(probaAposteriori[,k]*t(sapply(listDparEspece1vaExpl,
			     function(M) colSums(t(apply(M,1,function(x) fvecteur(x,u)))))))
	}
#calcul de la hessien en theta_k (u) 
 hessienJ<-function(u,k){			
		l<-lapply(listDparEspece1vaExpl, function(M) apply(M,1,function(x) fvecteur2(x,u)))
		l<-lapply(l, function(M) rowSums(M))
		Mat<-sapply(l, function(x) x)
		Mat<- rowSums(t(probaAposteriori[,k]*t(Mat)))
	    matrix(Mat, nrow=L+1, ncol=L+1)
	}
#algorithme de Newton-Raphson
 NewtonRaph<-function(MPoint_zero){
	Point_plus<-MPoint_zero
	repeat{
		Point_moins<-Point_plus
		chaq_theta<-split(Point_plus,Point_plus[,1])
		Hessien<-lapply(chaq_theta,function(y) hessienJ(y[-1],y[1]))
		Nabla<-lapply(chaq_theta,function(y) nablaJ(y[-1],y[1]))
		
		if(norme(sapply(Nabla, function(x) x))<0.0001) break
		
		Point_plus<-cbind(1:K,t(apply(Point_plus,1, function(x) Point_moins[x[1],][-1]- solve(Hessien[[x[1]]],Nabla[[x[1]]],tol=1e-22))))

	}
	Point_plus[,-1]
 }

  #calcul des scalaires par espece
 fMatriceScalaire<-function(x,M){
		t(apply(x,1,function(y) c(y[1],y[2],M%*%c(1,y[4:(L+3)]))))
		}
 #calcul des vraisemblances par espece
 fVraisParEspece<-function(x){
		 EM_Pi*apply(t(apply(x,1, function(y) dbinom(y[1],y[2],prob=1/(1+exp(-y[3:(K+2)]))))),2,prod)
	}
 
	
	#norme
    norme<- function(x) max(abs(x))

 
    #Initialisation des parametres
	EM_theta<- t(replicate(K, rnorm(L+1)))
	EM_Pi<- diff(c(0,sort(runif(K-1)),1))
    #Séparation des donnees par espece
    listDparEspece <- lapply(split(Tableau, Tableau[,3]),function(x) matrix(x,ncol=L+3))
 
 	#liste des matrices des produits scalaires initiaux par espece
	listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,EM_theta))
		 
	#liste des vraisemblances initiales par espece
	listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
	#Initialisation de la fonction a maximiser
	NQ<- 0

				 
	#Mise en forme des vecteurs de donnees des especes avec 1 comme premiere composante
	listDparEspece1vaExpl<-lapply(listDparEspece, function(x) t(apply(x,1, function(y) c(y[c(1:2)],1,y[4:(L+3)]))))

	debut<-Sys.time()
	iteration<-0
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

			 
			 #optimisation des theta_k, avec theta_k,m 
			 #comme point de depart de l'optimisation
			 
			 #determination du maximum par la methode de Newton-Raphson
			 EM_theta<-cbind(1:K, EM_theta)
			 EM_theta<-NewtonRaph(EM_theta)
			 
		
		### calcul de la nouvelle valeur de la fonction a maximiser : Q(theta|theta(m)) 
						#liste des matrices de produits scalaires par espece
						listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,EM_theta))
		 
						#liste des vraisemblances par espece
						listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
						
						NQ<-sum(rowSums(probaAposteriori*log(t(sapply(listVraisParEspece, function(x) x)))))
		### Condition d'arret
			iteration<-iteration+1
			if(abs(NQ-AQ)<0.0001) break
	}
	
			#log-vraisemblance
			logVrais<-sum(sapply(listVraisParEspece, function(x) log(sum(x))))
			### Regroupement
			EM_groupe<- apply(probaAposteriori, 1, which.max)
			png("classification.png")
			plot(as.numeric(names(EM_groupe)), EM_groupe,main="regroupement",xlab="Espèce",ylab="groupe")
            abline(v=10.5)
            abline(v=20.5)
            dev.off()
		fin<-Sys.time()
		temps<-fin-debut
return(list(Classification=EM_groupe, Proportions=EM_Pi, Parametres=EM_theta, Iterations=iteration, 
			log_vrasemblance=logVrais ,Visualisation="classification.png", Temps=temps))
 	
}	
