
EM_Recrutement<-function(Tableau, K){

##########################################################################
###		Estimation d'un melange de regression de Poisson par EM        ###
##########################################################################

 #nombre de variable explicative
 L<-ncol(Tableau)-2
 
 #fonction utilisee dans l'optimisation des theta pour le calcul des réels qui 
 #pondèrent les vecteurs lignes de données dans l'expression du gradient
 fvecteur<-function(x,u){
			 v<- x[-1]
			(-x[1]+exp(v%*%u))*v
		}
#fonction utilisee dans l'optimisation des theta pour le calcul des réels qui 
#pondèrent les vecteurs lignes de données dans l'expression de la hessienne		
 fvecteur2<- function(x,u){
			 v <- x[-1]		
		     (exp(v%*%u)* v) %*% t(v)
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
 NewtonRaph<-function(u_zero,k){
	u_plus<-u_zero
	nabla_plus<-nablaJ(u_zero,k)
	repeat{
		u_moins<-u_plus
		u_plus<- u_moins - solve(hessienJ(u_moins,k), nablaJ(u_moins,k))
		if(norme(u_plus-u_moins)<0.000001) break
	}
	u_plus
 }
 #calcul des scalaires par espece
 fMatriceScalaire<-function(x,M){
		t(apply(x,1,function(y) c(y[1],M%*%c(1,y[-c(1,2)]))))
		}
 #calcul des vraisemblances par espece
 fVraisParEspece<-function(x){
		 EM_Pi*apply(t(apply(x,1, function(y) dpois(y[1],lambda=exp(y[-1])))),2,prod)
	}
 
 #recherche lineaire de Wolfe par interpolation cubique

  #calcul des valeurs de Q (fonction a optimiser) avec le pas
  yy_plus<-function(Point,Point_plus){
		######
		listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,Point))
		listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
		Q<--sum(rowSums(probaAposteriori*log(t(sapply(listVraisParEspece, function(x) x)))))

		######
		
		listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,Point_plus))
		listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
		Q_plus<--sum(rowSums(probaAposteriori*log(t(sapply(listVraisParEspece, function(x) x)))))
		
		c(Q, Q_plus)
	}

   # test des conditions de wolfe
	wolfe<-function(vrais, rho, nabla2, nablaPoint_plusNabla){
		eps1<-0.0001
		eps2<-0.99
		signal<-0
		wolf21 <-  nablaPoint_plusNabla
		wolf13 <- -eps1*rho*nabla2
		wolf22 <- -eps2*nabla2		
		
		if(vrais[2] <= vrais[1] + wolf13 && wolf21 >= wolf22){
			signal<-1
		}
		signal
	}
	
	# recherche lineaire de wolfe par interpolation cubique
	#calcul du polynome et du pas local
	polySplineCubique<-function(vrais, rho, nabla2, nablaPoint_plusNabla){
	
		A<- matrix(c(rho^2, 2*rho, rho^3, 3*rho^2), nrow=2, ncol=2)
		b<- c(vrais[2]- vrais[1]+ rho*nabla2, nablaPoint_plusNabla + nabla2)
		#coefficients du polynome d'interpolation cubique
		P<- c(vrais[1], -nabla2, solve(A,b))
		#extremums
		delta<- P[3]^2 - 3*P[2]*P[4]
		if(delta>=0){
			extremum <- c( (-P[3]-sqrt(delta))/(3*P[4]), (-P[3]+sqrt(delta))/(3*P[4]) )
			minim<- ifelse(P[2]*extremum[1]+P[3]*extremum[1]^2+ P[4]*extremum[1]^3 < P[2]*extremum[2]+P[3]*extremum[2]^2+ P[4]*extremum[2]^3, extremum[1], extremum[2])
		}else{
			minim<-rho/2
		}
		minim
	}
	#recherche lineaire
	rechercheLin<-function(nabla,Point){
		rho<-1
		Point_plus<- Point - rho * nabla
		vrais<-yy_plus(Point,Point_plus)
		while( max(vrais)==Inf){
			rho<- runif(1,min=rho/2, max=rho)
			Point_plus<- Point - rho * nabla
			vrais<-yy_plus(Point, Point_plus)
		}
		nabla2<-sum(t(apply(nabla, 1, function(x) x%*%x)))
		Point_plus<- cbind(1:K, Point_plus)
		nablaPoint_plusNabla <- sum(t(apply(Point_plus, 1, function(x) (-nablaJ(x[-1], x[1]))%*%nabla[x[1],])))
		signal<-wolfe(vrais, rho, nabla2, nablaPoint_plusNabla)
		while(signal!=1 && abs(vrais[1]-vrais[2])>0.000000001){
			rho<-polySplineCubique(vrais,rho,nabla2,nablaPoint_plusNabla)
			Point_plus<- Point - rho * nabla
			vrais<-yy_plus(Point, Point_plus)
			Point_plus<- cbind(1:K, Point_plus)
			nablaPoint_plusNabla <- sum(t(apply(Point_plus, 1, function(x) (-nablaJ(x[-1], x[1]))%*%nabla[x[1],])))
			signal<-wolfe(vrais, rho, nabla2, nablaPoint_plusNabla)
		}
		rho
	}
	# algorithme de gradient pour initialisation des parametres
	algoGradient<-function(){
		for(i in 1:5){
			EM_theta<-cbind(1:K, EM_theta)
			nabla<- t(apply(EM_theta, 1, function(x) nablaJ( x[-1], x[1])))
			rho<- rechercheLin(nabla,EM_theta[,-1])
			EM_theta<-EM_theta[,-1] - rho * nabla
		}
		EM_theta
	}

 
 #Initialisation des parametres
	EM_theta<- t(replicate(K, rnorm(L+1)))
	EM_Pi<- diff(c(0,sort(runif(K-1)),1))
    #Séparation des donnees par espece
    listDparEspece <- lapply(split(Tableau, Tableau[,2]),function(x) matrix(x,ncol=L+2))
 
 	#liste des matrices des produits scalaires initiaux par espece
	listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,EM_theta))
		 
	#liste des vraisemblances initiales par espece
	listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
	#Initialisation de la fonction a maximiser
	NQ<- 0
	
	#Mise en forme des vecteurs de donnees des especes avec 1 comme premiere composante
	listDparEspece1vaExpl<-lapply(listDparEspece, function(x) t(apply(x,1, function(y) c(y[1],1,y[-c(1,2)]))))

	#debut<-Sys.time()
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
			 
			 #On commence par quelques iterations de la methode du gradient avec recherche lineaire de Wolfe par interpolation
			 #spline cubique afin garantir la convergence de la methode de Newton-Raphson, utilisee par la suite
			 if(iteration==0){
					EM_theta<- algoGradient()
			 }
			 #determination du maximum par la methode de Newton-Raphson
			 EM_theta<-cbind(1:K, EM_theta)
			 EM_theta<-t(apply(EM_theta,1,function(x) NewtonRaph(x[-1], x[1])))
			 
		
		### calcul de la nouvelle valeur de la fonction a maximiser : Q(theta|theta(m)) 
						#liste des matrices de produits scalaires par espece
						listMatriceScalaire<-lapply(listDparEspece, function(x) fMatriceScalaire(x,EM_theta))
		 
						#liste des vraisemblances par espece
						listVraisParEspece<- lapply(listMatriceScalaire, fVraisParEspece)
						
						NQ<-sum(rowSums(probaAposteriori*log(t(sapply(listVraisParEspece, function(x) x)))))
		### Condition d'arret
			iteration<-iteration+1
			if(abs(NQ-AQ)<0.000000001) break
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
				 
	return(list(Classification=EM_groupe, Proportions=EM_Pi, Parametres=EM_theta, Iterations=iteration, log_vrasemblance=logVrais ,Visualisation="classification.png"))
 

	
	#fin<-Sys.time()
			
	#(fin-debut)
}