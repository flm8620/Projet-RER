# Projet-RER

##Mode d'emploi

1. lancer le programme sans data.txt  
      Model.exe doit apparaître
2. Mettre data.txt à côté de Model.exe  
      data.txt contient les données d'observation  
      **ATTENTION**  
      les données fournis dans data est la résultat
      Repartition(uTest),  
      où uTest={indice Destination= 5 15 22  
                Proportion Desti= 0.33 0.33 0.33  
                Proportion Sorti= 0.5 0.5  
                }  
      il faut changer les données selon notre comptage  

3. indiquer u0 dans main()  
4. préciser les paramètre de model  
      changer les paramètres  
      dans StationModel::StationModel(...){...}  
5. (facultatif) préciser les paramètre pour l'algorithme  
      5.1 préciser double seuil=???;  
         dans Optimisation::minimiser(...)  

      5.2 préciser const int maxIter=???;  
               et  double seuil=???;  
         dans Optimisation::minimiserSurPropo(...){...}  
      5.3 préciser double lambda=???;  
         dans Optimisation::unPasPourProportions(...){...}  
6. lancer le programme et attendre  


##PROBLEME d'OPTIMISTION

Le pb de optimisation avec contrainte :

                           inf    J(u) ,   K dans V    
                         u dans K  
u = {IndiceDesti[i],ProportionDesti[i],ProportionSorti[i]}  
J(u) = || reparti(u)-observation ||_l^2 = Sigma_i( ( reparti(u)[i]-observation[i] )^2 )  
K={u dans V | 0<=IndiceDesti[i]<nbPortes , Sigma_i(ProportionDesti[i])=1 , Sigma_i(ProportionSorti[i])=1 }  

##Référence du cours CALCUL SC  
regarder Page 69 du poly Calcul scientifique   
3.4.3 Algorithme de gradient (à pas fix) avec projection  
ou Page 17/21 du poly du cours "méthode numérique pour l'optimisation"  

##Algorithme de ce programme:

####1. Indiquer u_0 pour commencer  
    il faut préciser:  
        indicePorteDestination  
        proportion destination   
        proportion sortie  
    ATTENTION, un changement:  
        indicePorteSortie est défini dans la constructeur de StationModel:  
        stationmodel.cpp: StationModel::StationModel(...){...}  
   
####2. Calcul la direction de descente pour indicePorteDesti et l'appliquer  
      uGauche=u_k;  
      uGauche.indiceDesti[i]-=1;  
      uDroite=u_k;  
      uDroite.indiceDesti[i]+=1;  
      
#####2.1 Minimiser J(u_k) en changeant seulement les variables de proportion(variable réelle):  
      
            uPropo_k=u_k;  
      
######2.1.1 calculer la direction de descente pour les Proportions   
         
            du= - grad_Propo ( J(uPropo_k) )  
######2.1.2 appliquer un pas de gradient non-projeté   
         
            uPropo_k+1 = uPropo_k + lambda * du  
######2.1.3 projeter l'état proposé (contrainte: Sigma_i(ProportionDesti[i])=1 , Sigma_i(ProportionSorti[i])=1)   
            uPropo_k+1 = projeter(uPropo_k+1)  
######2.1.4 test de convergence:
         
                |J(uPropo_k+1)-J(uPropo_k)|  
            si --------------------------- > epsilon  
                         J(u_k)  
                        
               uPropo_k=uPropo_k+1;  
               goto 2.1.1
            sinon, c-a-d solution trouvée  
               u_k = u_Propo_k+1  
      
#####2.2 Pareil, minimiser J(uGauche) et J(uDroite) comme ci-dessus  
      
#####2.3 appliquer un pas sur indiceDesti[i]  
         on compare les trois valeurs: J(u_k), J(uGauche) et J(uDroite):  
         si J(u_k) est le minimum, on change rien, on note direction[i]=0  
         si J(uGauche) est le minimum on note: direction[i]= -1  
         si J(uDroite) est le minimum on fait: direction[i]= +1  

####3. On applique les pas sur indiceDesti[i]:  
      u_k+1=u_k;  
      Pour chaque destination indiceDesti[i]:  
         u_k+1.indiceDesti[i] = u_k.indiceDesti[i] + direction[i]  
         
####4. test de convergence:  
          |J(u_k+1)-J(u_k)|  
      si ------------------ > epsilon  
            J(u_k)  
                         
            u_k=u_k+1;  
            goto 2.  
      sinon, c-a-d solution trouvée  
            on s'arrête et affiche la solution u_k+1  
