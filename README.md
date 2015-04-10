# Projet-RER

Mettre data.txt à côté de Model.exe
data.txt contient les données d'observation

---------------PROBLEME d'OPTIMISTION-----------------
Le pb de optimisation avec contrainte :

   inf    J(u) ,   K dans V
 u dans K
u = {IndiceDesti[i],ProportionDesti[i],ProportionSorti[i]}
J(u) = || reparti(u)-observation ||_l^2 = Sigma_i( ( reparti(u)[i]-observation[i] )^2 )
K={u dans V | 0<=IndiceDesti[i]<nbPortes , Sigma_i(ProportionDesti[i])=1 , Sigma_i(ProportionSorti[i])=1 }

-----------référence du cours CALCUL SC-----------------------
regarder Page 69 du poly Calcul scientifique 
3.4.3 Algorithme de gradient (à pas fix) avec projection
ou Page 17/21 du poly du cours "méthode numérique pour l'optimisation"


------------Algorithme de ce programme:------------

1. Indiquer un u0 pour commencer 
    il faut préciser:
        indicePorteDestination
        proportion destination 
        proportion sortie
    ATTENTION, un changement:
        indicePorteSortie est défini dans la constructeur de StationModel:
        stationmodel.cpp: StationModel::StationModel(...){...}
2. Calcul la direction de descente pour indicePorteDesti
    
