#include "stationmodel.h"
#include <iostream>
#include "normaldistribution.h"
using namespace std;

void StationModel::loiNormal()
{
    NormalDistribution normal;
    for (int i = 0; i < nbDestinations; i++){
        for (int j = 0; j < nbPortes; j++){
            double xSup=(j+0.5-destinations[i])/sigma;
            double xInf=(j-0.5-destinations[i])/sigma;
            double value=nbPrevoyants*proportionDestination[i] * (normal.getPhi(xSup) - normal.getPhi(xInf));
            s[j] += value;
        }
    }
}

void StationModel::loiExponentiel()
{
    //boucle sur toutes les sorties (1 loi par sortie)
    for (int i = 0; i < nbSorties; i++){
        //Cas sortie au début du quai
        if (sorties[i] == 0){
            for (int j = 0; j < nbPortes; j++){
                double A=(1 - exp(-lambda)) / (1 - exp(-lambda*nbPortes));
                //A est une constante pour normer la loi
                s[j] += A*exp(-lambda*j)   *   tempsStation*debitEntree*propStresse*proportionSorties[i];
                //la loi est normée puis elle est multpliée par le nombre de passager stressé et en retard entrant par la sortie i.
            }
        }
        //Cas sortie à la fin du quai
        if (sorties[i] == nbPortes - 1){
            for (int j = 0; j < nbPortes; j++){
                double A=(1 - exp(-lambda)) / (1 - exp(-lambda*nbPortes));
                s[j] += A*exp(-lambda*(nbPortes - 1 - j))*tempsStation*debitEntree*propStresse*proportionSorties[i];
            }
        }
        //Cas sortie au milieu du quai
        if (sorties[i]>0 && sorties[i] < nbPortes-1){
            for (int j = 0; j < sorties[i]; j++){
                double A=(1 - exp(-lambda)) / (1 + exp(1) - exp(-lambda*sorties[i]) - exp(-lambda*(nbPortes - sorties[i])));
                s[j] += A*exp(-lambda*(sorties[i] - 1 - j))*tempsStation*debitEntree*propStresse*proportionSorties[i];
            }
            for (int j = sorties[i]; j < nbPortes; j++){
                double A=(1 - exp(-lambda)) / (2 - exp(-lambda*sorties[i]) - exp(-lambda*(nbPortes - sorties[i])));
                s[j] += A*exp(-lambda*(j + 1 - sorties[i]))*tempsStation*debitEntree*propStresse*proportionSorties[i];
            }
            //cette expression est la synthèse de trois étapes
            //1/ à la porte la plus proche de la sortie, la fonction a pour valeur 1 et elle décroit de chaque côté selon une même loi exponentielle de paramètre -lambda
            //2/ je norme cette fonction sur le quai
            //3/ je multiplie par le nombre de passager stressé en retard entrant par la sortie i.
        }
    }
}

void StationModel::loiUniforme()
{
    //boucle sur toutes les sorties
    for (int i = 0; i < nbSorties; i++){
        //boucle sur les destinations
        for (int j = 0; j < nbDestinations; j++){
            //boucle sur les portes entre la sortie et les destinations
            int nbInterval=abs(destinations[j]-sorties[i])+1;
            for (int k = min(sorties[i],destinations[j]); k < max(sorties[i], destinations[j])+1; k++){
                s[k] += tempsStation*debitEntree*(1 - propStresse)*proportionSorties[i] * proportionDestination[j]/nbInterval;
            }
        }
    }
}

void StationModel::gotominilocal()
{
    //NEED TO BE FIXED
    for (int i = 0; i < nbSorties; i++){
        for (int j = 0; j < proportionSorties[i] * nbConfort; j++){
            if (sorties[i] == 0){
            //  O->----------------------
                for (int k = 0; k < nbPortes-1; k++){
                    if (s[k + 1] >= s[k]){
                        s[k] += 1;
                        continue;
                    }
                }
                //derniere porte
                s[nbPortes-1]+=1;
            }else if (sorties[i] == nbPortes-1){
            //  ----------------------<-O
                for (int k = nbPortes-1; k > 0; k--){
                    if (s[k - 1] >= s[k]){
                        s[k] += 1;
                        continue;
                    }
                }
                //derniere porte
                s[0]+=1;
            }else if (sorties[i] < nbPortes && sorties[i]>0){
            //  ---------<-O->-----------
                for (int k_left=sorties[i], k_right = sorties[i]; k_left < max( nbPortes - sorties[i], sorties[i] ); k_left++,k_right++){
                    //WAIT TO BE FINISHED
                    //the passenger goes in two direction in the same time?
                }
            }
        }
    }
}

void StationModel::initSortiesDesti(std::vector<int> &indicesPortes, std::vector<double> &proportions)
{
    int n1=indicesPortes.size();
    int n2=proportions.size();
    if(n1!=nbSorties+nbDestinations){
        cout<<"error: size of indicesPortes incorrect"<<endl;
        exit(0);
    }
    if(n2!=nbSorties+nbDestinations){
        cout<<"error: size of proportion incorrect"<<endl;
        exit(0);
    }
    //assigner les indice de porte
    for(int i=0;i<nbDestinations;i++){
        destinations[i]=indicesPortes[i];
    }
    for(int i=0;i<nbSorties;i++){
        sorties[i]=indicesPortes[i+nbDestinations];
    }
    // assurer que la sommation de proportion = 1
    double S=0;
    for (int i = 0; i < nbDestinations; i++){
        S+=proportions[i];
    }
    for (int i = 0; i < nbDestinations; i++){
        proportionDestination[i] =proportions[i]/S;
    }

    S=0;
    for (int i = 0; i < nbSorties; i++){
        S+=proportions[i+nbDestinations];
    }
    for (int i = 0; i < nbSorties; i++){
        proportionSorties[i] =proportions[i+nbDestinations]/S;
    }




}

StationModel::StationModel(int nbPortes, int nbDestinations, int nbSorties)
{
    this->nbPortes=nbPortes;
    this->nbDestinations=nbDestinations;
    this->nbSorties=nbSorties;
    s.assign(nbPortes,0.0);
    destinations.assign(nbDestinations,0);
    sorties.assign(nbSorties,0);
    //par default, proportion uniforme
    proportionDestination.assign(nbDestinations,1.0);
    proportionSorties.assign(nbSorties,1.0);

    //par default
    sigma=5.0;
    lambda=1;
    nbPrevoyants=1000;
    nbConfort=500;
    tempsStation=30;
    propStresse=0.8;
    debitEntree=50;
    ecartType=1;



}

StationModel::~StationModel()
{

}

std::vector<double> StationModel::getRepartition(std::vector<int> &indicesPortes, std::vector<double> &proportions)
{
    //init_quai
    s.assign(nbPortes,0.0);
    initSortiesDesti(indicesPortes,proportions);
    loiNormal();
    loiExponentiel();
    loiUniforme();
    //gotominilocal(); ----> NEED TO BE FIXED
    return s;
}

