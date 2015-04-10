#include "stationmodel.h"
#include <iostream>
#include "normaldistribution.h"
using namespace std;

void StationModel::normalizeProportion()
{
    // assurer que la sommation de proportion = 1
    double S=0;
    for (int i = 0; i < nbDestinations; i++){
        S+=proportionDestination[i];
    }
    if(S==0){
        cout<<"error: sum of proportion =0"<<endl;
        exit(0);
    }
    for (int i = 0; i < nbDestinations; i++){
        proportionDestination[i] /=S;
    }
    S=0;
    for (int i = 0; i < nbSorties; i++){
        S+=proportionSorties[i];
    }
    if(S==0){
        cout<<"error: sum of proportion =0"<<endl;
        exit(0);
    }
    for (int i = 0; i < nbSorties; i++){
        proportionSorties[i] /=S;
    }
}

void StationModel::loiNormal()
{
    NormalDistribution normal;
    vector<double> normalDiscret,normalDiscretNormee;
    normalDiscret.assign(nbPortes,0);
    normalDiscretNormee.assign(nbPortes,0);
    for (int i = 0; i < nbDestinations; i++){
        for (int j = 0; j < nbPortes; j++){
            double xSup=(j+0.5-destinations[i])/sigma;
            double xInf=(j-0.5-destinations[i])/sigma;
            normalDiscret[j]=(normal.getPhi(xSup) - normal.getPhi(xInf));
        }
        double S=0;
        for(int j=0;j<nbPortes;j++){
            S+=normalDiscret[j];
        }
        for(int j=0;j<nbPortes;j++){
            normalDiscretNormee[j]=normalDiscret[j]/S;

        }
        for(int j=0;j<nbPortes;j++){

        double value=nbPrevoyants*proportionDestination[i] * normalDiscretNormee[j];
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
    for (int i = 0; i < nbSorties; i++){
        for (int j = 0; j < proportionSorties[i] * nbConfort; j++){
            int k,k_left,k_right;
            if (sorties[i] == 0){
            //  O->----------------------
                for (k = 0; k < nbPortes-1; k++){
                    if (s[k + 1] >= s[k]){
                        s[k] += 1;
                        break;
                    }
                }
                //derniere porte
                if(k==nbPortes-1)s[nbPortes-1]+=1;
            }else if (sorties[i] == nbPortes-1){
            //  ----------------------<-O
                for (k = nbPortes-1; k > 0; k--){
                    if (s[k - 1] >= s[k]){
                        s[k] += 1;
                       break;
                    }
                }
                //derniere porte
                if(k==0)s[0]+=1;
            }else if (sorties[i] < nbPortes && sorties[i]>0){
            //  ---------<-O->-----------

                //Il y a deux choix, soit vers gauche, soit droite
                //On suppose que la proportion de choix gauche et de droite est
                //proportionnelle que la longueur de quai à gauche et à droite
                double partieVersGauche,partieVersDroite;
                partieVersGauche=double(sorties[i]-0)/(nbPortes-1);
                partieVersDroite=double(nbPortes-sorties[i])/(nbPortes-1);
                //Gauche+Droite=1

                //vers gauche
                for (k_left=sorties[i]; k_left >0; k_left--){
                    if(s[k_left-1]>=s[k_left]){
                        s[k_left]+=partieVersGauche;
                        break;
                    }
                }
                //derniere porte
                if(k==0)s[0]+=partieVersGauche;

                //vers droite
                for (k_right=sorties[i]; k_right <nbPortes-1; k_right++){
                    if(s[k_right-1]>=s[k_right]){
                        s[k_right]+=partieVersDroite;
                        break;
                    }
                }
                //derniere porte
                if(k==nbPortes-1)s[nbPortes-1]+=partieVersDroite;
            }
        }
    }
}

void StationModel::initSortiesDesti(std::vector<int> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti)
{
    int nDesti=indicesDesti.size();
    int npDesti=proportionDesti.size();
    int npSorti=proportionSorti.size();
    if(nDesti!=nbDestinations){
        cout<<"error: size of indicesDestination incorrect"<<endl;
        exit(0);
    }
    if(npDesti!=nbDestinations){
        cout<<"error: size of proportionDestination incorrect"<<endl;
        exit(0);
    }
    if(npSorti!=nbSorties){
        cout<<"error: size of proportionSortie incorrect"<<endl;
        exit(0);
    }
    //assigner les indice de porte destination
    for(int i=0;i<nbDestinations;i++){
        destinations[i]=indicesDesti[i];
    }

    //assigner les proportions
    for (int i = 0; i < nbDestinations; i++){
        proportionDestination[i] =proportionDesti[i];
    }
    for (int i = 0; i < nbSorties; i++){
        proportionSorties[i] =proportionSorti[i];
    }





}

StationModel::StationModel(int nbPortes, int nbDestinations)
{
    this->nbPortes=nbPortes;
    this->nbDestinations=nbDestinations;

    //par default
    this->nbSorties=2;


    s.assign(nbPortes,0.0);
    destinations.assign(nbDestinations,0);
    sorties.assign(nbSorties,0);

    //par default

    sorties[0]=0;
    sorties[1]=nbPortes-1;

    //par default, proportion uniforme
    proportionDestination.assign(nbDestinations,1.0);
    proportionSorties.assign(nbSorties,1.0);

    //par default
    sigma=2.0;
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

std::vector<double> StationModel::getRepartition(std::vector<int> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti)
{
    //init_quai
    s.assign(nbPortes,0.0);
    initSortiesDesti(indicesDesti,proportionDesti,proportionSorti);
    normalizeProportion();
    loiNormal();
    loiExponentiel();
    loiUniforme();
    gotominilocal();
    return s;
}

std::vector<double> StationModel::getRepartitionNonNormalizePropo(std::vector<int> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti)
{
    //init_quai
    s.assign(nbPortes,0.0);
    initSortiesDesti(indicesDesti,proportionDesti,proportionSorti);
    loiNormal();
    loiExponentiel();
    loiUniforme();
    gotominilocal();
    return s;
}

