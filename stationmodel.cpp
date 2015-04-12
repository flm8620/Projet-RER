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

    for (int i = 0; i < nbDestinations; i++){
        double nbVoyageur=nbPrevoyants*proportionDestination[i];
        normalDiscret.assign(nbPortes,0);
        normalDiscretNormee.assign(nbPortes,0);
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
            s[j] +=nbVoyageur * normalDiscretNormee[j];
        }
    }

}

void StationModel::loiExponentiel()
{


    //                             centre
    //          deux loi EXP:         _
    //                               /|\
    //                              / | \
    //                            _/  |  \_
    //                          _/|   |    \_
    //                      __--  |   |     |--__
    //                ___---      |   |     |    ---____
    //         ___----  |         |   |     |         | ----______
    //         |        |         |   |     |         |         | -------|
    //      -------|---------|---------|---------|---------|---------|------>  x
    // x :         0         1         2         3         4         5
    //          _______   _______   _______   _______   _______   _______
    //train :  /_______|-|_______|-|_______|-|_______|-|_______|-|_______\
    //boucle sur toutes les sorties (1 loi par sortie)
    vector<double> ExpDiscret,ExpDiscretNormee;

    for (int i = 0; i < nbSorties; i++){
        ExpDiscret.assign(nbPortes,0);
        ExpDiscretNormee.assign(nbPortes,0);
        double centre=sorties[i];
        int Icentre=floor(centre+0.5);
        Icentre=min(nbPortes-1,Icentre);//if centre=29.5, Icentre will be 30, out of index
        double nbVoyageur=nbRetardStresse*proportionSorties[i];
        // if(k==centre):
        ExpDiscret[Icentre]+=exp(-lambda*0)-exp(-lambda*(centre-(Icentre-0.5)))
                + exp(-lambda*0) - exp(-lambda*(Icentre+0.5-centre));
        for(int k=0;k<Icentre;k++){
            ExpDiscret[k]+=exp(-lambda*(centre-(k+0.5)))-exp(-lambda*(centre-(k-0.5)));
        }
        for(int k=Icentre+1;k<nbPortes;k++){
            ExpDiscret[k]+=exp(-lambda*((k-0.5)-centre))-exp(-lambda*((k+0.5)-centre));
        }
        double S=0;
        for(int j=0;j<nbPortes;j++){
            S+=ExpDiscret[j];
        }
        for(int j=0;j<nbPortes;j++){
            ExpDiscretNormee[j]=ExpDiscret[j]/S;
        }
        for(int j=0;j<nbPortes;j++){
            s[j] += nbVoyageur * ExpDiscretNormee[j];
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
            double distanceDestiSorti=abs(destinations[j]-sorties[i]);
            double left=min(sorties[i],destinations[j]);
            double right=max(sorties[i], destinations[j]);
            int Ibegin=floor(left+0.5);
            int Iend=floor(right+0.5);
            Ibegin=min(nbPortes-1,Ibegin);
            Iend=min(nbPortes-1,Iend);
            double nbVoyageur=nbRetardNonStresse*proportionSorties[i] * proportionDestination[j];
            //
            //
            //
            //            sortie[i]                             desti[i]
            //               |                                    |
            //               |                                    |
            //          _____V_   _______   _______   _______   __V____   _______
            //train :  /_______|-|_______|-|_______|-|_______|-|_______|-|_______\
            // k :         0         1         2         3         4         5
            //             ^                                       ^
            //          Ibegin                                   Iend
            //
            //longeur:       ---|---------|---------|---------|----
            //
            //         ^
            //         |     --------------------------------------  ->constant = 1/distaceDestiSorti
            //s[k]:    |     |  |         |         |         |   |
            //         |-----------------------------------------------------------> X
            if(Ibegin==Iend){
                s[Ibegin]+=nbVoyageur;
                continue;
            }
            for (int k = Ibegin; k <= Iend; k++){
                double longeur;
                if(left>=k-0.5&&left<k+0.5){
                    longeur=k+0.5-left;
                }else if(right>=k-0.5&&right<k-0.5){
                    longeur=right-(k-0.5);
                }else{
                    longeur=1;
                }
                s[k] += longeur/distanceDestiSorti*nbVoyageur;
            }
        }
    }
}

void StationModel::gotominilocal()
{
    double unPasPourJ=0.1;
    for (int i = 0; i < nbSorties; i++){
        for (double j = 0; j < proportionSorties[i] * nbConfort; j+=unPasPourJ){
            int k,k_left,k_right;


            //Il y a deux choix, soit vers gauche, soit droite
            //On suppose que la proportion de choix gauche et de droite est
            //proportionnelle que la longueur de quai à gauche et à droite
            double partieVersGauche,partieVersDroite;
            partieVersGauche=(sorties[i]-0.0)/(nbPortes-1);
            partieVersGauche=max(min(30.0,partieVersGauche),0.0);
            partieVersDroite=(nbPortes-1-sorties[i])/(nbPortes-1);
            partieVersGauche=max(min(30.0,partieVersGauche),0.0);
            //Gauche+Droite=1
            int Icentre=floor(sorties[i]+0.5);
            Icentre=min(nbPortes-1,Icentre);
            //vers gauche
            for (k_left=Icentre; k_left >0; k_left--){
                if(s[k_left-1]>=s[k_left]){
                    s[k_left]+=partieVersGauche*unPasPourJ;
                    break;
                }
            }
            //derniere porte
            if(k==0)s[0]+=partieVersGauche*unPasPourJ;

            //vers droite
            for (k_right=Icentre; k_right <nbPortes-1; k_right++){
                if(s[k_right+1]>=s[k_right]){
                    s[k_right]+=partieVersDroite*unPasPourJ;
                    break;
                }
            }
            //derniere porte
            if(k==nbPortes-1)s[nbPortes-1]+=partieVersDroite*unPasPourJ;

        }
    }
}

void StationModel::initSortiesDesti(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti)
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
    sorties[0]=0;            // au debut de wagon 0
    sorties[1]=(nbPortes-1)+0.5;// au bout de wagon 29

    //par default, proportion uniforme
    proportionDestination.assign(nbDestinations,1.0);
    proportionSorties.assign(nbSorties,1.0);

    //par default
    sigma=2.5;
    lambda=0.4;
    nbVoyageurTotal=1000;
    //nbPrevoyants=500;
    //nbConfort=200;
    //nbPrevoyants=0;
    //nbConfort=0;
    //propStresse=0.8;
    //ecartType=1;

}

StationModel::~StationModel()
{

}

std::vector<double> StationModel::getRepartition(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti,std::vector<double> propoVoyageur)
{
    //init_quai
    s.assign(nbPortes,0.0);
    nbPrevoyants=nbVoyageurTotal*propoVoyageur[0];
    nbConfort=nbVoyageurTotal*propoVoyageur[1];
    nbRetardStresse=nbVoyageurTotal*propoVoyageur[2];
    nbRetardNonStresse=nbVoyageurTotal*propoVoyageur[3];
    initSortiesDesti(indicesDesti,proportionDesti,proportionSorti);
    normalizeProportion();
    loiNormal();
    gotominilocal();
    loiExponentiel();
    loiUniforme();

    return s;
}

std::vector<double> StationModel::getRepartitionNonNormalizePropo(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti,std::vector<double> propoVoyageur)
{
    //init_quai
    s.assign(nbPortes,0.0);
    nbPrevoyants=nbVoyageurTotal*propoVoyageur[0];
    nbConfort=nbVoyageurTotal*propoVoyageur[1];
    nbRetardStresse=nbVoyageurTotal*propoVoyageur[2];
    nbRetardNonStresse=nbVoyageurTotal*propoVoyageur[3];
    initSortiesDesti(indicesDesti,proportionDesti,proportionSorti);
    loiNormal();
    gotominilocal();
    loiExponentiel();
    loiUniforme();

    return s;
}

