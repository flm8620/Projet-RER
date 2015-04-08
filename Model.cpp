#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <cmath>
#include "normaldistribution.h"
using namespace Imagine;
using namespace std;


const int nbPortes = 30;
const int nbPrevoyants = 1000;
const int nbConfort = 500;
const int tempsStation = 30;
const double propStresse = 0.8;
const int debitEntree = 50;
const double ecartType = 1;
const int nbDestinations = 3;
const double sigma = 5.0;
//tableau contenant le nombre de passagers dans le train
double s[30];
int destination[nbDestinations];
double proportionDestination[nbDestinations];
const int nbSorties=2;
int sorties[nbSorties];
//proportion des passagers qui arrive par chaque sortie
double proportionSorties[nbSorties];
//paramètre de la loi exponentielle
const int lambda = 1;


void init_sorties(){
    /*int nbsorties;
    cout << "Donnez le nombre de sortie de la gare (entier entre 1 et 30)" << endl;
    cin >> nbsorties;
    int sorties[nbsorties];*/
    cout << "Designez la porte la plus pres de la sortie de la gare" << endl;
    cout << "Donnez la proportion de gens qui arrive par la sortie" << endl;
    double S=0;
    for (int i = 0; i < nbSorties; i++){
        cout << "pour la sortie" << i+1 << " (entier entre 0 et 29)" << endl;
        cin >> sorties[i];
        cout << "un double quelconque" << endl;
        cin >> proportionSorties[i];
        S += proportionSorties[i];
    }
    for (int i = 0; i < nbSorties; i++){
        proportionSorties[i] /= S;
    }
}
void init_Destination(){
    cout << "Designez la porte la plus pres de chaque destination favorite" << endl;
    cout << "Donnez la proportion de gens qui veulent y aller" << endl;
    double S=0;
    for (int i = 0; i < nbDestinations; i++){
        cout << "pour la porte" << i + 1 << " (un entier entre 0 et 29)" << endl;
        cin >> destination[i];
        cout << "un double quelconque" << endl;
        cin >> proportionDestination[i];
        S+=proportionDestination[i];
    }
    for (int i = 0; i < nbDestinations; i++){
        proportionDestination[i] /= S;
    }
}

void gotominlocal(){
    //some problem of index detected, need to fix --Leman
    for (int i = 0; i < nbSorties; i++){
        for (int j = 0; j < proportionSorties[i] * nbConfort; j++){
            if (sorties[i] == 1){
                for (int k = 1; k < nbPortes; k++){
                    if (s[k + 1] >= s[k]){ s[k] += 1; }
                }
            }
            if (sorties[i] == nbPortes){
                for (int k = nbPortes-1; k > 0; k--){
                    if (s[k - 1] >= s[k]){ s[k] += 1; }
                }
            }
            if (sorties[i] < nbPortes && sorties[i]>0){
                for (int k = sorties[i]; k < max( nbPortes - sorties[i], sorties[i] ); k++){
                    if (s[sorties[i] + k + 1] >= s[sorties[i] + k]){ s[sorties[i] + k] += 1; }
                    else if (s[sorties[i] - k - 1] >= s[sorties[i] - k]){ s[sorties[i] - k] += 1; }
                }
            }
        }
    }
}

void loiexponentiel(int sorties[nbSorties]){
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

void loiuniforme(int sorties[nbSorties],int destination[nbDestinations]){
    //boucle sur toutes les sorties
    for (int i = 0; i < nbSorties; i++){
        //boucle sur les destinations
        for (int j = 0; j < nbDestinations; j++){
            //boucle sur les portes entre la sortie et les destinations
            int nbInterval=abs(destination[j]-sorties[i])+1;
            for (int k = min(sorties[i],destination[j]); k < max(sorties[i], destination[j])+1; k++){
                s[k] += tempsStation*debitEntree*(1 - propStresse)*proportionSorties[i] * proportionDestination[j]/nbInterval;
            }
        }
    }
}

void loinormale(int i, int j){
    //i - No.destination
    //j - No.porte
    NormalDistribution normal;
    double xSup=(j+0.5-destination[i])/sigma;
    double xInf=(j-0.5-destination[i])/sigma;
    double value=nbPrevoyants*proportionDestination[i] * (normal.getPhi(xSup) - normal.getPhi(xInf));
    s[j] += value;
}


void init_quai(){
    for (int i = 0; i < nbPortes; i++){
        s[i] = 0;
    }
}


int main2(){
    init_quai();

    init_Destination();
    init_sorties();

    for (int i = 0; i < nbDestinations; i++){
        for (int j = 0; j < nbPortes; j++){
            loinormale(i, j);
        }
    }
    /**gotominlocal();
    for (int i = 0; i < nbsorties; i++){
        for (int j = 0; j < proportionsorties[i] * propstresse*tempsstation / debitentree; j++){
            s[sorties[i]] += 1;
        }
    }
    for (int i = 0; i < nbsorties; i++){
        for (int j = 0; j < nbdestinations; j++){
            for (int k = 0; k < abs(sorties[i] - Destination[j]); j++){
                s[k + sorties[i]] += proportionsorties[i] * (1 - )*proportiondestination[j] * tempsstation / debitentree;
            }
        }
    }**/
    loiexponentiel(sorties);
    loiuniforme(sorties,destination);
    openWindow(620, 600);
    for (int i = 0; i < nbPortes; i++){
        drawLine(20 * i+10, 599, 20 * i+10, 599.-s[i], MAGENTA);
        cout << i<<'\t'<< s[i] << endl;
    }
    system("pause");
    return 0;
}
