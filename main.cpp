#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include "stationmodel.h"
#include "optimisation.h"
#include <vector>
#include <fstream>
using namespace Imagine;
using namespace std;


void printOutResultat(vector<double> R){

    double sum=0;
    for(int i=0; i<R.size();i++){
        drawLine(20 * i+10, 599, 20 * i+10, 599.-R[i], MAGENTA);
        cout << i<<'\t'<< R[i] << endl;
        sum+=R[i];
    }
    cout<<"sommation: "<<sum<<endl;

}
std::vector<double> readFile(const char* fileName,int nombrePortes){
    ifstream file(fileName);
    vector<double> observations(nombrePortes);
    cout<<"Lire ficher 'data.txt' ..."<<endl;
    if(!file){
        cout<<"can't open file 'data.txt' !"<<endl;
        exit(0);
    }
    for(int i=0;i<nombrePortes;i++){
        file>>observations[i];
    }
    if(!file){
        cout<<"can't get enough data"<<endl;
        exit(0);
    }
    for(int i=0;i<nombrePortes;i++){
        cout<<"Porte "<<i<<" "<<observations[i]<<endl;
    }
    cout<<endl;
    return observations;
}

int main(){


    int nombrePortes=30;
    int nombreDesti=3;
    //lire les données dans 'data.txt'
    vector<double> observations=readFile("data.txt",nombrePortes);


    StationModel model(nombrePortes,nombreDesti);
    //
    // le pb de optimisation avec contrainte :
    //
    //     inf    J(u) ,   K dans V
    //  u dans K
    //
    //    u = {IndiceDesti[i],ProportionDesti[i],ProportionSorti[i]}
    //    J(u) = || reparti(u)-observation ||_l^2 = Sigma_i( ( reparti(u)[i]-observation[i] )^2 )
    //    K={u dans V | 0<=IndiceDesti[i]<nbPortes , Sigma_i(ProportionDesti[i])=1 , Sigma_i(ProportionSorti[i])=1 }
    //
    //
    Optimisation optim;
    //faire optim savoir le model qu'on utilise
    optim.setModel(model);
    //faire optim savoir les données qu'on a lit
    optim.setObservation(observations);

    //Créer une structure de variable
    //choisir u0 dans K
    //IL FAUT ESSAYER DIFFERENT u0 POUR EVITER MIN LOCAL
    Variables uStart;
    uStart.desti.assign(3,0.0);
    uStart.desti[0]=10;
    uStart.desti[1]=15;
    uStart.desti[2]=20;
    uStart.propoDesti.assign(3,0.0);
    uStart.propoDesti[0]=0.3;
    uStart.propoDesti[1]=0.3;
    uStart.propoDesti[2]=0.3;
    uStart.propoSorti.assign(2,0.0);
    uStart.propoSorti[0]=0.5;
    uStart.propoSorti[1]=0.5;


    openWindow(620, 700);

    //minimisation
    //regarder Page 69 du poly Calcul Scientifique
    optim.minimiser(uStart);

    endGraphics();
}

int main2(){
    //30 portes, 3 destinations 2 sorties
    //pour changer les autre parametre, voir la constructeur de StationModel
    int nombrePortes=30;
    int nombreDesti=3;
    //int nombreSorti=2;
    StationModel model(nombrePortes,nombreDesti);

    vector<int>    indicesDesti(3);
    vector<double> proportionDesti(3);
    vector<double> proportionSorti(2);


//indices de portes pour Destinations et Sorties
//    Desti:  10  20  25           Sortie: 0  29
    indicesDesti[0]=5;    indicesDesti[1]=15;    indicesDesti[2]=22;


//proportions pour Destinations et Sorties
//    Desti: 1/3  1/3  1/3         Sortie: 1/2  1/2
//Ils vont être normalisés après.
    proportionDesti[0]=1.0;
    proportionDesti[1]=1.0;
    proportionDesti[2]=1.0;

    proportionSorti[0]=1.0;
    proportionSorti[1]=1.0;






    vector<double> resultat;
    resultat=model.getRepartition(indicesDesti,proportionDesti,proportionSorti);
    openWindow(620, 600);
    printOutResultat(resultat);
    endGraphics();
    return 0;

}
