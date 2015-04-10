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


int main(){

    ifstream file("data.txt");
    int nombrePortes=30;
    int nombreDesti=3;
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
    openWindow(620, 600);
    StationModel model(nombrePortes,nombreDesti);
    Optimisation optim;
    optim.setModel(model);
    optim.setObservation(observations);
    //start point:
    Variables uStart;
    uStart.desti.assign(3,0.0);
    uStart.desti[0]=12;
    uStart.desti[1]=14;
    uStart.desti[2]=16;
    uStart.propoDesti.assign(3,0.0);
    uStart.propoDesti[0]=1;
    uStart.propoDesti[1]=1;
    uStart.propoDesti[2]=1;
    uStart.propoSorti.assign(2,0.0);
    uStart.propoSorti[0]=1;
    uStart.propoSorti[1]=1;
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
    indicesDesti[0]=10;    indicesDesti[1]=20;    indicesDesti[2]=25;


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
