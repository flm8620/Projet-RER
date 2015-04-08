#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include "stationmodel.h"
#include <vector>
using namespace Imagine;
using namespace std;


void printOutResultat(vector<double> R){
    openWindow(620, 600);
    double sum=0;
    for(int i=0; i<R.size();i++){
        drawLine(20 * i+10, 599, 20 * i+10, 599.-R[i], MAGENTA);
        cout << i<<'\t'<< R[i] << endl;
        sum+=R[i];
    }
    cout<<"sommation: "<<sum<<endl;
    endGraphics();
}

int main(){
    //30 portes, 3 destinations 2 sorties
    //pour changer les autre parametre, voir la constructeur de StationModel
    int nombrePortes=30;
    int nombreDesti=3;
    int nombreSorti=2;
    StationModel model(nombrePortes,nombreDesti,nombreSorti);

    vector<int>    A(5);
    vector<double> B(5);


//A: indices de portes pour Destinations et Sorties
//    Desti:  10  20  25           Sortie: 0  29
    A[0]=10;    A[1]=20;    A[2]=25;    A[3]=0;     A[4]=29;


//B: proportions pour Destinations et Sorties
//    Desti: 1/3  1/3  1/3         Sortie: 1/2  1/2
//Ils vont être normalisés après.
    B[0]=1.0;   B[1]=1.0;   B[2]=1.0;   B[3]=1.0;   B[4]=1.0;






    vector<double> resultat;
    //gotominilocal NEED TO BE FIXED
    resultat=model.getRepartition(A,B);
    printOutResultat(resultat);
    return 0;

}
