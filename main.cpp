#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include "stationmodel.h"
#include "optimisation.h"
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <time.h>
using namespace Imagine;
using namespace std;
double randomAB(double A=0.0,double B=1.0){
    return A+double(rand())/RAND_MAX*(B-A);
}

void printOutResultat(vector<double> R){

    double sum=0;
    for(int i=0; i<R.size();i++){
        drawLine(20 * i+10, 599, 20 * i+10, 599.-R[i], MAGENTA);
        cout<< R[i] << endl;
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
int main2(){
    srand (time(NULL));
    int nombrePortes=30;
    int nombreDesti=3;
    vector<double> observations=readFile("data.txt",nombrePortes);
    StationModel model(nombrePortes,nombreDesti);
    Optimisation optim;
    optim.setModel(model);
    optim.setObservation(observations);

    Variables uStart;
    uStart.propoVoyageur.assign(4,0);
    double r1=randomAB(),r2=randomAB(),r3=randomAB(),r4=randomAB();
    double S=r1+r2+r3+r4;
    r1/=S;r2/=S;r3/=S;r4/=S;
    uStart.propoVoyageur[0]=r1;    //loi normal
    uStart.propoVoyageur[1]=r2;    //Confort gotominilocal
    uStart.propoVoyageur[2]=r3;   //loi exp,retard stresee
    uStart.propoVoyageur[3]=r4;   //loi uniform retard non-stress

    uStart.desti.assign(nombreDesti,0.0);
    uStart.desti[0]=randomAB(-0.5,nombrePortes-0.5);
    uStart.desti[1]=randomAB(-0.5,nombrePortes-0.5);
    uStart.desti[2]=randomAB(-0.5,nombrePortes-0.5);
    uStart.propoDesti.assign(nombreDesti,0.0);
    r1=randomAB();r2=randomAB();r3=randomAB();
    S=r1+r2+r3;
    r1/=S;r2/=S;r3/=S;
    uStart.propoDesti[0]=r1;
    uStart.propoDesti[1]=r2;
    uStart.propoDesti[2]=r3;
    uStart.propoSorti.assign(2,0.0);
    r1=randomAB();r2=randomAB();
    S=r1+r2;
    r1/=S;r2/=S;
    uStart.propoSorti[0]=r1;
    uStart.propoSorti[1]=r2;

    uStart.sigma=randomAB(1,5);
    uStart.lambda=randomAB(0.5,1.5);


    openWindow(620, 700);

    //minimisation
    //regarder Page 69 du poly Calcul Scientifique
    optim.minimiser(uStart);

    endGraphics();
}

int main(){
    int nombrePortes=30;
    int nombreDesti=3;
    vector<double> observations=readFile("data.txt",nombrePortes);
    StationModel model(nombrePortes,nombreDesti);
    Optimisation optim;
    optim.setModel(model);
    optim.setObservation(observations);

    double a[14]={
        0,
        14.5,
        29,
        0.33,
        0.33,
        0.33,
        0.5,
        0.5,
        0.5,
        0,
        0.2,
        0.3,
        4.5,
        1.2,




    };
    Variables u;
    u.propoVoyageur.assign(4,0);
    u.propoVoyageur[0]=a[8]
           ;    //loi normal
    u.propoVoyageur[1]=a[9]
       ;    //Confort gotominilocal
    u.propoVoyageur[2]=a[10]
          ;   //loi exp,retard stresee
    u.propoVoyageur[3]=a[11]
          ;   //loi uniform retard non-stress

    u.desti.assign(nombreDesti,0.0);
    u.desti[0]=a[0];
    u.desti[1]=a[1];
    u.desti[2]=a[2];
    u.propoDesti.assign(nombreDesti,0.0);
    u.propoDesti[0]=a[3];
    u.propoDesti[1]=a[4];
    u.propoDesti[2]=a[5];
    u.propoSorti.assign(2,0.0);
    u.propoSorti[0]=a[6];
    u.propoSorti[1]=a[7];
    u.sigma=a[12];
    u.lambda=a[13];






    vector<double> resultat;
    resultat=model.getRepartitionNonNormalizePropo(u);
    openWindow(620, 700);
    printOutResultat(resultat);
    optim.printOutCompare(u);
    //optim.minimiser(u);
    endGraphics();
    return 0;

}
