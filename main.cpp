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
int main(){
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

    uStart.desti.assign(3,0.0);
    uStart.desti[0]=randomAB(-0.5,nombrePortes-0.5);
    uStart.desti[1]=randomAB(-0.5,nombrePortes-0.5);
    uStart.desti[2]=randomAB(-0.5,nombrePortes-0.5);
    uStart.propoDesti.assign(3,0.0);
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


    openWindow(620, 700);

    //minimisation
    //regarder Page 69 du poly Calcul Scientifique
    optim.minimiser(uStart);

    endGraphics();
}

int main2(){


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
    uStart.propoVoyageur.assign(4,0);
    uStart.propoVoyageur[0]=0.5;    //loi normal
    uStart.propoVoyageur[1]=0.3;    //Confort gotominilocal
    uStart.propoVoyageur[2]=0.16;   //loi exp,retard stresee
    uStart.propoVoyageur[3]=0.05;   //loi uniform retard non-stress

    uStart.desti.assign(3,0.0);
    uStart.desti[0]=5.3;
    uStart.desti[1]=18;
    uStart.desti[2]=25;
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

int main3(){
    int nombrePortes=30;
    int nombreDesti=3;
    vector<double> observations=readFile("data.txt",nombrePortes);
    StationModel model(nombrePortes,nombreDesti);
    Optimisation optim;
    optim.setModel(model);
    optim.setObservation(observations);

    double a[12]={
        27.88,16.46,25.51,0.2713,0.7287,0.000,0.4269,0.5731,0.1059,0.1222,0.2971,0.4748
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

    u.desti.assign(3,0.0);
    u.desti[0]=a[0];
    u.desti[1]=a[1];
    u.desti[2]=a[2];
    u.propoDesti.assign(3,0.0);
    u.propoDesti[0]=a[3];
    u.propoDesti[1]=a[4];
    u.propoDesti[2]=a[5];
    u.propoSorti.assign(2,0.0);
    u.propoSorti[0]=a[6];
    u.propoSorti[1]=a[7];






    vector<double> resultat;
    resultat=model.getRepartitionNonNormalizePropo(u.desti,u.propoDesti,u.propoSorti,u.propoVoyageur);
    openWindow(620, 700);
    printOutResultat(resultat);
    optim.printOutCompare(u);
    //optim.minimiser(u);
    endGraphics();
    return 0;

}
