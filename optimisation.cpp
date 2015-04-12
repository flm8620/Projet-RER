#include "optimisation.h"
#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <cmath>
#include <iomanip>
using namespace std;
double Optimisation::fonctionObjectif(Variables u)
{
    vector<double> R;
    R=modelUtilise->getRepartitionNonNormalizePropo(u.desti,u.propoDesti,u.propoSorti,u.propoVoyageur);

    //norme l^2
    double l2=0;
    if(R.size()!=observation.size()){
        cout<<"error: size observation is wrong!"<<endl;
    }
    for(int i=0;i<R.size();i++){
        l2+=(R[i]-observation[i])*(R[i]-observation[i]);
    }
    //return sqrt(l2);
    return l2;
}

void Optimisation::projectionSousContrainte(Variables &u)
{
    //Contrainte:
    //1.
    //les indices de destination,  >=0 et <=nbPorte
    int nbPortes=modelUtilise->getNbPortes();
    for(int i=0;i<u.desti.size();i++){
        u.desti[i]=max(min(u.desti[i],nbPortes-1+0.5),-0.5);
    }

    //2.
    //propoDesti[1]+propoDesti[2]+... = 1
    //propoDesti[i]>0


    //propoDesti[i]>0
    for(int i=0;i<u.propoDesti.size();i++){

        u.propoDesti[i]=max(0.0,u.propoDesti[i]);
    }
    double S=0;
    //propoDesti[1]+propoDesti[2]+... = 1
    for(int i=0;i<u.propoDesti.size();i++){
        S+=u.propoDesti[i];
    }
    for(int i=0;i<u.propoDesti.size();i++){
        u.propoDesti[i]/=S;
    }

    //3.
    //propoSorti[1]+propoSorti[2]+... = 1
    //propoSorti[i]>0


    //propoSorti[i]>0
    for(int i=0;i<u.propoSorti.size();i++){

        u.propoSorti[i]=max(0.0,u.propoSorti[i]);
    }
    S=0;
    //propoSorti[1]+propoSorti[2]+... = 1
    for(int i=0;i<u.propoSorti.size();i++){
        S+=u.propoSorti[i];
    }
    for(int i=0;i<u.propoSorti.size();i++){
        u.propoSorti[i]/=S;
    }

    //4.
    //propoVoyageur[1]+propoVoyageur[2]+... = 1
    //propoVoyageur[i]>0

    //propoVoyageur[i]>0
    S=0;
    for(int i=0;i<u.propoVoyageur.size();i++){

        u.propoVoyageur[i]=max(0.0,u.propoVoyageur[i]);
    }
    S=0;
    //propoSorti[1]+propoSorti[2]+... = 1
    for(int i=0;i<u.propoVoyageur.size();i++){
        S+=u.propoVoyageur[i];
    }
    for(int i=0;i<u.propoVoyageur.size();i++){
        u.propoVoyageur[i]/=S;
    }
}





Variables Optimisation::calcGradient(Variables u)
{
    Variables gradU,uNext,uPrev;
    gradU=u;
    //gradient d'indice de destination
    double dx=0.01;
    for(int i=0; i<u.propoDesti.size();i++){
        uNext=u;
        uPrev=u;
        uNext.desti[i]+=dx;
        uPrev.desti[i]-=dx;
        // df/dx = ( f(x+dx)-f(x-dx) ) / 2dx
        gradU.desti[i]=( fonctionObjectif(uNext)-fonctionObjectif(uPrev) )/2/dx;
    }
    //gradient de proportion de destination
    dx=0.00001;
    for(int i=0; i<u.propoDesti.size();i++){
        uNext=u;
        uPrev=u;
        uNext.propoDesti[i]+=dx;
        uPrev.propoDesti[i]-=dx;
        // df/dx = ( f(x+dx)-f(x-dx) ) / 2dx
        gradU.propoDesti[i]=( fonctionObjectif(uNext)-fonctionObjectif(uPrev) )/2/dx;
    }

    //gradient de proportion de sortie
    for(int i=0; i<u.propoSorti.size();i++){
        uNext=u;
        uPrev=u;
        uNext.propoSorti[i]+=dx;
        uPrev.propoSorti[i]-=dx;
        // df/dx = ( f(x+dx)-f(x-dx) ) / 2dx
        gradU.propoSorti[i]=( fonctionObjectif(uNext)-fonctionObjectif(uPrev) )/2/dx;
    }
    //gradient de proportion de voyageur
    for(int i=0; i<u.propoVoyageur.size();i++){
        uNext=u;
        uPrev=u;
        uNext.propoVoyageur[i]+=dx;
        uPrev.propoVoyageur[i]-=dx;
        // df/dx = ( f(x+dx)-f(x-dx) ) / 2dx
        gradU.propoVoyageur[i]=( fonctionObjectif(uNext)-fonctionObjectif(uPrev) )/2/dx;
    }
    return gradU;
}




Variables Optimisation::unPas(Variables &u, Variables gradientU)
{
    Variables u_du=u;

    double lambda=0.0000001;
    double du;
    for(int i=0;i<gradientU.desti.size();i++){
        du=gradientU.desti[i]*lambda*10;
        u_du.desti[i]-=min(max(du,-0.01),0.01);
        if(abs(du)>0.01){cout<<"lambda for desti too large"<<endl;}
    }
    for(int i=0;i<gradientU.propoDesti.size();i++){
        du=gradientU.propoDesti[i]*lambda*0.5;
        u_du.propoDesti[i]-=min(max(du,-0.01),0.01);
        if(abs(du)>0.01){cout<<"lambda for propoDesti too large"<<endl;}
    }
    for(int i=0;i<gradientU.propoSorti.size();i++){
        du=gradientU.propoSorti[i]*lambda;
        u_du.propoSorti[i]-=min(max(du,-0.01),0.01);
        if(abs(du)>0.01){cout<<"lambda for propoSorti too large"<<endl;}
    }
    for(int i=0;i<gradientU.propoVoyageur.size();i++){
        du=gradientU.propoVoyageur[i]*lambda*0.5;
        u_du.propoVoyageur[i]-=min(max(du,-0.01),0.01);
        if(abs(du)>0.01){cout<<"lambda for propoVoyageur too large"<<endl;}
    }
    return u_du;
}

double Optimisation::TesterConvergence(Variables u_k, Variables u_k_1, Variables u0)
{
    return abs(fonctionObjectif(u_k)-fonctionObjectif(u_k_1))/fonctionObjectif(u0);
}

void Optimisation::printCompare(Variables u)
{
    vector<double> R=modelUtilise->getRepartitionNonNormalizePropo(u.desti,u.propoDesti,u.propoSorti,u.propoVoyageur);
    double S1=0,S2=0;
    cout<<"Comparer les résultats:"<<endl;
    cout<<"\t\tObservation\tSimulation\tDifférence"<<endl;
    for(int i=0;i<R.size();i++){
        cout<<"\t\t"<<observation[i]<<"\t\t"<<R[i]<<"\t\t"<<observation[i]-R[i]<<endl;
        S1+=observation[i];
        S2+=R[i];
    }
    cout<<endl;
    cout<<"somme\t\t"<<S1<<"\t\t"<<S2<<endl;
    cout<<"Norme l^2: || obser-simu ||_l^2 =\t"<<fonctionObjectif(u)<<endl;
}

Optimisation::Optimisation()
{
    modelUtilise=0;
}

Optimisation::~Optimisation()
{

}

void Optimisation::setObservation(std::vector<double> observation)
{
    if(modelUtilise==0){
        cout<<"Need to set Model in use firstly"<<endl;
        exit(0);
    }
    if(observation.size()!=modelUtilise->getNbPortes()){
        cout<<"Attention, the nb of doors in model and in observation is different!"<<endl;
        exit(0);
    }
    this->observation=observation;
}

void Optimisation::setModel(StationModel &model)
{
    modelUtilise=&model;
}

Variables Optimisation::minimiser(Variables uStart)
{
    Variables u0,u_k,u_k_1,u_solu,gradU;
    bool trouvee=false;
    u0=u_k=uStart;

    //seuil pour la condition d'arrêt
    double seuil=0.0000001;
    cout<<"u0 = "<<endl;
    u0.print();
    cout<<"J(u0) = "<<fonctionObjectif(u0)<<endl<<endl<<"Commencer l'optimisation"<<endl;
    printOutCompare(u0);
    //DEBUG
    //printCompare(u0);
    //Imagine::endGraphics();
    //exit(0);
    for(int i=0;i<20000;i++){
        //----output------
        u_k.print();
        cout<<"J(uk)="<<fonctionObjectif(u_k)<<endl;
        printOutCompare(u_k);

        gradU=calcGradient(u_k);
        u_k_1=unPas(u_k,gradU);
        projectionSousContrainte(u_k_1);

        //condition d'arrêt
        if(TesterConvergence(u_k,u_k_1,u0)<seuil){
            u_solu=u_k_1;
            trouvee=true;
            break;
        }
        u_k=u_k_1;
    }

    if(trouvee){
        cout<<"solution trouvée!"<<endl;
        cout<<"u = "<<endl;
        u_solu.print();
        cout<<endl;
        printCompare(u_solu);

    }else{
        cout<<"Je trouve pas la solution, désolé. la solution la plus proche:"<<endl;
        cout<<"u = "<<endl;
        u_k_1.print();
        cout<<endl;
        printCompare(u_k_1);

    }
    return u_solu;
}

void Optimisation::printOutCompare(Variables u)
{
    vector<double> R=modelUtilise->getRepartitionNonNormalizePropo(u.desti,u.propoDesti,u.propoSorti,u.propoVoyageur);
    Imagine::setBackGround(Imagine::WHITE);
    for(int i=0; i<R.size();i++){
        Imagine::fillRect(20 * i+10, 599.-R[i],5,R[i], Imagine::MAGENTA);
        Imagine::fillRect(20 * i+15, 599.-observation[i],5,observation[i], Imagine::BLACK);

    }

}



void Variables::print()
{
    //cout<<"Des ";
    for(int i=0;i<desti.size();i++){
        cout<<desti[i]<<",";
    }
    //cout<<"proD ";
    for(int i=0;i<propoDesti.size();i++){
        cout<<setiosflags(ios::showpoint)<<setprecision(4)<<propoDesti[i]<<",";
    }
    //cout<<"proS ";
    for(int i=0;i<propoSorti.size();i++){
        cout<<setiosflags(ios::showpoint)<<setprecision(4)<<propoSorti[i]<<",";
    }
    //cout<<"proV ";
    for(int i=0;i<propoVoyageur.size();i++){
        cout<<setiosflags(ios::showpoint)<<setprecision(4)<<propoVoyageur[i]<<",";
    }

}
