#include "optimisation.h"
#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <cmath>
#include <iomanip>
using namespace std;
double Optimisation::fonctionObjectif(Variables u)
{
    vector<double> R;
    R=modelUtilise->getRepartitionNonNormalizePropo(u.desti,u.propoDesti,u.propoSorti);

    //norme l^2
    double l2=0;
    if(R.size()!=observation.size()){
        cout<<"error: size observation is wrong!"<<endl;
    }
    for(int i=0;i<R.size();i++){
        l2+=(R[i]-observation[i])*(R[i]-observation[i]);
    }
    return l2;
}

void Optimisation::projectionSousContrainte(Variables &u)
{
    //Contrainte:
    //1.
    //les indices de destination,  >=0 et <=nbPorte
    //déjà réaliser dans calcGradient();

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

}

Variables Optimisation::takeBestMoveOfIndex(Variables u)
{
    Variables directionU=u,uNext,uPrev;

    // minimiserSurPropo(u):
    // On fix la partie indicePortes de u
    // on fait minimisation sur proportionDesti et proportion Sorti
    double J_min_u=minimiserSurPropo(u);

    //boucle sur les Destinations
    for(int i=0;i<u.desti.size();i++){
        uNext=u;
        uPrev=u;
        if(u.desti[i]==0){
            uNext.desti[i]+=1;
            if(minimiserSurPropo(uNext)<J_min_u){
                directionU.desti[i]=1;
            }else{
                directionU.desti[i]=0;
            }
        }else if(u.desti[i]==modelUtilise->getNbPortes()-1){
            uPrev.desti[i]-=1;
            if(minimiserSurPropo(uPrev)<J_min_u){
                directionU.desti[i]=-1;
            }else{
                directionU.desti[i]=0;
            }
        }else{
            uNext.desti[i]+=1;
            uPrev.desti[i]-=1;
            double left,right;
            left=minimiserSurPropo(uPrev);
            right=minimiserSurPropo(uNext);
            if(left<J_min_u){
                if(right<left){
                    directionU.desti[i]=1;
                }
                else{
                    directionU.desti[i]=-1;
                }
            }
            else if(right<J_min_u){
                directionU.desti[i]=1;
            }else{
                directionU.desti[i]=0;
            }
        }
    }

    //appliquer directionU sur u
    Variables u_du=u;
    for(int i=0;i<directionU.desti.size();i++){
        u_du.desti[i]+=directionU.desti[i];
    }
    minimiserSurPropo(u_du);//ici u_du change
    return u_du;

}

double Optimisation::minimiserSurPropo(Variables& u)
{
    Variables gradUPropo,u0,u_k,u_k_1,u_solu;
    bool trouvee=false;
    const int maxIter=1000;
    u0=u_k=u;
    //cout<<"sous pb starts:   Min sur Propo:"<<endl;
    //cout<<"J(u0) = "<<fonctionObjectif(u0)<<endl;
    double seuil=0.001;
    for(int i=0;i<maxIter;i++){
        //u_k.print();
        //cout<<"J(uk)="<<fonctionObjectif(u_k)<<endl;
        gradUPropo=calcGradientPourPropo(u_k);
        u_k_1=unPasPourProportions(u_k,gradUPropo);
        projectionSousContrainte(u_k_1);
        if(  TesterConvergence(u_k,u_k_1,u0) < seuil      ){
            u_solu=u_k_1;
            trouvee=true;
            break;
        }
        u_k=u_k_1;
    }
    if(trouvee){
        //cout<<"sous pb solution trouvée : "<<endl;
        //u_solu.print();cout<<endl;
        //cout<<"SOUS PROB ENDS"<<endl;
        u=u_solu;//assigner la valeur
        return fonctionObjectif(u_solu);//J(u_solu)
    }else{
        cout<<"sous pb échec: max itération atteint :"<<maxIter<<endl;
        cout<<"programme termine"<<endl;
        exit(0);
    }


}

Variables Optimisation::calcGradientPourPropo(Variables u)
{
    Variables gradU,uNext,uPrev;
    gradU=u;//pour avoir la meme taille avec u

    //gradient de proportion de destination
    double dx=0.00001;
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
    return gradU;


}

Variables Optimisation::unPasPourProportions(Variables &u,Variables gradientU)
{
    Variables u_du=u;

    double lambda=0.00000005;
    double du;
    for(int i=0;i<gradientU.propoDesti.size();i++){
        du=gradientU.propoDesti[i]*lambda;
        u_du.propoDesti[i]-=min(max(du,-0.01),0.01);
    }
    for(int i=0;i<gradientU.propoSorti.size();i++){
        du=gradientU.propoSorti[i]*lambda;
        u_du.propoSorti[i]-=min(max(du,-0.01),0.01);
    }
    return u_du;

}

Variables Optimisation::unPasPourIndiceDesti(Variables &u, Variables gradientU)
{
    Variables u_du=u;
    for(int i=0;i<gradientU.desti.size();i++){
        u_du.desti[i]+=gradientU.desti[i];
    }
    return u_du;
}

double Optimisation::TesterConvergence(Variables u_k, Variables u_k_1, Variables u0)
{
    return abs(fonctionObjectif(u_k)-fonctionObjectif(u_k_1))/fonctionObjectif(u0);
}

void Optimisation::printCompare(Variables u)
{
    vector<double> R=modelUtilise->getRepartition(u.desti,u.propoDesti,u.propoSorti);
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

Optimisation::setObservation(std::vector<double> observation)
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

Optimisation::setModel(StationModel &model)
{
    modelUtilise=&model;
}

Variables Optimisation::minimiser(Variables uStart)
{
    Variables u0,u_k,u_k_1,u_solu;
    bool trouvee=false;
    u0=u_k=uStart;

    //seuil pour la condition d'arrêt
    double seuil=0.001;
    cout<<"u0 = "<<endl;
    u0.print();
    cout<<"J(u0) = "<<fonctionObjectif(u0)<<endl<<endl<<"Commencer l'optimisation"<<endl;
    for(int i=0;i<2000;i++){

        //----output------
        u_k.print();
        cout<<"J(uk)="<<fonctionObjectif(u_k)<<endl;
        printOutCompare(u_k);

        //faire un pas sur les indiceDesti
        // y compris le gradient, le pas, le projection
        u_k_1=takeBestMoveOfIndex(u_k);

        //condition d'arrêt
        if(TesterConvergence(u_k,u_k_1,u_k)<seuil){
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
        cout<<"Je trouve pas la solution, désolé. Changer uStart?"<<endl;
        exit(0);
    }
    return u_solu;
}

void Optimisation::printOutCompare(Variables u)
{
    vector<double> R=modelUtilise->getRepartition(u.desti,u.propoDesti,u.propoSorti);
    Imagine::setBackGround(Imagine::WHITE);
    for(int i=0; i<R.size();i++){
        Imagine::drawLine(20 * i+10, 599, 20 * i+10, 599.-R[i], Imagine::MAGENTA);
        Imagine::drawLine(20 * i+12, 599, 20 * i+12, 599.-observation[i], Imagine::BLACK);
    }

}



void Variables::print()
{
    //cout<<"Des ";
    for(int i=0;i<desti.size();i++){
        cout<<desti[i]<<" ";
    }
    //cout<<"proD ";
    for(int i=0;i<propoDesti.size();i++){
        cout<<setiosflags(ios::showpoint)<<setprecision(4)<<propoDesti[i]<<"  ";
    }
    //cout<<"proS ";
    for(int i=0;i<propoSorti.size();i++){
        cout<<setiosflags(ios::showpoint)<<setprecision(4)<<propoSorti[i]<<"  ";
    }

}
