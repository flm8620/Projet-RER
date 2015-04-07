#include <Imagine/Images.h>
#include <Imagine/Graphics.h>
#include <cmath>
#include "normaldistribution.h"
using namespace Imagine;
using namespace std;


const int nbportes = 30;
const int nbprevoyants = 1000;
const int nbconfort = 500;
const int tempsstation = 30;
const double propstresse = 0.8;
const int debitentree = 50;
const double ecarttype = 1;
const int nbdestinations = 3;
const double sigma = 5.0;
//tableau contenant le nombre de passagers dans le train
double s[30];
int Destination[nbdestinations];
double proportiondestination[nbdestinations];
const int nbsorties=2;
int sorties[nbsorties];
//proportion des passagers qui arrive par chaque sortie
double proportionsorties[nbsorties];
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
	for (int i = 0; i < nbsorties; i++){
		cout << "pour la sortie" << i+1 << " (entier entre 1 et 30)" << endl;
		cin >> sorties[i];
		cout << "un double quelconque" << endl;
		cin >> proportionsorties[i];
		S += proportionsorties[i];
	}
	for (int i = 0; i < nbsorties; i++){
		proportionsorties[i] = proportionsorties[i] / S;
	}
}
void init_Destination(){
	cout << "Designez la porte la plus pres de chaque destination favorite" << endl;
	cout << "Donnez la proportion de gens qui veulent y aller" << endl;
	for (int i = 0; i < nbdestinations; i++){
        cout << "pour la porte" << i + 1 << " (un entier entre 0 et 29)" << endl;
		cin >> Destination[i];
		cout << "un double quelconque" << endl;
		cin >> proportiondestination[i];
	}
}

void gotominlocal(){
	for (int i = 0; i < nbsorties; i++){
		for (int j = 0; j < proportionsorties[i] * nbconfort; j++){
			if (sorties[i] == 1){
				for (int k = 1; k < nbportes; k++){
					if (s[k + 1] >= s[k]){ s[k] += 1; }
				}
			}
			if (sorties[i] == nbportes){
				for (int k = nbportes-1; k > 0; k--){
					if (s[k - 1] >= s[k]){ s[k] += 1; }
				}
			}
			if (sorties[i] < nbportes && sorties[i]>0){
                for (int k = sorties[i]; k < max( nbportes - sorties[i], sorties[i] ); k++){
					if (s[sorties[i] + k + 1] >= s[sorties[i] + k]){ s[sorties[i] + k] += 1; }
					else if (s[sorties[i] - k - 1] >= s[sorties[i] - k]){ s[sorties[i] - k] += 1; }
				}
			}
		}
	}
}

void loiexponentiel(int sorties[nbsorties]){
	//boucle sur toutes les sorties (1 loi par sortie)
	for (int i = 0; i < nbsorties; i++){
		//Cas sortie au début du quai
		if (sorties[i] == 0){
			for (int j = 0; j < nbportes; j++){
                double A=(1 - exp(-lambda)) / (1 - exp(-lambda*nbportes));
                //A est une constante pour normer la loi
                s[j] += A*exp(-lambda*j)   *   tempsstation*debitentree*propstresse*proportionsorties[i];
				//la loi est normée puis elle est multpliée par le nombre de passager stressé et en retard entrant par la sortie i.
            }
		}
		//Cas sortie à la fin du quai
		if (sorties[i] == nbportes - 1){
			for (int j = 0; j < nbportes; j++){
                double A=
                s[j] += (1 - exp(-lambda)) / (1 - exp(-lambda*nbportes))*exp(-lambda*(nbportes - j-1))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
		}
		//Cas sortie au milieu du quai
        if (sorties[i]>0 && sorties[i] < nbportes-1){
			for (int j = 0; j < sorties[i]; j++){
                double A=(1 - exp(-lambda)) / (1 + exp(1) - exp(-lambda*sorties[i]) - exp(-lambda*(nbportes - sorties[i])));
                s[j] += A*exp(-lambda*(sorties[i] - 1 - j))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
			for (int j = sorties[i]; j < nbportes; j++){
                double A=(1 - exp(-lambda)) / (2 - exp(-lambda*sorties[i]) - exp(-lambda*(nbportes - sorties[i])));
                s[j] += A*exp(-lambda*(j + 1 - sorties[i]))*tempsstation*debitentree*propstresse*proportionsorties[i];
			}
			//cette expression est la synthèse de trois étapes
            //1/ à la porte la plus proche de la sortie, la fonction a pour valeur 1 et elle décroit de chaque côté selon une même loi exponentielle de paramètre -lambda
			//2/ je norme cette fonction sur le quai
			//3/ je multiplie par le nombre de passager stressé en retard entrant par la sortie i.
		}
	}
}

void loiuniforme(int sorties[nbsorties]){
	//boucle sur toutes les sorties
	for (int i = 0; i < nbsorties; i++){
		//boucle sur les destinations
		for (int j = 0; j < nbdestinations; j++){
			//boucle sur les portes entre la sortie et les destinations
			for (int k = min(i, j); k < max(i, j); k++){
                s[k] += tempsstation*debitentree*(1 - propstresse)*proportionsorties[i] * proportiondestination[j]/(j-i);
			}
		}
	}
}

void loinormale(int i, int j){
    //i - No.destination
    //j - No.porte
    NormalDistribution normal;
    double xSup=(j+1-Destination[i])/sigma;
    double xInf=(j - Destination[i])/sigma;
    double value=nbprevoyants*proportiondestination[i] * (normal.getPhi(xSup) - normal.getPhi(xInf));
    s[j] += value;
}


void init_quai(){
	for (int i = 0; i < nbportes; i++){
		s[i] = 0;
	}
}


int main(){
	init_quai();

	init_Destination();
	init_sorties();

	for (int i = 0; i < nbdestinations; i++){
		for (int j = 0; j < nbportes; j++){
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
	loiuniforme(sorties);


    openWindow(620, 600);
        for (int i = 0; i < nbportes; i++){
            drawLine(20 * i+10, 599, 20 * i+10, 599.-s[i], MAGENTA);
            cout << s[i] << endl;
        }
	system("pause");
	return 0;
}
