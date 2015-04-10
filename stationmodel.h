#ifndef STATIONMODEL_H
#define STATIONMODEL_H
#include <vector>
#include <cmath>

class StationModel
{
private:
    int nbPortes;
    int nbPrevoyants;
    int nbConfort;
    int tempsStation;
    double propStresse;
    int debitEntree;
    double ecartType;
    int nbDestinations;
    double sigma;
    //tableau contenant le nombre de passagers dans le train
    std::vector<double> s;
    std::vector<int> destinations;
    std::vector<double> proportionDestination;
    int nbSorties;
    std::vector<int> sorties;
    //proportion des passagers qui arrive par chaque sortie
    std::vector<double> proportionSorties;
    //paramètre de la loi exponentielle
    int lambda;

    void normalizeProportion();
    void loiNormal();
    void loiExponentiel();
    void loiUniforme();
    void gotominilocal();
    void initSortiesDesti(std::vector<int> indicesDesti,std::vector<double> proportionDesti, std::vector<double>proportionSorti);
public:
    StationModel(int nbPortes, int nbDestinations);
    ~StationModel();
    std::vector<double> getRepartition(std::vector<int> indicesDesti, std::vector<double> proportionDesti,std::vector<double> proportionSorti);
    std::vector<double> getRepartitionNonNormalizePropo(std::vector<int> indicesDesti, std::vector<double> proportionDesti,std::vector<double> proportionSorti);
    int getNbPortes(){return nbPortes;}
};

#endif // STATIONMODEL_H
