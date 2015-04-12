#ifndef STATIONMODEL_H
#define STATIONMODEL_H
#include <vector>
#include <cmath>

class StationModel
{
private:
    int nbPortes;
    double nbPrevoyants;
    double nbConfort;
    double nbRetardStresse;
    double nbRetardNonStresse;
    int nbVoyageurTotal;
    double ecartType;
    int nbDestinations;
    double sigma;
    //tableau contenant le nombre de passagers dans le train
    std::vector<double> s;
    std::vector<double> destinations;
    std::vector<double> proportionDestination;
    int nbSorties;
    std::vector<double> sorties;
    //proportion des passagers qui arrive par chaque sortie
    std::vector<double> proportionSorties;
    //paramètre de la loi exponentielle
    double lambda;

    void normalizeProportion();
    void loiNormal();
    void loiExponentiel();
    void loiUniforme();
    void gotominilocal();
    void initSortiesDesti(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double>proportionSorti);
public:
    StationModel(int nbPortes, int nbDestinations);
    ~StationModel();
    std::vector<double> getRepartition(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti,std::vector<double> propoVoyageur);
    std::vector<double> getRepartitionNonNormalizePropo(std::vector<double> indicesDesti, std::vector<double> proportionDesti, std::vector<double> proportionSorti, std::vector<double> propoVoyageur);
    int getNbPortes(){return nbPortes;}
};

#endif // STATIONMODEL_H
