#ifndef OPTIMISATION_H
#define OPTIMISATION_H
#include <vector>
#include <iostream>

#include "stationmodel.h"
struct Variables{
    std::vector<double>desti;
    std::vector<double>propoDesti;
    std::vector<double>propoSorti;
    std::vector<double>propoVoyageur;
    void print();
};

class Optimisation
{
private:
    std::vector<double> observation;
    StationModel* modelUtilise;
    double fonctionObjectif(Variables u);
    void projectionSousContrainte(Variables& u);
    Variables calcGradient(Variables u);
    Variables unPas(Variables &u,Variables gradientU);
    double TesterConvergence(Variables u_k,Variables u_k_1,Variables u0);
    void printCompare(Variables u);
public:
    Optimisation();
    ~Optimisation();
    void setObservation(std::vector<double>observation);
    void setModel(StationModel& model);
    Variables minimiser(Variables uStart);
    void printOutCompare(Variables u);
};

#endif // OPTIMISATION_H
