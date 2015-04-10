#ifndef OPTIMISATION_H
#define OPTIMISATION_H
#include <vector>
#include <iostream>

#include "stationmodel.h"
struct Variables{
    std::vector<int>desti;
    std::vector<double>propoDesti;
    std::vector<double>propoSorti;
    void print();
};

class Optimisation
{
private:
    std::vector<double> observation;
    StationModel* modelUtilise;
    double fonctionObjectif(Variables u);
    void projectionSousContrainte(Variables& u);
    Variables takeBestMoveOfIndex(Variables u);
    double minimierSurPropo(Variables &u);
    Variables calcGradientPourPropo(Variables u);
    Variables unPasPourProportions(Variables &u, Variables gradientU);
    Variables unPasPourIndiceDesti(Variables &u, Variables gradientU);
    double TesterConvergence(Variables u_k,Variables u_k_1,Variables u0);
    void printCompare(Variables u);
public:
    Optimisation();
    ~Optimisation();
    setObservation(std::vector<double>observation);
    setModel(StationModel& model);
    Variables minimiser(Variables uStart);
    void printOutCompare(Variables u);
};

#endif // OPTIMISATION_H
