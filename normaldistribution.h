#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H
#include <iostream>
#include <cmath>


class NormalDistribution
{
private:
    static double T[300];
    static bool isInitialized;
public:
    NormalDistribution();
    ~NormalDistribution();
    double getPhi(double x);
    double getPhi2(double x);
    double getPhi1(double x);
};

#endif // NORMALDISTRIBUTION_H
