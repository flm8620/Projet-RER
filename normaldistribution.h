#ifndef NORMALDISTRIBUTION_H
#define NORMALDISTRIBUTION_H


class NormalDistribution
{
private:
    double T[300];
public:
    NormalDistribution();
    ~NormalDistribution();
    double getPhi(double x);
};

#endif // NORMALDISTRIBUTION_H
