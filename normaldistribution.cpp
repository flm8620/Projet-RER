#include "normaldistribution.h"



NormalDistribution::NormalDistribution()
{

}

NormalDistribution::~NormalDistribution()
{

}


double NormalDistribution::getPhi(double x)
{
    //http://www.johndcook.com/blog/cpp_phi/

    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign;
    if (x < 0){
        sign = -1;
        x = -x/1.414213562373;//sqrt(2)
    }else{
        sign = 1;
        x = x/1.414213562373;//sqrt(2)
    }

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*std::exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

