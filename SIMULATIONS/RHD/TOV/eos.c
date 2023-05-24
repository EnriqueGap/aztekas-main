#include"main.h"
double eos(double valor, int i)
{
    //double Gamma=1.6666666666, k=3.119E12;
    //double Gamma=1.33333333, k=4.883E14;
    double Gamma=1.5, k=3.7E13;
    if(i==PRE){
        return k*pow(valor,Gamma);
    }
    if(i==RHO){
        return pow(valor/k,1/Gamma);
    }
}