#include"main.h"
float eos(float valor, int i)
{
    float mass_density;
    if(i==PRE){//Compute the pressure given a mass density
        return(polytropicK*pow(valor,polytropicExp));
    }
    if(i==RHO){//Compute the energy density given the pressure
        mass_density = pow(valor/polytropicK,1/polytropicExp);
        return(polytropicK*pow(mass_density,polytropicExp)/(polytropicExp-1) + mass_density*c_cgs*c_cgs);
    }
    if(i==RHOB){//Compute the energy density given the mass density
        return(polytropicK*pow(valor,polytropicExp)/(polytropicExp-1) + valor*c_cgs*c_cgs);
    }
}