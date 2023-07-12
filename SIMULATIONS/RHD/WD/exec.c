#include"main.h"
float eos(float valor, int j)
{
    float pressure, energy_den, mass_den;
    int i;

    if(j==PRE){//If we know the mass density and we want the pressure
        if (valor>=coreparam->rho0){
            if(valor>=coreparam->rho[1]) i=2;
            else if(valor>=coreparam->rho[0]) i=1;
            else i=0;
            print("\n%.8E\n", coreparam->rho[i]);
            energy_den = coreparam->k[i]*pow(valor,coreparam->gamma[i])/(coreparam->gamma[i]-1) + (1+coreparam->a[i])*valor - coreparam->lambda[i];
            pressure   = coreparam->k[i]*pow(valor,coreparam->gamma[i]) + coreparam->lambda[i];
        }
        else{
            if(valor>=crustparam->rho[4]) i=4;
            else if(valor>=crustparam->rho[3]) i=3;
            else if(valor>=crustparam->rho[2]) i=2;
            else if(valor>=crustparam->rho[1]) i=1;
            else if(valor>=crustparam->rho[0]) i=0;
            print("\n%.8E\n", crustparam->rho[i]);
            energy_den = crustparam->k[i]*pow(valor,crustparam->gamma[i])/(crustparam->gamma[i]-1) + (1+crustparam->a[i])*valor - crustparam->lambda[i];
            pressure   = crustparam->k[i]*pow(valor,crustparam->gamma[i]) + crustparam->lambda[i];
        }
        return pressure*pow(c_cgs,2);
    }

    if(j==RHO){//If we know the pressure and we want the mass density
        pressure = valor/pow(c_cgs,2);
        if (pressure>=coreparam->pre0){
            if(pressure>=coreparam->pre[1]) i=2;
            else if(pressure>=coreparam->pre[0]) i=1;
            else i=0;
            mass_den = pow((pressure - coreparam->lambda[i])/coreparam->k[i], 1/coreparam->gamma[i]);
        }
        else{
            if(pressure>=crustparam->pre[4]) i=4;
            else if(pressure>=crustparam->pre[3]) i=3;
            else if(pressure>=crustparam->pre[2]) i=2;
            else if(pressure>=crustparam->pre[1]) i=1;
            else if(pressure>=crustparam->pre[0]) i=0;
            mass_den = pow((pressure - crustparam->lambda[i])/crustparam->k[i], 1/crustparam->gamma[i]);
        }
        return mass_den;
    }
}