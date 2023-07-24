#include"main.h"
double eos(double valor, int i)
{
    if(i==PRE){
        return polytropicK*pow(valor,polytropicExp);
    }
    if(i==RHO){
        float eos_mass_density = pow(valor/polytropicK,1/polytropicExp);
        return eos_mass_density*pow(c_cgs, 2) + valor/(polytropicExp -1);
    }
    if(i==RHOB){
        return valor*pow(c_cgs, 2) + polytropicK*pow(valor,polytropicExp)/(polytropicExp -1);
    }
}
/*float eos(float valor, int j)
{
    float pressure, energy_den, mass_den;
    int i;

    if(j==PRE){//If we know the mass density and we want the pressure
        if (valor>=coreparam->rho0){
            if(valor>=coreparam->rho[1]) i=2;
            else if(valor>=coreparam->rho[0]) i=1;
            else i=0;
            pressure   = coreparam->k[i]*pow(valor,coreparam->gamma[i]) + coreparam->lambda[i];
        }
        else{
            if(valor>=crustparam->rho[4]) i=4;
            else if(valor>=crustparam->rho[3]) i=3;
            else if(valor>=crustparam->rho[2]) i=2;
            else if(valor>=crustparam->rho[1]) i=1;
            else i=0;
            pressure   = crustparam->k[i]*pow(valor,crustparam->gamma[i]) + crustparam->lambda[i];
        }
        return pressure*pow(c_cgs,2);
    }
    if(j==RHOB){//If we know the mass density and we want the energy density
        if (valor>=coreparam->rho0){
            if(valor>=coreparam->rho[1]) i=2;
            else if(valor>=coreparam->rho[0]) i=1;
            else i=0;
            energy_den = coreparam->k[i]*pow(valor,coreparam->gamma[i])/(coreparam->gamma[i]-1) + (1+coreparam->a[i])*valor - coreparam->lambda[i];
        }
        else{
            if(valor>=crustparam->rho[4]) i=4;
            else if(valor>=crustparam->rho[3]) i=3;
            else if(valor>=crustparam->rho[2]) i=2;
            else if(valor>=crustparam->rho[1]) i=1;
            else i=0;
            energy_den = crustparam->k[i]*pow(valor,crustparam->gamma[i])/(crustparam->gamma[i]-1) + (1+crustparam->a[i])*valor - crustparam->lambda[i];
        }
        return energy_den*pow(c_cgs,2);
    }
    if(j==RHO){//If we know the pressure and we want the energy density
        pressure = valor/pow(c_cgs,2);
        if (pressure>=coreparam->pre0){
            if(pressure>=coreparam->pre[1]) i=2;
            else if(pressure>=coreparam->pre[0]) i=1;
            else i=0;
            mass_den = pow((pressure - coreparam->lambda[i])/coreparam->k[i], 1/coreparam->gamma[i]);
            energy_den = coreparam->k[i]*pow(mass_den,coreparam->gamma[i])/(coreparam->gamma[i]-1) + (1+coreparam->a[i])*mass_den - coreparam->lambda[i];
        }
        else{
            if(pressure>=crustparam->pre[4]) i=4;
            else if(pressure>=crustparam->pre[3]) i=3;
            else if(pressure>=crustparam->pre[2]) i=2;
            else if(pressure>=crustparam->pre[1]) i=1;
            else i=0;
            mass_den = pow((pressure - crustparam->lambda[i])/crustparam->k[i], 1/crustparam->gamma[i]);
            energy_den = crustparam->k[i]*pow(mass_den,crustparam->gamma[i])/(crustparam->gamma[i]-1) + (1+crustparam->a[i])*mass_den - crustparam->lambda[i];
        }
        return energy_den*pow(c_cgs,2);
    }
}*/
void makeCrust()
{
    sly4crust.rho[4]=5.317E11, sly4crust.rho[3]=3.350E11, sly4crust.rho[2]=1.826E8, sly4crust.rho[1]=6.285E5, sly4crust.rho[0]=0;

    sly4crust.k[4]=1.746E-8, sly4crust.k[3]=-7.957E29, sly4crust.k[2]=1.662E-6, sly4crust.k[1]=5.726E-8, sly4crust.k[0]=5.214E-9;

    sly4crust.gamma[4]=1.382, sly4crust.gamma[3]=-1.842, sly4crust.gamma[2]=1.269, sly4crust.gamma[1]=1.440, sly4crust.gamma[0]=1.611;

    sly4crust.lambda[4]=7.077E8, sly4crust.lambda[3]=1.193E9, sly4crust.lambda[2]=-6.025E3, sly4crust.lambda[1]=-1.354, sly4crust.lambda[0]=0;

    sly4crust.a[4]=8.208E-3, sly4crust.a[3]=1.035E-2, sly4crust.a[2]=-5.278E-4, sly4crust.a[1]=-1.861E-5, sly4crust.a[0]=0;

    for (int i=0; i<5; i++) sly4crust.pre[i]=sly4crust.k[i]*pow(sly4crust.rho[i],sly4crust.gamma[i]) + sly4crust.lambda[i];
}
void makeCore()
{
    apr.rho0=pow(10,14.040),	apr.k[0]=pow(10,-33.210),	apr.gamma[0]=3.169,	apr.gamma[1]=3.452,	apr.gamma[2]=3.310;
    completeCore(&apr);
    
    bhf.rho0=pow(10,14.130),	bhf.k[0]=pow(10,-35.016),	bhf.gamma[0]=3.284,	bhf.gamma[1]=2.774,	bhf.gamma[2]=2.618;
    completeCore(&bhf);
    
    fps.rho0=pow(10,14.087),	fps.k[0]=pow(10,-32.985),	fps.gamma[0]=3.147,	fps.gamma[1]=2.652,	fps.gamma[2]=2.199;
    completeCore(&fps);
    
    h4.rho0=pow(10,13.499),	h4.k[0]=pow(10,-23.310),	h4.gamma[0]=2.514,	h4.gamma[1]=2.333,	h4.gamma[2]=1.562;
    completeCore(&h4);
    
    kde0v.rho0=pow(10,13.978),	kde0v.k[0]=pow(10,-30.250),	kde0v.gamma[0]=2.967,	kde0v.gamma[1]=2.835,	kde0v.gamma[2]=2.803;
    completeCore(&kde0v);
    
    kde0v1.rho0=pow(10,13.929),	kde0v1.k[0]=pow(10,-29.232),	kde0v1.gamma[0]=2.900,	kde0v1.gamma[1]=2.809,	kde0v1.gamma[2]=2.747;
    completeCore(&kde0v1);
    
    mpa1.rho0=pow(10,14.088),	mpa1.k[0]=pow(10,-40.301),	mpa1.gamma[0]=3.662,	mpa1.gamma[1]=3.057,	mpa1.gamma[2]=2.298;
    completeCore(&mpa1);
    
    ms1.rho0=pow(10,13.657),	ms1.k[0]=pow(10,-30.170),	ms1.gamma[0]=2.998,	ms1.gamma[1]=2.123,	ms1.gamma[2]=1.955;
    completeCore(&ms1);
    
    ms1b.rho0=pow(10,13.795),	ms1b.k[0]=pow(10,-33.774),	ms1b.gamma[0]=3.241,	ms1b.gamma[1]=2.136,	ms1b.gamma[2]=1.963;
    completeCore(&ms1b);
    
    qhc19.rho0=pow(10,14.102),	qhc19.k[0]=pow(10,-36.879),	qhc19.gamma[0]=3.419,	qhc19.gamma[1]=2.760,	qhc19.gamma[2]=2.017;
    completeCore(&qhc19);
    
    rsCore.rho0=pow(10,13.641), rsCore.k[0]=pow(10,-25.150),	rsCore.gamma[0]=2.636,	rsCore.gamma[1]=2.677,	rsCore.gamma[2]=2.647;
    completeCore(&rsCore);
    
    sk255.rho0=pow(10,13.679),	sk255.k[0]=pow(10,-25.990),	sk255.gamma[0]=2.693,	sk255.gamma[1]=2.729,	sk255.gamma[2]=2.667;
    completeCore(&sk255);
    
    sk272.rho0=pow(10,13.732),	sk272.k[0]=pow(10,-27.597),	sk272.gamma[0]=2.804,	sk272.gamma[1]=2.793,	sk272.gamma[2]=2.733;
    completeCore(&sk272);
    
    ski2.rho0=pow(10,13.552),	ski2.k[0]=pow(10,-24.202),	ski2.gamma[0]=2.575,	ski2.gamma[1]=2.639,	ski2.gamma[2]=2.656;
    completeCore(&ski2);
    
    ski3.rho0=pow(10,13.660),	ski3.k[0]=pow(10,-26.457),	ski3.gamma[0]=2.729,	ski3.gamma[1]=2.680,	ski3.gamma[2]=2.708;
    completeCore(&ski3);
    
    ski4.rho0=pow(10,13.907),	ski4.k[0]=pow(10,-31.008),	ski4.gamma[0]=3.029,	ski4.gamma[1]=2.759,	ski4.gamma[2]=2.651;
    completeCore(&ski4);
    
    ski5.rho0=pow(10,13.438),	ski5.k[0]=pow(10,-23.109),	ski5.gamma[0]=2.505,	ski5.gamma[1]=2.708,	ski5.gamma[2]=2.727;
    completeCore(&ski5);
    
    ski6.rho0=pow(10,13.902),	ski6.k[0]=pow(10,-31.089),	ski6.gamma[0]=3.035,	ski6.gamma[1]=2.762,	ski6.gamma[2]=2.653;
    completeCore(&ski6);
    
    skmp.rho0=pow(10,13.763),	skmp.k[0]=pow(10,-27.116),	skmp.gamma[0]=2.766,	skmp.gamma[1]=2.741,	skmp.gamma[2]=2.698;
    completeCore(&skmp);
    
    skop.rho0=pow(10,13.761),	skop.k[0]=pow(10,-26.089),	skop.gamma[0]=2.693,	skop.gamma[1]=2.660,	skop.gamma[2]=2.579;
    completeCore(&skop);
    
    sly2.rho0=pow(10,13.967),	sly2.k[0]=pow(10,-31.070),	sly2.gamma[0]=3.026,	sly2.gamma[1]=2.871,	sly2.gamma[2]=2.760;
    completeCore(&sly2);
    
    sly230a.rho0=pow(10,14.021),    sly230a.k[0]=pow(10,-33.385),	sly230a.gamma[0]=3.184,	sly230a.gamma[1]=2.895,	sly230a.gamma[2]=2.588;
    completeCore(&sly230a);
    
    sly4.rho0=pow(10,13.980),   sly4.k[0]=pow(10,-31.350),	sly4.gamma[0]=3.045,	sly4.gamma[1]=2.884,	sly4.gamma[2]=2.773;
    completeCore(&sly4);
    
    sly9.rho0=pow(10,13.899),   sly9.k[0]=pow(10,-30.657),	sly9.gamma[0]=3.005,	sly9.gamma[1]=2.796,	sly9.gamma[2]=2.652;
    completeCore(&sly9);
    
    wff1.rho0=pow(10,14.133),	wff1.k[0]=pow(10,-34.394),	wff1.gamma[0]=3.240,	wff1.gamma[1]=3.484,	wff1.gamma[2]=3.695;
    completeCore(&wff1);
}

void completeCore(core *eos)
{
    eos->rho[2]=massDensity, eos->rho[1]=pow(10,14.99), eos->rho[0]=pow(10,14.87);

    eos->lambda[0]=sly4crust.k[4]*pow(eos->rho0,sly4crust.gamma[4])+sly4crust.lambda[4]-eos->k[0]*pow(eos->rho0,eos->gamma[0]);
    double aux=sly4crust.k[4]*pow(eos->rho0,sly4crust.gamma[4])/(sly4crust.gamma[4]-1) + (1+sly4crust.a[4])*eos->rho0 - sly4crust.lambda[4];
    eos->a[0]=aux/eos->rho0 - (eos->k[0]*pow(eos->rho0,eos->gamma[0]-1)/(eos->gamma[0]-1)) + (eos->lambda[0]/eos->rho0) -1;

    for (int i=0; i<2; i++){
        eos->k[i+1]=eos->k[i]*(eos->gamma[i]/eos->gamma[i+1])*pow(eos->rho[i],eos->gamma[i]-eos->gamma[i+1]);
        eos->lambda[i+1]=eos->lambda[i] + (1 - (eos->gamma[i]/eos->gamma[i+1]))*eos->k[i]*pow(eos->rho[i],eos->gamma[i]);
        eos->a[i+1]=eos->a[i] + (eos->gamma[i]*eos->k[i]*pow(eos->rho[i],eos->gamma[i]-1)*(eos->gamma[i+1] - eos->gamma[i])/( (eos->gamma[i+1]-1)*(eos->gamma[i]-1) ));
    }
    for (int i=0; i<3; i++) eos->pre[i]=eos->k[i]*pow(eos->rho[i],eos->gamma[i]) + eos->lambda[i];
    eos->pre0=eos->k[0]*pow(eos->rho0,eos->gamma[0]) + eos->lambda[0];
}