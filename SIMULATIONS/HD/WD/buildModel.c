#include"main.h"

int buildModel(){
    for (int i=gc; i<=Nx1-gc; i++){

        if (i==gc){
            ODE_Integration(i, 3, taylor_ode);
        }
        else{
            ODE_Integration(i, 3, interior_ode);
        }

        if(U(RHO,i+1)!=U(RHO,i+1)){
            U(RHO,i+1)=0;
            U(PRE,i+1)=0;
            U(MB,i+1)=U(MB,i);
        }
    }
}