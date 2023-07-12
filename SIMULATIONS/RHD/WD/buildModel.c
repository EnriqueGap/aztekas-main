#include"main.h"

int buildModel(){
    //gc = ghost cells
    for (int i=gc; i<=Nx1-gc; i++){

        if (i==gc){
            printf("\n\n\n%d\n\n\n",i);
            //When i==gc, radius=0
            //solve the ODE using Taylor Expansion for dP/dr
            ODE_Integration(i, 3, taylor_ode);
        }
        else{
            ODE_Integration(i, 3, interior_ode);
        }

        if(U(RHO,i+1)!=U(RHO,i+1)){
            printf("\n\n\n%d\n\n\n",i);
            system("read -p 'Press Enter to continue...' var");
            U(RHO,i+1)=0;
            U(PRE,i+1)=0;
            U(MB,i+1)=U(MB,i);
        }
    }
}