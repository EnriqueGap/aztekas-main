#include "main.h"

double  interior_ode(float x, double y[], int i){
    if (i==RHO) return 0;                               //no ode for rho
    if (i==PRE) return(-y[RHO]*G_cgs*y[MB]/pow(x,2));	//dP/dr = - rho*G*m/r^2
    if (i==MB) 	return(4*PI*y[RHO]*pow(x,2));	        //dm/dr = 4pi*rho*r^2
}

double  taylor_ode(float x, double y[], int i){
    if (i==RHO) return 0;                                       //no ode for rho
    if (i==PRE) return(-y[RHO]*G_cgs*4*PI*massDensity*x/3);	//dP/dr ~ - rho*G*4pi*rho_c*r/3
    if (i==MB) 	return(4*PI*y[RHO]*pow(x,2));	                //dm/dr = 4pi*rho*r^2
}