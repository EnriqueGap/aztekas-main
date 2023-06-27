#include "main.h"

float  interior_ode(float x, float y[], int i){
    
    float numerator, denominator;
    
    numerator = (4*PI*y[PRE]*pow(x,3)/pow(c_cgs,2)) + y[MB];
    denominator = x*((pow(c_cgs,2)*x/G_cgs) -2*y[MB]);
    
    if (i==RHO) return 0;                                           //no ode for E
    if (i==PRE) return(-(y[RHO]+y[PRE])*numerator/denominator);	    //dP/dr = -(E+P)(4pi*P*r^3/c^2 + m)/r(c^2*r/G -2m)
    if (i==MB) 	return(4*PI*y[RHO]*pow(x,2)/pow(c_cgs,2));	        //dm/dr = 4pi*E*r^2/c^2
}

float  taylor_ode(float x, float y[], int i){

    float numerator, denominator, central_density;
    
    central_density = eos(massDensity,RHOB);
    numerator = 4*PI*(3*y[PRE]*x + central_density*x);
    denominator = 3*pow(c_cgs,4)/G_cgs - (8*PI*central_density*pow(x,2));

    if (i==RHO) return 0;                                           //no ode for rho
    if (i==PRE) return(-(y[RHO]+y[PRE])*numerator/denominator);	    //dP/dr = -(E+P)(4pi)(3P*r + E_c*r)/3c^2*(c^2/G -8pi*E_c*r^2/3c^2))
    if (i==MB) 	return(4*PI*y[RHO]*pow(x,2)/pow(c_cgs,2));	        //dm/dr = 4pi*rho*r^2/c^2
}

