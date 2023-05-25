/**
 * @file /integration/ODE/integration.c
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Main function for the time integration in the conservative variables
 * \f$ \mathbf{Q} \f$.
 */

#include"main.h"

void ODE_Integration(int i, int neq, float (*ODESystem)(float x, float y[], int i))
{
   rk_order=2;
   rk_ rk;
   rk.h=(x1max-x1min)/(Nx1-2*gc);
   float y[neq], f[neq], u1[neq], u2[neq];
   int k;
   for(k=0; k<neq; k++) y[k]=U(k,i);
   //Calculo de F con y i-ésima
   for(k=0; k<neq; k++) f[k]=(*ODESystem)(grid.X1[i], y, k);
   //Para calcular y + rk.h*f hacemos
   for(k=0; k<neq; k++){
      //y i-ésima
      rk.u0=y[k];
      //F(y_i)
      rk.f=f[k];
      //Llamado para calcular y_i + rk.h*f
      Runge_Kutta(&rk, 1);
      //Guardado
      u1[k]=rk.u1;
   }
   u1[RHO]=eos(u1[PRE], RHO);
   //F(y+ rk.h+f)
   for(k=0; k<neq; k++) f[k]=(*ODESystem)(grid.X1[i]+rk.h, u1, k);
   //Calculo del paso final
   for(k=0; k<=VX1; k++){
      //y-iésima
      rk.u0=y[k];
      //y + rk.h*f
      rk.u1=u1[k];
      //F(y+ rk.h+f)
      rk.f=f[k];
      //Llamado para calcular y_i+1
      Runge_Kutta(&rk, 2);
      //Guardado
      u2[k]=rk.u2;
   }
   u2[RHO]=eos(u2[PRE], RHO);
   //Paso a U i+1
   for(k=0; k<=VX1; k++) U(k,i+1)=u2[k];
}
