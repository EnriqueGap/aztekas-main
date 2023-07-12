/*
 * File Name : initial.c
 * Description : aztekas initial module for Shock-Tube
 * Creation Date : 26-09-2019
 * Last Modified : 18-02-2020 09:43:01
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"main.h"

void Initial()
{
   //Initialize grid.time
   grid.time = 0.0;

   //Initialize dt
   dt = 0.0;

   makeCrust();
   makeCore();
   coreparam=&ms1b;                     //EOS to use in the Core
   crustparam=&sly4crust;              //EOS to use in the crust    

#if DIM == 1 

   /////////////////////////////
   //-------Riemann-1D--------//
   /////////////////////////////
   for(int i = 0; i <= Nx1; i++)
   {
      U(RHO,i) = massDensity;
      U(PRE,i) = eos(massDensity,PRE);
      U(MB,i) = 0;
   }
#endif
}
