#include"main.h"

int buildModel(){
    rk_order=2;
    rk_ rk;
    int j,k;
    rk.h=(x1max-x1min)/(Nx1-2*gc);
    double y[VX1+1], f[VX1+1], u1[VX1+1], u2[VX1+1];

    for (int i=gc; i<=Nx1-gc; i++){
        //Mi contador "extra" en 0 y no en ghost cells
        j=i-gc;
        
        //Alternativa al grid
        //grid.X1[i]=x1min+rk.h*j

        //Que U empiece en i=gc es el primer valor en ser impreso en .dat
        //Es decir U(gc) es la condición inicial pues el valor impreso de x es x1min
        //y i-ésima
        for(k=0; k<=VX1; k++) y[k]=U(k,i);
        //Calculo de F con y i-ésima
        for(k=0; k<=VX1; k++) f[k]=interior_ode(grid.X1[i], y, k);
        //Para calcular y + rk.h*f hacemos
        for(k=0; k<=VX1; k++){
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
        for(k=0; k<=VX1; k++) f[k]=interior_ode(grid.X1[i]+rk.h, u1, k);
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
        if (u2[RHO]!=u2[RHO]){
            u2[RHO]=0;
            u2[PRE]=0;
            u2[MB]=U(MB,i);
        }
        //Paso a U i+1
        for(k=0; k<=VX1; k++) U(k,i+1)=u2[k];
    }
}
/*
#define RHO       0
#define PRE       1
#define VX1       2
#define VX2       3
#define VX3       4
*/