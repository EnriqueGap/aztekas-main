#include"main.h"
double function(double x, double y);
double f(int eq);
int ode_solver();

int test_module(){
    rk_order=2;
    rk_ rk;
    int j;
    rk.h=(x1max-x1min)/Nx1;
    for (int i=gc; i<=Nx1-gc; i++){
        j=i-gc;
        rk.u0=U(RHO,i);
        rk.f=function(grid.X1[i], rk.u0);
        Runge_Kutta(&rk, 1);
        rk.f=function(grid.X1[i]+rk.h, rk.u1);
        Runge_Kutta(&rk, 2);
        U(RHO,i+1)=rk.u2;
    }
}

double function(double x, double y)
{
    return cos(x);
}

int ode_solver()
{
    rk_order =2;
    rk_ rk;
    rk.h=(x1max-x1min)/Nx1;
    double aux0[eq], aux1[eq];
    for (int i =gc; i<=Nx1-gc; i++){
        for (int n=0; n<eq; n++){
            aux0[n]=U(n,i);
            rk.f=f(n); //esto está mal
            Runge_Kutta(&rk, 1);
            aux1[n]=rk.u1;
            rk.f=f(n); //también esto está mal
            Runge_Kutta(&rk, 2);
            U(n,i+1)=rk.u2;
        }
    }
}
double f(int eq)
{
    printf("Hello World");
}