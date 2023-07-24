/**
 * @file integration.h
 *
 * @author Alejandro Aguayo-Ortiz
 *
 * @brief Integration functions and variables
 *
 */

typedef struct
{
   double u0;
   double u1;
   double u2;
   double u3;
   double u4;
   double k1;
   double k2;
   double k3;
   double k4;
   double h;
   double f;
}rk_;

int rk_order;

/* Define pointers */                                                           
double *U, *U0, *U1, *U2, *U3;                                                  
double *Q, *Q0, *Q1, *Q2, *Q3;                                                  
                                                                                
double *U1p, *U1m;                                                              
double *U2p, *U2m;                                                              
                                                                                
/* Define number file */                                                        
int itprint;                                                                    
                                                                                
/* Define freq. output dt and time */                                           
double dtprint, tprint;  

void Equation_System_Solver();
 
void Hyperbolic_Integration();

void ODE_Integration(int i, int neq, double (*ODESystem)(float x, double y[], int i));
 
void Runge_Kutta(rk_ *rk, int order);
 
void Method_of_Lines(int order);
 
int MxV(double *M, double *V, double *L);
