/*
 * File Name : user_param.h
 * Description : aztekas user parameters header file for Relativistic Shock-Tube
 * Creation Date : 28-09-2019
 * Last Modified : 19-02-2020 13:20:54
 * Created By : Alejandro Aguayo-Ortiz
 */

#include"macros.h"

#define HORIZONTAL         0
#define VERTICAL           1
#define DIAGONAL           2

#define outflow_x1max      TRUE
#define outflow_x1min      TRUE
#define outflow_x2max      TRUE
#define outflow_x2min      TRUE

#define RECONST            MINMOD
#define FLUX               HLL
#define GRID               UNIFORM

#define INTERFACE          DIAGONAL

#define PI 3.141592653589
#define MB 2 //barionic mass
#define RHOB 3

int buildModel();
double  interior_ode(float x, double y[], int i);
double  taylor_ode(float x, double y[], int i);

double eos(double valor, int i);

float massDensity, polytropicExp, polytropicK;

typedef struct{
	double rho[5];
	double pre[5];
	double k[5];
	double gamma[5];
	double lambda[5];
	double a[5];
}crust;

crust sly4crust, *crustparam;

typedef struct{
	double rho[3];
	double pre[3];
	double k[3];
	double gamma[3];
	double lambda[3];
	double a[3];
	double rho0;
	double pre0;
}core;

core apr, bhf, fps, h4, kde0v, *coreparam,
	kde0v1, mpa1, ms1, ms1b, qhc19,
	rsCore, sk255, sk272, ski2, ski3,
	ski4, ski5, ski6, skmp, skop,
	sly2, sly230a, sly4, sly9, wff1;

void makeCore();
void makeCrust();
void completeCore(core *eos);