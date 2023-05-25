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

int buildModel();
double  interior_ode(double x, double y[], int i);
double  taylor_ode(double x, double y[], int i);

double eos(double valor, int i);

float massDensity, polytropicExp, polytropicK;