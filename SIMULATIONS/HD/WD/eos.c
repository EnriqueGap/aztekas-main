#include"main.h"
double eos(double valor, int i)
{
    if(i==PRE){
        return polytropicK*pow(valor,polytropicExp);
    }
    if(i==RHO){
        return pow(valor/polytropicK,1/polytropicExp);
    }
}