#pragma once
#include <iostream>
static double _root3(double);
double root3(double);
int SolveP3(double*, double, double, double);
int solve(double, double, double, double);
double* getKoef(double**);

void input(double**);
void outMatr(double**, int);
double** generateMatrix(double(*test)[3]);
double** genTempMatr(double**, double);
void outEquals(double, double, double);
void outKorni(double*);
int gauss(double**, double*, double*, const int);