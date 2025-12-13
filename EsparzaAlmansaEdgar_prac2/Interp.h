#ifndef INTERP_H
#define INTERP_H

// Interpolación de Lagrange
int lagrange_eval(double z, int n, double *x, double *y, double *pz);
// Método de Newton
double horner(double z, int n, double *x, double *c);
int difdiv(int n, double *x, double *y);
// Integración numérica
double trap_full(int n, double a, double b);
double trap_refined(double Tn, int n, double a, double b);
double f(double x);
double trap_refined(double Tn, int n, double a, double b);

#endif 