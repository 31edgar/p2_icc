#include "Interp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int lagrange_eval(double z, int n, double *x, double *y, double *pz) {
    // Creem una variable on anirem sumant cada f(xi) li(z) que calculem
    double sum = 0;

    // Comencem a calcular el polinomi de Lagrange
    for (int i=0; i<n+1; i++) {
        // Variable on inicialment guardarem Li(z)
        double producte = 1;

        // Ara, calcularem Li(z)
        for (int j=0; j<n+1; j++) {
            // Ens assegurem que j sigui diferent de i
            if (j != i) {
                // Calculem el denominador
                double denominador = (x[i]-x[j]);

                // Ens assegurem que sigui suficientment gran
                if (fabs(denominador) < 1e-12)
                    return -1;

                // Actualitzem el producte
                producte *= (z-x[j]) / denominador;
            }
        }

        // Una vegada tenim Li(z) guardat a producte, el multipliquem per f(xi)
        producte *= y[i];

        // I per últim, sumem el producte final a la suma de productes
        sum += producte;
    }

    // Una vegada hem acabat, guardarem el resultat (la suma) a pz
    *pz = sum;

    return 0;
}