#include "Interp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Interpolación Polinómica
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

double horner(double z, int n, double *x, double *c) {
    // Començarem pel coeficient de major grau
    double resultat = c[n];

    // Iterarem desde i=n-1 fins a 0
    for (int i=n-1; i>=0; i--) {
        resultat = c[i] + (z - x[i]) * resultat;
    }

    return resultat;
}

int difdiv(int n, double *x, double *y) {
    // Tolerància
    const double TOLERANCIA = 1e-12; 

    // Recorrem totes les columnes (ordres) de la taula de diferencies dividides
    for (int i = 1; i <= n; i++) {
        
        // Recorrem tots els elements de la columna n i realitzem el càlcul de la diferencia dividida d'ordre i
        // Però anirem de n fins a i, per que sino operariem amb dades incorrectes ja que estem treballant i escrivint a la mateixa taula
        for (int j = n; j >= i; j--) {
            
            // Calculem el denominador (diferencia dels extrems)
            double denominador = x[j] - x[j - i];
            
            // Comprovem si el denominador es suficientment gran
            if (fabs(denominador) < TOLERANCIA) {
                return -1;
            }
            
            // Calculem la nova diferencia dividida
            // y[j] i y[j-1] son les diferencies d'ordre (i-1) de la iteració anterior
            y[j] = (y[j] - y[j - 1]) / denominador;
        }
        
        // D'aquesta manera, anirem emmagatzemant a cada element de y el valor que toca:
        // EXEMPLE: y[0] = f(x0) | y[1] = f[xo,x1] | y[2] = f[x0,x1,x2] | ... | y[n] = f[x0,x1,x2,...,xn]
    }
    
    return 0;
}

// Integración Numérica
double trap_full(int n, double a, double b) {
    // Definim h i calculem alguns termes del resultat
    double h = (b-a)/n;
    double resultat = 0.5 * f(a) + 0.5 * f(b);

    // Ara ens faltarà calcular el sumatori des de 1 fins a n-1 de f(xj)
    for (int j=1; j<=n-1; j++) {
        resultat += f(a + j*h);
    }

    // Multipliquem tot per h
    resultat *= h;

    return resultat;
}

double f(double x) {
    return exp(x);
}

double trap_refined(double Tn, int n, double a, double b) {
    // n ==> intervals per Tn
    // 2m ==> intervals per T2n
    
    // Calculem la nova h
    double h_nou = (b - a) / (2.0 * n); 

    // Sumem f(x) només en els nous punts (els imparells)
    double suma_nous = 0.0;
    
    // El bucle s'executa n vegades, una per cada nou punt afegit
    // K va de 0 a n-1, és a dir, els indexs imparells
    for (int k = 0; k < n; k++) {
        double x_nous = a + (2.0 * k + 1.0) * h_nou; 
        
        suma_nous += f(x_nous);
    }

    // Apliquem la formula de recurrencia
    // T2n = (1/2) * Tn + h_nou * (Sumatori dels nous f(x))
    double T2n = 0.5 * Tn + h_nou * suma_nous;
    
    return T2n;
}