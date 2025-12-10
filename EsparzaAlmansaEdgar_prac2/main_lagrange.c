#include "Interp.h"

#include <stdio.h>
#include <stdlib.h>

double fun(double x);

double calculaX(int eqOcheb, int n, int j, double a, double b) {
    double pi = acos(-1.0); // Guardem el valor de pi (he volgut que sigui el més exacte possible)

    if (eqOcheb) { // Si escollim equidistant
        return a + j*((b-a)/n);
    } else { // Si escollim de Chebyshev
        return ((a+b)/2) + ((b-a)/2) * cos((2*j+1)*pi/(2*(n+1)));
    }
}

int main(void) {
    // Demanem totes les dades
    int n, eqOcheb; // Si eq0cheb = 1 EQUIDISTANTS | si = 0 CHEBYSHEV
    double a, b;

    printf("Introduce el grado del polinomio interpolador:\n");
    scanf("%d\n", &n);

    printf("Introduce los extremos a i b del intervalo [a,b]:\n");
    scanf("%lf %lf\n", &a, &b);

    printf("Elige entre nodos equidistantes (0) o de Chebyshev (1):\n");
    scanf("%d", &eqOcheb);

    // Reservem espai a memòria
    int N_NODOS = n + 1;

    double *x_j = (double *)malloc(N_NODOS * sizeof(double));
    double *y_j = (double *)malloc(N_NODOS * sizeof(double));

    if (x_j == NULL || y_j == NULL) {
        perror("Error de asignacion de memoria");
        return 1;
    }

    // Creem els n+1 nodes
    for (int j=0; j<n+1; j++) {
        double x = calculaX(eqOcheb, n, j, a, b);

        x_j[j] = x;
        y_j[j] = fun(x);
    }

    // Evaluem els polinomis

    return 0;
}