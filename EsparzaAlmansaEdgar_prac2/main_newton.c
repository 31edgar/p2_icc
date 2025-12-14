#include "Interp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    scanf("%d", &n);

    printf("Introduce los extremos a i b del intervalo [a,b]:\n");
    scanf("%lf %lf", &a, &b);

    printf("Elige entre nodos equidistantes (1) o de Chebyshev (0):\n");
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

    // Obrim l'arxiu de sortida
    FILE *archivo_salida = fopen("resultados_newton.txt", "w");
    if (archivo_salida == NULL) {
        perror("Error al abrir el archivo de salida 'resultados_newton.txt'");
        free(x_j);
        free(y_j);
        return -1;
    }

    // Hem d'evaluar en 1000 punts z_j a dins l'interval [a,b]
    const int NUM_PUNTOS_EVAL = 1000;

    // Per tant, calculem quina distància hi haurà entre cada punt consecutiu
    double h_eval = (b - a) / (NUM_PUNTOS_EVAL - 1.0);

    // Booleà per saber si hem pogut calcular el polinomi o no
    int error_eval;
    error_eval = difdiv(n, x_j, y_j);

    if (error_eval == -1) {
        // Alliberem la memoria
        printf("Error al executar la funció difdiv");
        fclose(archivo_salida);
        free(x_j);
        free(y_j);
        return -1;
    }

    // Bucle d'evaluació
    for (int j=0; j<NUM_PUNTOS_EVAL; j++) {
        // Calculem el punt z_j
        double z = a + j * h_eval;

        // Obtenim els valors necessaris
        double f_zj = fun(z); // Valor real de la funció
        double Pn_zj;   // Valor del polinomi interpolador
        
        // Calculem i guardem el valor del polinomi en el punt z_j
        Pn_zj = horner(z, n, x_j, y_j);

        double error_abs = fabs(f_zj - Pn_zj);

        fprintf(archivo_salida, "%.15e\t%.15e\t%.15e\t%.15e\n", z, f_zj, Pn_zj, error_abs);
    }

    // FINAL
    printf("Evaluacion completada. Resultados guardados en 'resultados_newton.txt'.\n");
    fclose(archivo_salida);
    free(x_j);
    free(y_j);
    return 0;
}