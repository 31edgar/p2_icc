#include "Interp.h" // Asumimos que contiene la definicion de fun(x)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    // Dades
    const int N_MAX_ITER = 10; // Nombre màxim d'iteracions, n=2^N_MAX_ITER. Si N_MAX=1024, N_MAX_ITER=10
    const double TOLERANCIA = 1e-6;

    const double A = 0.0;
    const double B = 1.0;
    const double I_EXACTE = exp(1.0) - 1.0;

    // Preparem l'arxiu de sortida
    FILE *archivo_salida = fopen("trapez_2n.dat", "w");
    if (archivo_salida == NULL) {
        perror("Error al abrir trapez_2n.dat");
        return 1;
    }

    // Variables on guardarem la informació
    double Tn = 0.0; // T_n
    double T2n = 0.0; // T_{2n}
    double dif_abs = 0.0;
    double error_abs = 0.0;
    int n = 1; // N de subintervals al pas anterior
    char *razon_finalizacion = "Límite de iteraciones (n) superado.";

    // Pas base (n=1)
    // T1 = (h/2) * (f(a) + f(b)), on h = (b-a)/1
    double h_base = B - A;
    Tn = (h_base / 2.0) * (f(A) + f(B)); // T_1

    // Capçalera de l'arxiu de sortida
    fprintf(archivo_salida, "2n\tT2n\t|T2n - Tn|\t|I - T2n|\n");

    // Bucle principal: s'itera per calcular T2, T4, T8, ...
    for (int i = 1; i <= N_MAX_ITER; i++) {
        // Calculem la nova aproximació T_{2n} a partirde T_n
        T2n = trap_refined(Tn, n, A, B);

        // Calculem la diferencia i error absolut
        dif_abs = fabs(T2n - Tn);
        error_abs = fabs(I_EXACTE - T2n);

        // Escrivim al arxiu de sortida (el nou n es 2*n)
        fprintf(archivo_salida, "%d\t%.15e\t%.15e\t%.15e\n", 
                2 * n, T2n, dif_abs, error_abs);
        
        // Criteri de parada
        if (dif_abs < TOLERANCIA) {
            razon_finalizacion = "Criterio de convergencia (|T2n - Tn| < 1e-6) alcanzado.";
            break; 
        }

        // Preparem la següent iteració
        Tn = T2n; // T_n actual pasa a ser T_n anterior
        n *= 2;   // n es duplica per a la següent crida a trap_refined
    }
    
    // Si el bucle acaba per iteracions, n=1024, T2n ja es l'ultim valor calculat
    fclose(archivo_salida);

    printf("\n----------------------------------------------------\n");
    printf("Proceso de integración numérica (Trapecio Refinado) finalizado.\n");
    printf("Razón de finalización: %s\n", razon_finalizacion);
    printf("----------------------------------------------------\n");
    printf("Integral: I = ∫_%.1f^%.1f e^x dx = e - 1\n", A, B);
    printf("Valor Exacto (I):\t\t%.15e\n", I_EXACTE);
    printf("Valor Aproximado (T2n):\t%.15e (con 2n=%d subintervalos)\n", T2n, 2 * n);
    printf("Error Absoluto Final:\t%.15e\n", error_abs);
    printf("----------------------------------------------------\n");

    return 0;
}