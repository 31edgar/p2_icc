#include "Interp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    // Dades
    const int N_MAX = 1024;
    const double TOLERANCIA = 1e-6;

    const double A = 0.0;
    const double B = 1.0;
    const double I_EXACTE = exp(1.0) - 1.0;

    // Preparem l'arxiu de sortida
    FILE *archivo_salida = fopen("trapez_naive.dat", "w");
    if (archivo_salida == NULL) {
        perror("Error al abrir trapez_naive.dat");
        return 1;
    }

    // Variables on guardarem la informació
    double Tn_anterior = 0.0;
    double Tn_actual = 0.0;
    double diferencia_sucesiva = 0.0;
    double error_absolut = 0.0;
    int n = 1; // N de subintervals
    char *razon_finalizacion = "Límite de iteraciones (n) superado.";

    // Capçalera de l'arxiu de sortida
    fprintf(archivo_salida, "n\tTn\t|Tn - Tn-1|\t|I - Tn|\n");

    // Bucle principal
    while (n <= N_MAX) {
        // Calculem l'aproximació
        Tn_actual = trap_full(n, A, B);

        // Calculem la diferencia i l'error
        if (n == 1) {
            // A la primera iteració no podem calcular-la
            diferencia_sucesiva = 0.0;
        } else {
            diferencia_sucesiva = fabs(Tn_actual - Tn_anterior);
        }

        error_absolut = fabs(I_EXACTE - Tn_actual);
        
        // Escribim a l'arxiu de sortida
        fprintf(archivo_salida, "%d\t%.15e\t%.15e\t%.15e\n", 
                n, Tn_actual, diferencia_sucesiva, error_absolut);
        
        // Criteri de parada - Si la diferencia és més petita que 10⁻6 parem
        if (n > 1 && diferencia_sucesiva < TOLERANCIA) {
            razon_finalizacion = "Criterio de convergencia (|Tn - Tn-1| < 1e-6) alcanzado.";
            break; 
        }

        // Preparem per a la següent iteració
        Tn_anterior = Tn_actual;
        n *= 2; // Dupliquem el nombre de subintervals
    }

    // Si el bucle ha acabat perquè n ha superat N_MAX
    if (n > N_MAX) {
        // Li assignem el valor de N_MAX a n
        n = N_MAX;
        // Recalculem els valors finals si hem duplicat n en excès
        Tn_actual = trap_full(n, A, B);
        error_absolut = fabs(I_EXACTE - Tn_actual);
    }
    
    fclose(archivo_salida);

    printf("\n----------------------------------------------------\n");
    printf("Proceso de integración numérica (Trapecio Compuesto) finalizado.\n");
    printf("Razón de finalización: %s\n", razon_finalizacion);
    printf("----------------------------------------------------\n");
    printf("Integral: I = ∫_%.1f^%.1f e^x dx = e - 1\n", A, B);
    printf("Valor Exacto (I):\t\t%.15e\n", I_EXACTE);
    printf("Valor Aproximado (Tn):\t%.15e (con n=%d subintervalos)\n", Tn_actual, n);
    printf("Error Absoluto Final:\t%.15e\n", error_absolut);
    printf("----------------------------------------------------\n");

    return 0;
}