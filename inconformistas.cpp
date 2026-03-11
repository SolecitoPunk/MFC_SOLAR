#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>


using namespace std;

// --- PARÁMETROS DEL MODELO ---
const double p0 = 0.01;
const double b = 0.02;
const double r = 0.1;
const double rb = r * b; // 0.002

// --- FUNCIONES AUXILIARES ---

// Solución Exacta: p(t) = 1 - (1 - p0) * e^(-rbt)
double get_exact(double t) {
    return 1.0 - (1.0 - p0) * exp(-rb * t);
}

// Funciones de iteración (Un solo paso)
double step_euler(double p, double h) {
    return p + h * rb * (1.0 - p);
}

double step_taylor(double p, double h) {
    // p_{n+1} = p_n + (hrb - (hrb)^2/2)(1 - p_n)
    double factor = (h * rb) - (pow(h * rb, 2) / 2.0);
    return p + factor * (1.0 - p);
}

double step_trapecio(double p, double h) {
    // p_{n+1} = [p_n(1 - hrb/2) + hrb] / (1 + hrb/2)
    double term = (h * rb) / 2.0;
    return (p * (1.0 - term) + h * rb) / (1.0 + term);
}

int main() {
    // ==========================================
    // PARTE 1: ANÁLISIS NUMÉRICO (Tabla y Errores)
    // ==========================================
    double h_sim = 1.0;  // Paso de 1 año
    double t_max = 50.0;
    int steps = (int)(t_max / h_sim);

    // Inicialización
    double t = 0;
    double p_e = p0, p_ta = p0, p_tr = p0;
    double p_ex = p0;

    cout << "=================================================================================" << endl;
    cout << "                       TABLA DE VALORES (h = " << h_sim << " anio)" << endl;
    cout << "=================================================================================" << endl;
    cout << left << setw(6) << "t" 
         << setw(15) << "Exacta" 
         << setw(15) << "Euler" 
         << setw(15) << "Taylor 2" 
         << setw(15) << "Trapecio" << endl;
    cout << string(66, '-') << endl;

    // Imprimir estado inicial
    cout << fixed << setprecision(6);
    cout << left << setw(6) << t 
         << setw(15) << p_ex 
         << setw(15) << p_e 
         << setw(15) << p_ta 
         << setw(15) << p_tr << endl;

    // Bucle de simulación paso a paso
    for(int i = 0; i < steps; i++) {
        // Calcular pasos
        p_e  = step_euler(p_e, h_sim);
        p_ta = step_taylor(p_ta, h_sim);
        p_tr = step_trapecio(p_tr, h_sim);
        
        t += h_sim;
        p_ex = get_exact(t);

        // Imprimir fila (mostramos cada 5 años para no saturar, o quita el 'if' para ver todo)
        if ((int)t % 5 == 0 || t == t_max) {
             cout << left << setw(6) << (int)t 
                  << setw(15) << p_ex 
                  << setw(15) << p_e 
                  << setw(15) << p_ta 
                  << setw(15) << p_tr << endl;
        }
    }

    // CÁLCULO DE ERRORES EN t = 50
    double err_e  = abs(p_ex - p_e);
    double err_ta = abs(p_ex - p_ta);
    double err_tr = abs(p_ex - p_tr);

    cout << endl;
    cout << "--- ANALISIS DE ERRORES EN t = 50 ---" << endl;
    cout << "Valor Exacto p(50): " << p_ex << endl;
    cout << "Error Euler:        " << scientific << err_e << fixed << endl;
    cout << "Error Taylor 2:     " << scientific << err_ta << fixed << endl;
    cout << "Error Trapecio:     " << scientific << err_tr << fixed << endl;
    cout << "-------------------------------------" << endl << endl;


    // ==========================================
    // PARTE 2: COMPARACIÓN DE EFICIENCIA (Benchmark)
    // ==========================================
    long long N_bench = 10000000; // 10^7 iteraciones
    double h_bench = t_max / N_bench;
    
    cout << "Iniciando Benchmark de Eficiencia (" << N_bench << " iteraciones)..." << endl;
    
    // Euler Benchmark
    auto start = chrono::high_resolution_clock::now();
    double p_temp = p0;
    double factor_e = h_bench * rb;
    for(long long i=0; i<N_bench; i++) {
        p_temp = p_temp + factor_e * (1.0 - p_temp);
    }
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> t_euler = end - start;
    if (p_temp > 1000.0) cout << "Dummy print: " << p_temp << endl;
    
    // Taylor Benchmark
    start = chrono::high_resolution_clock::now();
    p_temp = p0;
    double factor_ta = (h_bench * rb) - (pow(h_bench * rb, 2)/2.0);
    for(long long i=0; i<N_bench; i++) {
        p_temp = p_temp + factor_ta * (1.0 - p_temp);
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> t_taylor = end - start;
    if (p_temp > 1000.0) cout << "Dummy print: " << p_temp << endl;

    // Trapecio Benchmark
    start = chrono::high_resolution_clock::now();
    p_temp = p0;
    double term_bench = (h_bench * rb) / 2.0;
    double num_add = h_bench * rb;
    double den_bench = 1.0 + term_bench;
    double num_mult = 1.0 - term_bench;
    
    for(long long i=0; i<N_bench; i++) {
         p_temp = (p_temp * num_mult + num_add) / den_bench;
    }
    end = chrono::high_resolution_clock::now();
    chrono::duration<double> t_trap = end - start;
    if (p_temp > 1000.0) cout << "Dummy print: " << p_temp << endl;

    cout << "Tiempo Euler:    " << t_euler.count() << " s" << endl;
    cout << "Tiempo Taylor:   " << t_taylor.count() << " s" << endl;
    cout << "Tiempo Trapecio: " << t_trap.count() << " s" << endl;

    return 0;
}
