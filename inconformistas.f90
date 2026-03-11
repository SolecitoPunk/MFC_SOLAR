program rashevsky_full
    implicit none

    ! Parámetros físicos
    real(8), parameter :: p0 = 0.01d0
    real(8), parameter :: r = 0.1d0
    real(8), parameter :: b = 0.02d0
    real(8), parameter :: rb = r * b
    real(8), parameter :: t_max = 50.0d0
    
    ! Variables para simulación numérica
    real(8) :: h_sim = 1.0d0
    integer :: steps, i
    real(8) :: t
    real(8) :: p_e, p_ta, p_tr, p_ex
    real(8) :: err_e, err_ta, err_tr
    
    ! Variables para Taylor y Trapecio (Pre-cálculo)
    real(8) :: term_taylor, term_trap

    ! Variables para Benchmark
    integer(8), parameter :: N_bench = 10000000
    real(8) :: h_bench
    real(8) :: t_start, t_end
    real(8) :: time_e, time_ta, time_tr
    real(8) :: p_temp
    integer(8) :: j
    real(8) :: factor_bench_e, factor_bench_ta
    real(8) :: tr_num_mult, tr_num_add, tr_den

    ! ==========================================
    ! PARTE 1: ANÁLISIS NUMÉRICO (h = 1.0)
    ! ==========================================
    steps = int(t_max / h_sim)
    
    ! Inicializar
    t = 0.0d0
    p_e = p0
    p_ta = p0
    p_tr = p0
    p_ex = p0
    
    print *, "=========================================================================="
    print *, "                       TABLA DE RESULTADOS (h=1.0)"
    print *, "=========================================================================="
    print '(A6, A15, A15, A15, A15)', "t", "Exacta", "Euler", "Taylor 2", "Trapecio"
    print '(A66)', "------------------------------------------------------------------"
    
    ! Imprimir t=0
    print '(I6, 4F15.6)', int(t), p_ex, p_e, p_ta, p_tr
    
    ! Bucle temporal
    do i = 1, steps
        ! 1. Euler
        p_e = p_e + h_sim * rb * (1.0d0 - p_e)
        
        ! 2. Taylor Orden 2
        term_taylor = (h_sim * rb) - ((h_sim * rb)**2 / 2.0d0)
        p_ta = p_ta + term_taylor * (1.0d0 - p_ta)
        
        ! 3. Trapecio
        term_trap = (h_sim * rb) / 2.0d0
        p_tr = (p_tr * (1.0d0 - term_trap) + h_sim * rb) / (1.0d0 + term_trap)
        
        ! Actualizar tiempo y exacta
        t = t + h_sim
        p_ex = 1.0d0 - (1.0d0 - p0) * exp(-rb * t)
        
        ! Imprimir cada 5 pasos o al final
        if (mod(i, 5) == 0) then
            print '(I6, 4F15.6)', int(t), p_ex, p_e, p_ta, p_tr
        end if
    end do
    
    ! Cálculo de Error
    err_e = abs(p_ex - p_e)
    err_ta = abs(p_ex - p_ta)
    err_tr = abs(p_ex - p_tr)
    
    print *
    print *, "--- ERRORES ABSOLUTOS EN t=50 ---"
    print '(A20, F15.8)', "Valor Exacto:", p_ex
    print '(A20, E15.4)', "Error Euler:", err_e
    print '(A20, E15.4)', "Error Taylor:", err_ta
    print '(A20, E15.4)', "Error Trapecio:", err_tr
    print *
    
    ! ==========================================
    ! PARTE 2: BENCHMARK (10^7 iteraciones)
    ! ==========================================
    print *, "Iniciando Benchmark de Eficiencia..."
    h_bench = t_max / dble(N_bench)
    
    ! Euler
    call cpu_time(t_start)
    p_temp = p0
    factor_bench_e = h_bench * rb
    do j = 1, N_bench
        p_temp = p_temp + factor_bench_e * (1.0d0 - p_temp)
    end do
    call cpu_time(t_end)
    time_e = t_end - t_start
    print *, "Resultado Euler (Check): ", p_temp
    
    ! Taylor
    call cpu_time(t_start)
    p_temp = p0
    factor_bench_ta = (h_bench * rb) - ((h_bench * rb)**2 / 2.0d0)
    do j = 1, N_bench
        p_temp = p_temp + factor_bench_ta * (1.0d0 - p_temp)
    end do
    call cpu_time(t_end)
    time_ta = t_end - t_start
    print *, "Resultado Taylor (Check): ", p_temp

    ! Trapecio
    call cpu_time(t_start)
    p_temp = p0
    term_trap = (h_bench * rb) / 2.0d0
    tr_num_mult = 1.0d0 - term_trap
    tr_num_add = h_bench * rb
    tr_den = 1.0d0 + term_trap
    do j = 1, N_bench
        p_temp = (p_temp * tr_num_mult + tr_num_add) / tr_den
    end do
    call cpu_time(t_end)
    time_tr = t_end - t_start
    print *, "Resultado Trap (Check): ", p_temp
    
    print '(A15, F10.4, A)', "Tiempo Euler:", time_e, " s"
    print '(A15, F10.4, A)', "Tiempo Taylor:", time_ta, " s"
    print '(A15, F10.4, A)', "Tiempo Trap:", time_tr, " s"

end program rashevsky_full

