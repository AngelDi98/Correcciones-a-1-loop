program compare_methods
    use polylog_module
    implicit none
    
    ! Constantes físicas
    double precision :: pi, alpha, e
    double precision :: m, m2, w, mf, w0
    real(8) :: B

    ! Resultados intermedios
    double precision :: lambda_val, old_lambda
    double precision :: delta_val, old_delta
    double precision :: mu_val, old_mu

    double precision :: B0M, B0V, DB0M, C0M, C0B, A0M
    double precision :: dM, MVa, MA, f, MTa
    double precision :: MVb, MTb, MS
    double precision :: dGa, dGCa, dGCb

    integer :: N
    integer :: unit_total, unit_conv

    ! --- Variables para el calculo de error ---

    double precision :: R, eps
    double precision :: sum_R2, sum_eps2
    double precision :: L2_norm, RMS_rel, max_rel

    ! --- Funciones externas de LoopTools ---
    double precision B0, B1, DB0, DB1, C0, A0
    external B0, B1, DB0, DB1, C0, A0
    external ltini

    double precision SetLambda
    external SetLambda

    double precision SetDelta
    external SetDelta

    double precision SetMudim
    external SetMudim

    ! --- Inicialización de constantes ---
    pi = 3.141592653589793d0
    alpha = 1d0/137.035999d0
    e = dsqrt(4d0*pi*alpha)

    ! --- Entradas de masa (MeV) ---
    m = 494d0
    print *, "La masa de los kaones se mantiene fija K^{+/-}=", m, "MeV"
    mf = 1019d0
    print *, "La masa del phi variará de 1 en 1 MeV, iniciando en:", mf, "MeV"
    m2 = m**2
    w0 = 5d0
    print *, "El umbral de detección de fotones suaves fijará en:", w0, "MeV"

    ! --- Lectura de parámetros ---
    print *, "Introduzca el valor de Lambda:"
    read *, lambda_val
    print *, "Introduzca el valor de Delta"
    read *, delta_val
    print *, "Introduzca el valor de Mu"
    read *, mu_val

    call ltini()

    old_lambda = SetLambda(lambda_val)
    old_delta  = SetDelta(delta_val)
    old_mu     = SetMudim(mu_val)

    ! --- Abrir archivos CSV ---
    unit_total = 10
    unit_conv  = 20

    open(unit=unit_total, file="total_method.csv", status="replace", action="write")
    open(unit=unit_conv,  file="convection_method.csv", status="replace", action="write")

    write(unit_total,'(A)') "mf,dGC_total"
    write(unit_conv ,'(A)') "mf,dGC_convection"

    sum_R2 = 0d0
    sum_eps2 = 0d0
    max_rel = 0d0
    N = 0

    do while(mf <= 3000d0)

        w = (1d0/2d0)*(mf**2 - 2d0*m2)
	B = dsqrt(1-4d0*m2/(mf**2))

        ! --- Funciones Passarino-Veltman ---
        A0M  = A0(m2)
        B0M  = B0(m2, 0d0, m2)
        B0V  = B0(mf**2, m2, m2)
        DB0M = DB0(m2, 0d0, m2)
        C0M  = C0(m2, m2, mf**2, m2, 0d0, m2)
        C0B  = C0(0d0, m2, m2, m2, m2, 0d0)

	! --- Factor de corrección soft-bremsstrahlung ---

	MS = ((2d0*B**2/(1d0-B**2))*(log(4*w0**2/lambda_val)-(1/B)*log((1+B)/(1-B)))-&
		log(w0**2/lambda_val)*log((1+B)/(1-B))+2d0*(Li2(B)-Li2(-B)))

        ! ==============================
        ! MÉTODO TOTAL
        ! ==============================

        dM = (e**2/(16d0*pi**2))*(B0M + 2d0*m2*DB0M)

        MVa = -(e**2/(16d0*pi**2))*(A0M/m2 - (4d0*w/(m2-w))*B0V + 4d0*w*C0M - (2d0*(m2-3d0*w)/(m2-w))*B0M)

        MA = (e**2/(16d0*pi**2))*(A0M/m2 - 4d0*B0M)

        f = MA + MVa + 2d0*dM
        MTa = 1d0 + 2d0*f + (alpha/pi)*MS

        dGa  = ((4.6d0**2d0)/(48d0*pi))*mf*((1d0 - (4d0*m2)/mf**2)**(3d0/2d0))
        dGCa = dGa * MTa

        ! ==============================
        ! CONVECTION TERM
        ! ==============================

        MVb = -(1/(4d0*pi))*(B0V - B0M + 4d0*w*C0M + 4d0*m2*C0B)
        MTb = 1d0 + 2d0*alpha*MVb - (alpha/(2d0*pi))*MS
        dGCb = dGa * MTb

        ! --- Guardar tuplas (mf, dGC) ---
        write(unit_total,'(ES18.10,",",ES18.10)') mf, dGCa
        write(unit_conv ,'(ES18.10,",",ES18.10)') mf, dGCb

	! --- Residuo absoluto ---
	R = dGCa - dGCb

	! --- Error relativo ---
	if (dGCa /= 0d0) then
    	    eps = R / dGCa
	else
    	    eps = 0d0
	end if

	! --- Acumuladores ---
	sum_R2 = sum_R2 + R**2
	sum_eps2 = sum_eps2 + eps**2

	if (abs(eps) > max_rel) then
    	    max_rel = abs(eps)
	end if

	N = N + 1
        mf = mf + 1d0

    end do

    close(unit_total)
    close(unit_conv)

    L2_norm = dsqrt(sum_R2 / dble(N))
    RMS_rel = dsqrt(sum_eps2 / dble(N))

    print *, "====================================="
    print *, "Resultados del análisis:"
    print *, "Norma L2 absoluta     =", L2_norm
    print *, "RMS error relativo    =", RMS_rel
    print *, "Error relativo máximo =", max_rel
    print *, "====================================="


end program compare_methods