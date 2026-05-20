program correccion_scalar
	use polylog_module
	implicit none
	
	! Constantes físicas
	double precision :: pi, alpha, e
	real(8) :: B
	double precision :: m, m2, w, mf, w0

	! Resultados intermedios
	double precision :: lambda_val,old_lambda, delta_val, old_delta, mu_val, old_mu
	double precision :: B0M, B0V, DB0M, C0M, A0M
	double precision :: dM, MV, MA, MT, dG, MB, f, dGC, MS

	! --- Funciones externas de LoopTools ---
	double precision B0, B1, DB0, DB1, C0, A0
	external B0, B1, DB0, DB1, C0, A0
	external ltini

	! --- Subrutina SetLambda ---
	double precision SetLambda
	external SetLambda

	! --- Subrutina SetDelta ---
	double precision SetDelta
	external SetDelta

	! --- Subrutina SetMudim
	double precision SetMudim
	external SetMudim

	! --- Inicialización de constantes ---
	pi = 3.141592653589793d0
	alpha = 1d0/137.035999d0
	e = dsqrt(4d0*pi*alpha)

	! --- Entradas de masa (MeV) ---
	print *, "Ingrese la masa de la partícula escalar dispersada (MeV):"
	read *, m
	print *, "Ingrese la masa de la partícula que decae (MeV): "
	read *, mf
	m2 = m**2
	w = (1d0/2d0)*(mf**2 - 2d0*m2)
	B = dsqrt(1-4d0*m2/(mf**2))

	! --- Lectura del valor de parámetros ---
	print *, "Introduzca el valor de Lambda:"
	read *, lambda_val
	print *, "Introduzca el valor de Delta"
	read *, delta_val
	print *, "Introduzca el valor de Mu"
	read *, mu_val
	print *, "Introduzca el umbral de detección de fotones"
	read *, w0

	! --- Inicializar LoopTools ---
	call ltini()

	! --- Cambio del valor de parámetros ---
	old_lambda = SetLambda(lambda_val)
	old_delta = SetDelta(delta_val)
	old_mu = SetMudim(mu_val)

	! --- Funciones de Passarino Veltman evaluadas en m^2 y mf^2 ---
	
	A0M = A0(m2)
	B0M = B0(m2, 0d0, m2)
	B0V = B0(mf**2, m2, m2)
	DB0M = DB0(m2, 0d0, m2)
	C0M = C0(m2, m2, mf**2, m2, 0d0, m2)

	! --- Ecuaciones de FeynCalc ---

	dM = (e**2/(16d0*pi**2))*(B0M + 2d0*m2*DB0M)
	
	MV = -(e**2/(16d0*pi**2))*(A0M/m2 - (4d0*w/(m2-w))*B0V + 4d0*w*C0M - (2d0*(m2-3d0*w)/(m2-w))*B0M)

	MA = (e**2/(16d0*pi**2d0))*(A0M/m2 - 4d0*B0M)

	MS = (e**2/(2d0*pi)**2)*((2d0*B**2/(1d0-B**2))*(log(4*w0**2/lambda_val)-(1/B)*log((1+B)/(1-B)))-&
		log(w0**2/lambda_val)*log((1+B)/(1-B))+2d0*(Li2(B)-Li2(-B)))

	f = MA + MV + 2d0*dM

	MT = 1d0 + 2d0*f + MS

	dG = ((4.6d0**2d0)/(48d0*pi))*mf*((1d0 - (4d0*m2)/mf**2)**(3d0/2d0))

	dGC = dG*MT

	print *, "--- Resultados ---"
	print *, "El ancho de decaimiento a tree level es:", dG 
	print *, "El ancho de decaimiento a NLO es:", dGC

end program correccion_scalar