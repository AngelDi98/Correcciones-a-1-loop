program convection_term
	use polylog_module
	implicit none
	
	! Constantes físicas
	double precision :: pi, alpha, e
	real(8) :: B
	double precision :: m, m2, w, mf, w0

	! Resultados intermedios
	double precision :: lambda_val,old_lambda, delta_val, old_delta, mu_val, old_mu
	double precision :: B0M, B0V, C0M, C0B
	double precision :: MV, MT, dG, dGC, MS

	! --- Funciones externas de LoopTools ---
	double precision B0, C0
	external B0, C0
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
	
	B0M = B0(m2, 0d0, m2)
	BOV = B0(mf**2, m2, m2)
	C0M = C0(m2, m2, mf**2, m2, 0d0, m2)
	C0B = C0(0d0, m2, m2, m2, m2, 0d0)

	! --- Ecuaciones de FeynCalc ---
	
	MV = -(1/(4d0*pi))*(B0V - B0M + 4d0*w*C0M + 4d0*m2*C0B)

	MS = -(1/(4d0*pi))*((2d0*B**2/(1d0-B**2))*(log(4*w0**2/lambda_val)-(1/B)*log((1+B)/(1-B)))-&
		log(w0**2/lambda_val)*log((1+B)/(1-B))+2d0*(Li2(B)-Li2(-B)))

	MT = 1d0 + (2d0*alpha)*(MV + MS)

	dG = ((4.6d0**2d0)/(48d0*pi))*(mf**2d0)*((1d0 - (4d0*m2)/mf**2)**(3d0/2d0))

	dGC = dG*MT

	print *, "--- Resultados ---"
	print *, "El ancho de decaimiento es a tree level es:", dG
	print *, "El ancho de decaimiento al NLO es:", dGC 

end program convection_term