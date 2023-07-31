TITLE DCFiber node channels

: 2020
: Leo Medina based on work by Cameron C. McIntyre (McIntyre et al, 2002)
:
: Slow K+, and Leakage currents
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX node_ikil
	NONSPECIFIC_CURRENT ik
	NONSPECIFIC_CURRENT il
	RANGE gkbar, gl, ek, el
	RANGE s_inf
	RANGE tau_s
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	gkbar   = 0.04 	(mho/cm2)
	gl	= 0.007 (mho/cm2)
	ek      = -90.0 (mV)
	el	= -90.0 (mV)
	vtraub=-80
	asA = 0.3
	asB = -27
	asC = -5
	bsA = 0.03
	bsB = 10
	bsC = -1

}

STATE {
	s
}

ASSIGNED {
    v (mV)
    celsius		(degC)
	ik      (mA/cm2)
	il      (mA/cm2)
	s_inf
	tau_s
	q10_1
	q10_2
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik   = gkbar * s * (v - ek)
	il   = gl * (v - el)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
    evaluate_fct(v)
	s' = (s_inf - s) / tau_s
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:
	q10_3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	s = s_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = q10_3*asA / (Exp((v2+asB)/asC) + 1)
	b = q10_3*bsA / (Exp((v2+bsB)/bsC) + 1)
	tau_s = 1 / (a + b)
	s_inf = a / (a + b)

}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON