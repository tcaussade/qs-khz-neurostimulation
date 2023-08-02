TITLE Motor Axon Juxtaparanode Fast Potassium Channel

: Jan 2015
: Leo Medina
:
: This model is described in detail in:
:
: McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
: mammalian nerve fibers: influence of afterpotentials on the recovery
: cycle. Journal of Neurophysiology 87:995-1006, 2002.

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX juxtaparanode_channels
	NONSPECIFIC_CURRENT ikf
	RANGE gkfbar, ek
	RANGE n_inf
	RANGE tau_n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
    	
    gkfbar  = 0.002   (mho/cm2)	
	ek      = -90.0 (mV)	
	vtraub=-80
	anA = 0.0462
    anB = 83.2
    anC = 1.1
    bnA = 0.0824
    bnB = 66.0
    bnC = 10.5
}

STATE {
	n
}

ASSIGNED {
	v (mV)
    celsius (degC)
    ikf     (mA/cm2)    
	n_inf
	tau_n
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ikf = gkfbar * n*n*n*n * (v - ek)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
    evaluate_fct(v)
	n' = (n_inf - n) / tau_n
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:

	q10_3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	a = q10_3*vtrap9(v) 
	b = q10_3*vtrap10(v)
	tau_n = 1 / (a + b)
	n_inf = a / (a + b)
}

FUNCTION vtrap9(x) {
	if (fabs((x+anB)/anC) < 1e-6) {
		vtrap9 = anA*anC 
	}else{
		vtrap9 = (anA*(x+anB)) / (1 - Exp(-(x+anB)/anC)) 
	}
}

FUNCTION vtrap10(x) {
	if (fabs((x+bnB)/bnC) < 1e-6) {
		vtrap10 = bnA*bnC 
	}else{
		vtrap10 = (bnA*(-(x+bnB))) / (1 - Exp((x+bnB)/bnC))
	}
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON