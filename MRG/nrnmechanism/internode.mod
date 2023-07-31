TITLE Hyperpolarization-activated current

: Jan 2015
: Leo Medina 
:
: This current is described in detail in:
:
: Howells, Trevillion, Bostock, Burke (2012). The voltage dependence of Ih in human myelinated axons. J Physiol 590(7)

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX internode_channels
	NONSPECIFIC_CURRENT ih
	RANGE ghbar, eh
	RANGE q_inf
	RANGE tau_q
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ghbar  = 0.000025   (mho/cm2)
	eh	= -52.7 (mV)
	    
    aqA = 0.00522
    aqB = -94.2
    aqC = -12.2
    bqA = 0.00522
    bqB = -94.2
    bqC = -12.2
    
}

STATE {
	q
}

ASSIGNED {
    v               (mV)
    celsius		(degC)		
	ih      (mA/cm2)
	q_inf
	tau_q
	q10_3
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ih = ghbar * q * (v - eh)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations
    evaluate_fct(v)
	q' = (q_inf - q) / tau_q
}

UNITSOFF

INITIAL {
:
:	Q10 adjustment
:

	q10_3 = 3.0 ^ ((celsius-36)/ 10 )

	evaluate_fct(v)
	q = q_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b

	a = q10_3*aqA*Exp((v-aqB)/aqC)
	b = q10_3*bqA/Exp((v-bqB)/bqC)
	tau_q = 1 / (a + b)
	q_inf = a / (a + b)
	
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	}else{
		Exp = exp(x)
	}
}

UNITSON