: This channel is implemented by Nathan Titus 2019. 
: M-type potassium channel which directly represents the contributions of
: Kv7.2 (n) and Kv7.3 (m). Kv7.5 was considered "similar enough" to
: Kv7.2 so as to be considered a part of that current (n). 
: Data From:

NEURON {
	SUFFIX kv12
	USEION k READ ek WRITE ik
	RANGE gbar, ek, ik
	RANGE ninf,tau_n
	RANGE ninfshift, ik, gp, g
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gbar 	(S/cm2) 
	q10 = 3
        
	ninfshift = 0 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ik	(mA/cm2)
	g	(S/cm2)
    tau_n  (ms)
    ninf
	gp
    ek	(mV)
	celsius (degC)
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gp = n^3
	g = gbar*gp
	ik = g * (v-ek)
}

INITIAL {
	: assume that equilibrium has been reached
    rates(v)    
    n=ninf :slow

}

DERIVATIVE states {
	rates(v)
    n' = (ninf - n)/tau_n
          
}

? rates
PROCEDURE rates(Vm (mV)) (/ms) {    
	LOCAL Q10
		TABLE ninf,tau_n DEPEND celsius FROM -120 TO 100 WITH 440
UNITSOFF
		Q10 = q10^((celsius-22)/10)
        ninf= (1/(1+exp((-1*(Vm + 21)/7))))^(1/3)
        tau_n = 0.4 + 103/(exp((Vm+76)/28)+exp(-1*(Vm+13)/23))
		
        tau_n=tau_n/Q10/2
UNITSON

}
