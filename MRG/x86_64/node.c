/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__node_channels
#define _nrn_initial _nrn_initial__node_channels
#define nrn_cur _nrn_cur__node_channels
#define _nrn_current _nrn_current__node_channels
#define nrn_jacob _nrn_jacob__node_channels
#define nrn_state _nrn_state__node_channels
#define _net_receive _net_receive__node_channels 
#define evaluate_fct evaluate_fct__node_channels 
#define states states__node_channels 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnapbar _p[0]
#define gnapbar_columnindex 0
#define gnabar _p[1]
#define gnabar_columnindex 1
#define gkbar _p[2]
#define gkbar_columnindex 2
#define gl _p[3]
#define gl_columnindex 3
#define ena _p[4]
#define ena_columnindex 4
#define ek _p[5]
#define ek_columnindex 5
#define el _p[6]
#define el_columnindex 6
#define inap _p[7]
#define inap_columnindex 7
#define ina _p[8]
#define ina_columnindex 8
#define ik _p[9]
#define ik_columnindex 9
#define il _p[10]
#define il_columnindex 10
#define mp_inf _p[11]
#define mp_inf_columnindex 11
#define m_inf _p[12]
#define m_inf_columnindex 12
#define h_inf _p[13]
#define h_inf_columnindex 13
#define s_inf _p[14]
#define s_inf_columnindex 14
#define hs_inf _p[15]
#define hs_inf_columnindex 15
#define tau_mp _p[16]
#define tau_mp_columnindex 16
#define tau_m _p[17]
#define tau_m_columnindex 17
#define tau_h _p[18]
#define tau_h_columnindex 18
#define tau_s _p[19]
#define tau_s_columnindex 19
#define tau_hs _p[20]
#define tau_hs_columnindex 20
#define mp _p[21]
#define mp_columnindex 21
#define m _p[22]
#define m_columnindex 22
#define h _p[23]
#define h_columnindex 23
#define s _p[24]
#define s_columnindex 24
#define hs _p[25]
#define hs_columnindex 25
#define Dmp _p[26]
#define Dmp_columnindex 26
#define Dm _p[27]
#define Dm_columnindex 27
#define Dh _p[28]
#define Dh_columnindex 28
#define Ds _p[29]
#define Ds_columnindex 29
#define Dhs _p[30]
#define Dhs_columnindex 30
#define q10_1 _p[31]
#define q10_1_columnindex 31
#define q10_2 _p[32]
#define q10_2_columnindex 32
#define q10_3 _p[33]
#define q10_3_columnindex 33
#define v _p[34]
#define v_columnindex 34
#define _g _p[35]
#define _g_columnindex 35
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_Exp(void);
 static void _hoc_evaluate_fct(void);
 static void _hoc_vtrap(void);
 static void _hoc_vtrap8(void);
 static void _hoc_vtrap7(void);
 static void _hoc_vtrap6(void);
 static void _hoc_vtrap2(void);
 static void _hoc_vtrap1(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_node_channels", _hoc_setdata,
 "Exp_node_channels", _hoc_Exp,
 "evaluate_fct_node_channels", _hoc_evaluate_fct,
 "vtrap_node_channels", _hoc_vtrap,
 "vtrap8_node_channels", _hoc_vtrap8,
 "vtrap7_node_channels", _hoc_vtrap7,
 "vtrap6_node_channels", _hoc_vtrap6,
 "vtrap2_node_channels", _hoc_vtrap2,
 "vtrap1_node_channels", _hoc_vtrap1,
 0, 0
};
#define Exp Exp_node_channels
#define vtrap vtrap_node_channels
#define vtrap8 vtrap8_node_channels
#define vtrap7 vtrap7_node_channels
#define vtrap6 vtrap6_node_channels
#define vtrap2 vtrap2_node_channels
#define vtrap1 vtrap1_node_channels
 extern double Exp( _threadargsprotocomma_ double );
 extern double vtrap( _threadargsprotocomma_ double );
 extern double vtrap8( _threadargsprotocomma_ double );
 extern double vtrap7( _threadargsprotocomma_ double );
 extern double vtrap6( _threadargsprotocomma_ double );
 extern double vtrap2( _threadargsprotocomma_ double );
 extern double vtrap1( _threadargsprotocomma_ double );
 /* declare global and static user variables */
#define HSbiasV HSbiasV_node_channels
 double HSbiasV = -15;
#define SIF SIF_node_channels
 double SIF = 10000;
#define asC asC_node_channels
 double asC = -5;
#define asB asB_node_channels
 double asB = -27;
#define asA asA_node_channels
 double asA = 0.3;
#define ahC ahC_node_channels
 double ahC = 11;
#define ahB ahB_node_channels
 double ahB = 114;
#define ahA ahA_node_channels
 double ahA = 0.062;
#define amC amC_node_channels
 double amC = 10.3;
#define amB amB_node_channels
 double amB = 21.4;
#define amA amA_node_channels
 double amA = 1.86;
#define ampC ampC_node_channels
 double ampC = 10.2;
#define ampB ampB_node_channels
 double ampB = 27;
#define ampA ampA_node_channels
 double ampA = 0.01;
#define bsC bsC_node_channels
 double bsC = -1;
#define bsB bsB_node_channels
 double bsB = 10;
#define bsA bsA_node_channels
 double bsA = 0.03;
#define bhC bhC_node_channels
 double bhC = 13.4;
#define bhB bhB_node_channels
 double bhB = 31.8;
#define bhA bhA_node_channels
 double bhA = 2.3;
#define bmC bmC_node_channels
 double bmC = 9.16;
#define bmB bmB_node_channels
 double bmB = 25.7;
#define bmA bmA_node_channels
 double bmA = 0.086;
#define bmpC bmpC_node_channels
 double bmpC = 10;
#define bmpB bmpB_node_channels
 double bmpB = 34;
#define bmpA bmpA_node_channels
 double bmpA = 0.00025;
#define vtraub vtraub_node_channels
 double vtraub = -80;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gnapbar_node_channels", "mho/cm2",
 "gnabar_node_channels", "mho/cm2",
 "gkbar_node_channels", "mho/cm2",
 "gl_node_channels", "mho/cm2",
 "ena_node_channels", "mV",
 "ek_node_channels", "mV",
 "el_node_channels", "mV",
 "inap_node_channels", "mA/cm2",
 "ina_node_channels", "mA/cm2",
 "ik_node_channels", "mA/cm2",
 "il_node_channels", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double hs0 = 0;
 static double h0 = 0;
 static double m0 = 0;
 static double mp0 = 0;
 static double s0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "vtraub_node_channels", &vtraub_node_channels,
 "ampA_node_channels", &ampA_node_channels,
 "ampB_node_channels", &ampB_node_channels,
 "ampC_node_channels", &ampC_node_channels,
 "bmpA_node_channels", &bmpA_node_channels,
 "bmpB_node_channels", &bmpB_node_channels,
 "bmpC_node_channels", &bmpC_node_channels,
 "amA_node_channels", &amA_node_channels,
 "amB_node_channels", &amB_node_channels,
 "amC_node_channels", &amC_node_channels,
 "bmA_node_channels", &bmA_node_channels,
 "bmB_node_channels", &bmB_node_channels,
 "bmC_node_channels", &bmC_node_channels,
 "ahA_node_channels", &ahA_node_channels,
 "ahB_node_channels", &ahB_node_channels,
 "ahC_node_channels", &ahC_node_channels,
 "bhA_node_channels", &bhA_node_channels,
 "bhB_node_channels", &bhB_node_channels,
 "bhC_node_channels", &bhC_node_channels,
 "asA_node_channels", &asA_node_channels,
 "asB_node_channels", &asB_node_channels,
 "asC_node_channels", &asC_node_channels,
 "bsA_node_channels", &bsA_node_channels,
 "bsB_node_channels", &bsB_node_channels,
 "bsC_node_channels", &bsC_node_channels,
 "HSbiasV_node_channels", &HSbiasV_node_channels,
 "SIF_node_channels", &SIF_node_channels,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"node_channels",
 "gnapbar_node_channels",
 "gnabar_node_channels",
 "gkbar_node_channels",
 "gl_node_channels",
 "ena_node_channels",
 "ek_node_channels",
 "el_node_channels",
 0,
 "inap_node_channels",
 "ina_node_channels",
 "ik_node_channels",
 "il_node_channels",
 "mp_inf_node_channels",
 "m_inf_node_channels",
 "h_inf_node_channels",
 "s_inf_node_channels",
 "hs_inf_node_channels",
 "tau_mp_node_channels",
 "tau_m_node_channels",
 "tau_h_node_channels",
 "tau_s_node_channels",
 "tau_hs_node_channels",
 0,
 "mp_node_channels",
 "m_node_channels",
 "h_node_channels",
 "s_node_channels",
 "hs_node_channels",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 36, _prop);
 	/*initialize range parameters*/
 	gnapbar = 0.005;
 	gnabar = 3;
 	gkbar = 0.04;
 	gl = 0.007;
 	ena = 50;
 	ek = -90;
 	el = -90;
 	_prop->param = _p;
 	_prop->param_size = 36;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _node_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 36, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 node_channels /hpc/home/ts490/dcfiber/nrnmechanism/node.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "DCFiber node channels";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int evaluate_fct(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[5], _dlist1[5];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   evaluate_fct ( _threadargscomma_ v ) ;
   Dmp = ( mp_inf - mp ) / tau_mp ;
   Dm = ( m_inf - m ) / tau_m ;
   Dh = ( h_inf - h ) / tau_h ;
   Ds = ( s_inf - s ) / tau_s ;
   Dhs = ( hs_inf - hs ) / tau_hs ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 evaluate_fct ( _threadargscomma_ v ) ;
 Dmp = Dmp  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_mp )) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_h )) ;
 Ds = Ds  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_s )) ;
 Dhs = Dhs  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_hs )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   evaluate_fct ( _threadargscomma_ v ) ;
    mp = mp + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_mp)))*(- ( ( ( mp_inf ) ) / tau_mp ) / ( ( ( ( - 1.0 ) ) ) / tau_mp ) - mp) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( m_inf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_h)))*(- ( ( ( h_inf ) ) / tau_h ) / ( ( ( ( - 1.0 ) ) ) / tau_h ) - h) ;
    s = s + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_s)))*(- ( ( ( s_inf ) ) / tau_s ) / ( ( ( ( - 1.0 ) ) ) / tau_s ) - s) ;
    hs = hs + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_hs)))*(- ( ( ( hs_inf ) ) / tau_hs ) / ( ( ( ( - 1.0 ) ) ) / tau_hs ) - hs) ;
   }
  return 0;
}
 
static int  evaluate_fct ( _threadargsprotocomma_ double _lv ) {
   double _la , _lb , _lv2 ;
 _la = q10_1 * vtrap1 ( _threadargscomma_ _lv ) ;
   _lb = q10_1 * vtrap2 ( _threadargscomma_ _lv ) ;
   tau_mp = 1.0 / ( _la + _lb ) ;
   mp_inf = _la / ( _la + _lb ) ;
   _la = q10_1 * vtrap6 ( _threadargscomma_ _lv ) ;
   _lb = q10_1 * vtrap7 ( _threadargscomma_ _lv ) ;
   tau_m = 1.0 / ( _la + _lb ) ;
   m_inf = _la / ( _la + _lb ) ;
   _la = q10_2 * vtrap8 ( _threadargscomma_ _lv ) ;
   _lb = q10_2 * bhA / ( 1.0 + Exp ( _threadargscomma_ - ( _lv + bhB ) / bhC ) ) ;
   tau_h = 1.0 / ( _la + _lb ) ;
   h_inf = _la / ( _la + _lb ) ;
   _la = q10_2 * vtrap8 ( _threadargscomma_ _lv + HSbiasV ) ;
   _lb = q10_2 * bhA / ( 1.0 + Exp ( _threadargscomma_ - ( _lv + HSbiasV + bhB ) / bhC ) ) ;
   tau_hs = SIF / ( _la + _lb ) ;
   hs_inf = _la / ( _la + _lb ) ;
   _lv2 = _lv - vtraub ;
   _la = q10_3 * asA / ( Exp ( _threadargscomma_ ( _lv2 + asB ) / asC ) + 1.0 ) ;
   _lb = q10_3 * bsA / ( Exp ( _threadargscomma_ ( _lv2 + bsB ) / bsC ) + 1.0 ) ;
   tau_s = 1.0 / ( _la + _lb ) ;
   s_inf = _la / ( _la + _lb ) ;
    return 0; }
 
static void _hoc_evaluate_fct(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 evaluate_fct ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap;
 if ( _lx < - 50.0 ) {
     _lvtrap = 0.0 ;
     }
   else {
     _lvtrap = bsA / ( Exp ( _threadargscomma_ ( _lx + bsB ) / bsC ) + 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap1 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap1;
 if ( fabs ( ( _lx + ampB ) / ampC ) < 1e-6 ) {
     _lvtrap1 = ampA * ampC ;
     }
   else {
     _lvtrap1 = ( ampA * ( _lx + ampB ) ) / ( 1.0 - Exp ( _threadargscomma_ - ( _lx + ampB ) / ampC ) ) ;
     }
   
return _lvtrap1;
 }
 
static void _hoc_vtrap1(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap1 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap2 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap2;
 if ( fabs ( ( _lx + bmpB ) / bmpC ) < 1e-6 ) {
     _lvtrap2 = bmpA * bmpC ;
     }
   else {
     _lvtrap2 = ( bmpA * ( - ( _lx + bmpB ) ) ) / ( 1.0 - Exp ( _threadargscomma_ ( _lx + bmpB ) / bmpC ) ) ;
     }
   
return _lvtrap2;
 }
 
static void _hoc_vtrap2(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap2 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap6 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap6;
 if ( fabs ( ( _lx + amB ) / amC ) < 1e-6 ) {
     _lvtrap6 = amA * amC ;
     }
   else {
     _lvtrap6 = ( amA * ( _lx + amB ) ) / ( 1.0 - Exp ( _threadargscomma_ - ( _lx + amB ) / amC ) ) ;
     }
   
return _lvtrap6;
 }
 
static void _hoc_vtrap6(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap6 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap7 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap7;
 if ( fabs ( ( _lx + bmB ) / bmC ) < 1e-6 ) {
     _lvtrap7 = bmA * bmC ;
     }
   else {
     _lvtrap7 = ( bmA * ( - ( _lx + bmB ) ) ) / ( 1.0 - Exp ( _threadargscomma_ ( _lx + bmB ) / bmC ) ) ;
     }
   
return _lvtrap7;
 }
 
static void _hoc_vtrap7(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap7 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap8 ( _threadargsprotocomma_ double _lx ) {
   double _lvtrap8;
 if ( fabs ( ( _lx + ahB ) / ahC ) < 1e-6 ) {
     _lvtrap8 = ahA * ahC ;
     }
   else {
     _lvtrap8 = ( ahA * ( - ( _lx + ahB ) ) ) / ( 1.0 - Exp ( _threadargscomma_ ( _lx + ahB ) / ahC ) ) ;
     }
   
return _lvtrap8;
 }
 
static void _hoc_vtrap8(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap8 ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double Exp ( _threadargsprotocomma_ double _lx ) {
   double _lExp;
 if ( _lx < - 100.0 ) {
     _lExp = 0.0 ;
     }
   else {
     _lExp = exp ( _lx ) ;
     }
   
return _lExp;
 }
 
static void _hoc_Exp(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  Exp ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  hs = hs0;
  h = h0;
  m = m0;
  mp = mp0;
  s = s0;
 {
   q10_1 = pow( 2.2 , ( ( celsius - 20.0 ) / 10.0 ) ) ;
   q10_2 = pow( 2.9 , ( ( celsius - 20.0 ) / 10.0 ) ) ;
   q10_3 = pow( 3.0 , ( ( celsius - 36.0 ) / 10.0 ) ) ;
   evaluate_fct ( _threadargscomma_ v ) ;
   mp = mp_inf ;
   m = m_inf ;
   h = h_inf ;
   s = s_inf ;
   hs = hs_inf ;
   }
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   inap = gnapbar * mp * mp * mp * ( v - ena ) ;
   ina = gnabar * m * m * m * h * hs * ( v - ena ) ;
   ik = gkbar * s * ( v - ek ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += inap;
 _current += ik;
 _current += il;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
 {   states(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = mp_columnindex;  _dlist1[0] = Dmp_columnindex;
 _slist1[1] = m_columnindex;  _dlist1[1] = Dm_columnindex;
 _slist1[2] = h_columnindex;  _dlist1[2] = Dh_columnindex;
 _slist1[3] = s_columnindex;  _dlist1[3] = Ds_columnindex;
 _slist1[4] = hs_columnindex;  _dlist1[4] = Dhs_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/hpc/home/ts490/dcfiber/nrnmechanism/node.mod";
static const char* nmodl_file_text = 
  "TITLE DCFiber node channels\n"
  "\n"
  ": 2016\n"
  ": Leo Medina based on work by Cameron C. McIntyre (McIntyre et al, 2002)\n"
  ":\n"
  ": Fast Na+ with slow inactivation, Persistent Na+, Slow K+, and Leakage currents\n"
  ": responsible for nodal action potential\n"
  ": Iterative equations H-H notation rest = -80 mV\n"
  ":\n"
  ": This model is described in detail in:\n"
  ":\n"
  "\n"
  "INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX node_channels\n"
  "	NONSPECIFIC_CURRENT ina\n"
  "	NONSPECIFIC_CURRENT inap\n"
  "	NONSPECIFIC_CURRENT ik\n"
  "	NONSPECIFIC_CURRENT il\n"
  "	RANGE gnapbar, gnabar, gkbar, gl, ena, ek, el\n"
  "	RANGE mp_inf, m_inf, h_inf, s_inf, hs_inf\n"
  "	RANGE tau_mp, tau_m, tau_h, tau_s, tau_hs\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "\n"
  "	gnapbar = 0.005	(mho/cm2)\n"
  "	gnabar	= 3.0	(mho/cm2)\n"
  "	gkbar   = 0.04 	(mho/cm2)\n"
  "	gl	= 0.007 (mho/cm2)\n"
  "	ena     = 50.0  (mV)\n"
  "	ek      = -90.0 (mV)\n"
  "	el	= -90.0 (mV)\n"
  "	vtraub=-80\n"
  "	ampA = 0.01\n"
  "	ampB = 27\n"
  "	ampC = 10.2\n"
  "	bmpA = 0.00025\n"
  "	bmpB = 34\n"
  "	bmpC = 10\n"
  "	amA = 1.86\n"
  "	amB = 21.4\n"
  "	amC = 10.3\n"
  "	bmA = 0.086\n"
  "	bmB = 25.7\n"
  "	bmC = 9.16\n"
  "	ahA = 0.062\n"
  "	ahB = 114.0\n"
  "	ahC = 11.0\n"
  "	bhA = 2.3\n"
  "	bhB = 31.8\n"
  "	bhC = 13.4\n"
  "	asA = 0.3\n"
  "	asB = -27\n"
  "	asC = -5\n"
  "	bsA = 0.03\n"
  "	bsB = 10\n"
  "	bsC = -1\n"
  "	HSbiasV = -15 : bias voltage in hs (Elliot, 1997)\n"
  "	SIF = 10000  : speed of hs with respecto to h (Elliot, 1997)\n"
  "\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	mp m h s hs\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "    v (mV)\n"
  "    celsius		(degC)\n"
  "	inap    (mA/cm2)\n"
  "	ina	(mA/cm2)\n"
  "	ik      (mA/cm2)\n"
  "	il      (mA/cm2)\n"
  "	mp_inf\n"
  "	m_inf\n"
  "	h_inf\n"
  "	s_inf\n"
  "	hs_inf\n"
  "	tau_mp\n"
  "	tau_m\n"
  "	tau_h\n"
  "	tau_s\n"
  "	tau_hs\n"
  "	q10_1\n"
  "	q10_2\n"
  "	q10_3\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	inap = gnapbar * mp*mp*mp * (v - ena)\n"
  "	ina = gnabar * m * m * m * h * hs *(v - ena)\n"
  "	ik   = gkbar * s * (v - ek)\n"
  "	il   = gl * (v - el)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {   : exact Hodgkin-Huxley equations\n"
  "       evaluate_fct(v)\n"
  "	mp'= (mp_inf - mp) / tau_mp\n"
  "	m' = (m_inf - m) / tau_m\n"
  "	h' = (h_inf - h) / tau_h\n"
  "	s' = (s_inf - s) / tau_s\n"
  "	hs'= (hs_inf - hs) / tau_hs\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  ":\n"
  ":	Q10 adjustment\n"
  ":\n"
  "\n"
  "	q10_1 = 2.2 ^ ((celsius-20)/ 10 )\n"
  "	q10_2 = 2.9 ^ ((celsius-20)/ 10 )\n"
  "	q10_3 = 3.0 ^ ((celsius-36)/ 10 )\n"
  "\n"
  "	evaluate_fct(v)\n"
  "	mp = mp_inf\n"
  "	m = m_inf\n"
  "	h = h_inf\n"
  "	s = s_inf\n"
  "	hs = hs_inf\n"
  "}\n"
  "\n"
  "PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2\n"
  "\n"
  "	a = q10_1*vtrap1(v)\n"
  "	b = q10_1*vtrap2(v)\n"
  "	tau_mp = 1 / (a + b)\n"
  "	mp_inf = a / (a + b)\n"
  "\n"
  "	a = q10_1*vtrap6(v)\n"
  "	b = q10_1*vtrap7(v)\n"
  "	tau_m = 1 / (a + b)\n"
  "	m_inf = a / (a + b)\n"
  "\n"
  "	a = q10_2*vtrap8(v)\n"
  "	b = q10_2*bhA / (1 + Exp(-(v+bhB)/bhC))\n"
  "	tau_h = 1 / (a + b)\n"
  "	h_inf = a / (a + b)\n"
  "\n"
  "	a = q10_2*vtrap8(v+HSbiasV)\n"
  "	b = q10_2*bhA / (1 + Exp(-(v+HSbiasV+bhB)/bhC))\n"
  "	tau_hs = SIF / (a + b)\n"
  "	hs_inf = a / (a + b)\n"
  "\n"
  "	v2 = v - vtraub : convert to traub convention\n"
  "\n"
  "	a = q10_3*asA / (Exp((v2+asB)/asC) + 1)\n"
  "	b = q10_3*bsA / (Exp((v2+bsB)/bsC) + 1)\n"
  "	tau_s = 1 / (a + b)\n"
  "	s_inf = a / (a + b)\n"
  "\n"
  "}\n"
  "\n"
  "FUNCTION vtrap(x) {\n"
  "	if (x < -50) {\n"
  "		vtrap = 0\n"
  "	}else{\n"
  "		vtrap = bsA / (Exp((x+bsB)/bsC) + 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION vtrap1(x) {\n"
  "	if (fabs((x+ampB)/ampC) < 1e-6) {\n"
  "		vtrap1 = ampA*ampC\n"
  "	}else{\n"
  "		vtrap1 = (ampA*(x+ampB)) / (1 - Exp(-(x+ampB)/ampC))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION vtrap2(x) {\n"
  "	if (fabs((x+bmpB)/bmpC) < 1e-6) {\n"
  "		vtrap2 = bmpA*bmpC : Ted Carnevale minus sign bug fix\n"
  "	}else{\n"
  "		vtrap2 = (bmpA*(-(x+bmpB))) / (1 - Exp((x+bmpB)/bmpC))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION vtrap6(x) {\n"
  "	if (fabs((x+amB)/amC) < 1e-6) {\n"
  "		vtrap6 = amA*amC\n"
  "	}else{\n"
  "		vtrap6 = (amA*(x+amB)) / (1 - Exp(-(x+amB)/amC))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION vtrap7(x) {\n"
  "	if (fabs((x+bmB)/bmC) < 1e-6) {\n"
  "		vtrap7 = bmA*bmC : Ted Carnevale minus sign bug fix\n"
  "	}else{\n"
  "		vtrap7 = (bmA*(-(x+bmB))) / (1 - Exp((x+bmB)/bmC))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION vtrap8(x) {\n"
  "	if (fabs((x+ahB)/ahC) < 1e-6) {\n"
  "		vtrap8 = ahA*ahC : Ted Carnevale minus sign bug fix\n"
  "	}else{\n"
  "		vtrap8 = (ahA*(-(x+ahB))) / (1 - Exp((x+ahB)/ahC))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION Exp(x) {\n"
  "	if (x < -100) {\n"
  "		Exp = 0\n"
  "	}else{\n"
  "		Exp = exp(x)\n"
  "	}\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif