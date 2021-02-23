//#include "TClassHolder.h"
#include <iostream>
#include <fstream>
#include <string>
//#include <gsl/gsl_rng.h>
//#include "globals.hh"
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <TGraphErrors.h>
#include <stdio.h>
#include <complex>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include "TFitter.h"
#include <TTree.h>
//#include "Randomize.hh"
#include <TRandom1.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TVector3.h>
#include "TLorentzVectorC.h"
#include "Math/GSLMinimizer.h"
#include "Math/Functor.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2D.h"
#include "TH1D.h"
#include "malloc.h"
#include "TGraph2D.h"
#include <TGenPhaseSpace.h>

using namespace std;
#ifndef Cmd3Generator2pi2pi0_ke_H
#define Cmd3Generator2pi2pi0_ke_H
#define mK 493.677
#define mKn 497.614
#define alpha1 0.00729735253
#define  pi 3.141592
//#include "Cmd3CLHEPGenerator.h"
//#include "Cmd3PrimaryGenerator.h"
//#include <TLorentzVector.h>
#include "TLorentzVectorC.h"
//#include "Cmd3CLHEPGenerator.h"

#define me 0.000511
#define Mpi 0.139570
#define mpiz 0.134976
#define Momega 0.78265
#define Womega 0.00849
#define Ma1 1.250
#define Wa1 0.40
#define Ma1_2 1.640
#define Wa1_2 0.20
#define Ma2 1.318
#define Wa2 0.150
#define Mrho 0.77526
#define Wrho 0.1491
#define MK 0.493
#define Msigma 0.5
#define Wsigma 0.564
#define Mf0 0.98
#define Wf0 0.07
#define Wf2 0.25
#define Mf2 1.175
#define Mh1 1.17
#define Wh1 0.36
#define Mpi1300 1.30
#define Wpi1300 0.36
#define pi 3.141592
#define Convers 0.389379292E9

double newfactor = 0.;

class Cmd3Generator2pi2pi0_ke// public Cmd3CLHEPGenerator
{

        public:
              Cmd3Generator2pi2pi0_ke();
              virtual ~Cmd3Generator2pi2pi0_ke() {};
              double cross_a1_sigmapi_pi(double *e, double *par);
              double cross_a1_rhopi_pi(double *e, double *par);
              double cross_omega_pi0(double *e, double *par);
              //void Start_Of_Run(CmdParam& p);
              //HepMC::GenEvent* Event(int NumEvent=0);
	      void change(int num, double value);
              void fill_parameters();
	      double majoranta = 200.;
              string results = "results.dat";
	      double cross_norm(double *e, double *par);
              string scross_norm = "cross_norm_a1pi.dat";
              double levicivita(int A, int B, int C, int D);
              double TLorentzVectorC_index(TLorentzVectorC P, int index);
              double svertka_leviciv(TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC sver_vec_leviciv(TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);

              
              TLorentzVectorC  H_ph_sp(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_ph_sp_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_ph_sp(double *e, double *par);

              double factor_phsp();
              double factor_phsp_ph();

              TLorentzVectorC  H_a1_sigmapi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_a1_sigmapi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double factor_a1pi_sigmapi();
              double factor_a1pi_sigmapi_ph();

              TLorentzVectorC  H_a1_rhopi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_a1_rhopi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double factor_a1pi_rhopi();
              double factor_a1pi_rhopi_ph();
              
              double cross_a1_rhopi_pi_II(double *e, double *par);
              TLorentzVectorC  H_a1_rhopi_pi_beforepi0sim_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_a1_rhopi_pi_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);              
              double factor_a1pi_rhopi_II();
              double factor_a1pi_rhopi_ph_II();
	      TLorentzVectorC H_a1_rhopi_pi_c_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);


              TLorentzVectorC  H_omega_pi0_init(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_omega_pi0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC H_omega_pi0_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double factor_ompi0();
              double factor_ompi0_ph();


              TLorentzVectorC  H_rho_f0_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_rho_f0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_rho_f0(double *e, double *par);
              double factor_rho_f0();
              double factor_rho_f0_ph();

              TLorentzVectorC  H_rho_f0_inter_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_rho_f0_inter(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_rho_f0_inter(double *e, double *par);
              double factor_rho_f0_inter();
              double factor_rho_f0_inter_ph();

              TLorentzVectorC  H_rho_sigma(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC  H_rho_sigma_beforesim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_rho_sigma(double *e, double *par);
              double factor_rho_sigma();
              double factor_rho_sigma_ph();


              TLorentzVectorC  H_rhop_rhom_beforpi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_rhop_rhom(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_rhop_rhom(double *e, double *par);
              double factor_rhop_rhom();
              double factor_rhop_rhom_ph();

	      double cross_rhop_rhom_II(double *e, double *par);
 	      TLorentzVectorC  H_rhop_rhom_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);

              double mod(TLorentzVectorC P){return P.X().real()*P.X().real() + P.X().imag()*P.X().imag() + P.Y().real()*P.Y().real() + P.Y().imag()*P.Y().imag();}
              TLorentzVectorC  H_a2_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_a2_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double cross_a2_pi(double *e, double *par);
              double factor_a2pi();
              double factor_a2pi_ph();

              double cross_h1_rhopi(double *e, double *par);
              TLorentzVectorC  H_h1_rhopi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_h1_rhopi_bef(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              TLorentzVectorC  H_h1_rhopi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double factor_h1_rhopi();
              double factor_h1_rhopi_ph();
              double cross_rho_f2(double *e, double *par);
	      TLorentzVectorC H_rho_f2(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              void readd(double minenergy, double maxenergy, string );
              void readd_bkgr(double minenergy, double maxenergy, string);
              void readd_c(double minenergy, double maxenergy, string);
	      TLorentzVectorC H_a1_rhopi_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_a1_sigmapi_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_rho_f0_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_rho_f0_inter_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_rho_sigma_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_ph_sp_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC H_a2_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      TLorentzVectorC matrix(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double matrix_squared(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double matrix_squared1(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4);
              void find_majoranta0(double E);
              void find_min0(double E);
              double cross_section();
              void fill_sample();
              double generate(int);
              double matrix_squared_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
              double matrix_squared1_c(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4);

	      double matrix_squared_m000(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      double matrix_squared_mmp0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4);
	      double matrix_squared_m0001(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4);
	      double matrix_squared_mmp01(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4);

	      double factor_ampl1(double enrgy, int number);
	      double factor_ampl(double enrgy, int number);

              void fill0_parameters();
	      void fill_omega_pi0();
              void fill_a1pi_rhopi();
              void fill_a1pi_sigmapi();
              void fill_rhof0();
              void fill_rhosigma();
              void fill_phsp();
              void fill_rhoprhom();
              void fill_h1pi();
              void fill_a2pi();
              double cross_s[100];
              double dcross_s[100];

              double Energy = 1.02;
              int isISR = 0;

              double omegapi0_ = 1.0;
              double fompi0 = 0.;
              double dompi0 = 0.;
              double omegapi0phase_ = 0.;
              
              double a1pi_rhopi_ =  0.;
              double a1pi_rhopi_phase = 0;
              double da1pi_rhopi = 0.;
              double fa1pi_rhopi = 0.;
              
              double a1pi_sigmapi_ = 0.0;
              double a1pi_sigmapi_phase = 0.0;
              double da1pi_sigmapi = 0.0;
              double fa1pi_sigmapi = 0.0;
              
              double rhoprhom_ = 0.;
              double rhoprhom_phase = 0.;
              double rhoprhom_II = 0.;
              double rhoprhom_phase_II = 0.;
              double drhoprhom = 0.;
              double frhoprhom = 0.;              
              
              double rhof0_ = 0.;
              double rhof0phase_ = 0.;
              double drhof0 = 0.;
              double frhof0 = 0.;              
              
              double rhosigma_ = 0.;
              double rhosigmaphase_ = 0.;
              double drhosigma = 0.;
              double frhosigma = 0.;              

              double rhof0_inter = 0.;
              double rhof0phase_inter = 0.;
              double rhof2 = 0.;
              double rhof2_phase = 0.;  

              double h1pi = 0;
              double h1pi_ph = 0;
              double dh1pi = 0;
              double fh1pi = 0;
              
              double a2pi_ = 0.;
              double a2pi_ph = 0.;
              double a1pi_rhopi_II =  0.;
              double a1pi_rhopi_phase_II = 0;
              double phsp_ = 0.;
              double phspphase_ = 0.;

              
              double cross_section_c();
              void fill_sample_c();


	      void print();
              double LIKE();
 	      double cross(double E);
	      double visible_cross_section(double E);
	      double generate_photon_energy(double E);
	      double generate_photon_angle(double E, double x);
	      int Mode2pi2pi0 = 1;
	      double cross2pi2pi0(double E);
	      double cross4pi(double E);

};
#endif


double Likexp[1000000];
double Likmc[1000000];
double Likexp_c[1000000];
double Likmc_c[1000000];

double OmgWidth(double);
double RhoWidth(double);
double A1Width(double *e, double *par);

complex<double>  propagator_omega(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Momega*Momega,sqrt(s)*OmgWidth(s));
  return numerator/denominator;
}

complex<double>  propagator_omegaprime(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-1.42*1.42,sqrt(s)*0.2);
  return numerator/denominator;
}


complex<double>  propagator_a1(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  double ene = sqrt(s);
  double factor = A1Width(&ene, &ene);
  if(sqrt(s) <= 2.*Mpi)factor = 1;
  complex<double> denominator(s-Ma1*Ma1,Ma1*Wa1*factor);
  return numerator/denominator;
}

complex<double>  propagator_a1_2(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Ma1_2*Ma1_2,Ma1_2*Wa1_2);
  return numerator/denominator;
}

complex<double>  propagator_a2(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Ma2*Ma2,sqrt(s)*Wa2);
  return numerator/denominator;
}

complex<double>  propagator_rho(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Mrho*Mrho,Mrho*RhoWidth(s));
  return numerator/denominator;
}


complex<double> M_sigma(0.572,0);
complex<double> M_f0(0.94,0);
complex<double> consttt = complex<double>(0.01,0.);
complex<double> g_f0pipi = complex<double>(-2.4,0.);

complex<double> P_pipi(complex<double> m, double g_Spipi){
	return 1.5*g_Spipi*g_Spipi*sqrt(1.-4.*Mpi*Mpi/m/m)/16./3.141592*(complex<double>(0.,1.) + log((1.-sqrt(1.-4*Mpi*Mpi/m/m))/(1.+sqrt(1.+4*Mpi*Mpi/m/m)))/3.141592);
}

complex<double> P_KK(complex<double> m, double g_SKK){
        complex<double> vklad1 = 2.*g_SKK*g_SKK*sqrt(1.-4.*0.494*0.494/m/m)/16./3.141592*(complex<double>(0.,1.) + log((1.-sqrt(1.-4.*0.494*0.494/m/m))/(1.+sqrt(1.+4.*0.494*0.494/m/m)))/3.141592);
	complex<double> vklad2 = -2.*g_SKK*g_SKK*sqrt(4.*0.494*0.494/m/m-1.)/16./3.141592*(1. - 2./3.1415*TMath::ATan(sqrt(4.*0.494*0.494/m.real()/m.real()-1.)));
	if(m.real() >= 2.*0.494)return vklad1;
	else return vklad2;
}

complex<double> P_etaeta(complex<double> m, double g_Setaeta){
	complex<double> vklad1 =  g_Setaeta*g_Setaeta*sqrt(1.-4.*0.548*0.548/m/m)/16./3.141592*(complex<double>(0.,1.) + log((1.-sqrt(1.-4*0.548*0.548/m/m))/(1.+sqrt(1.+4*0.548*0.548/m/m)))/3.141592);
	complex<double> vklad2 = -g_Setaeta*g_Setaeta*sqrt(4.*0.548*0.548/m/m-1.)/16./3.141592*(1.-2./3.141592*TMath::ATan(4*0.548*0.548/m.real()/m.real()-1));
	if(m.real() >= 2.*0.548)return vklad1;
	else return vklad2;
}

complex<double> P_sigma(complex<double> m){
	return P_pipi(m,3.)+P_KK(m,0.46)+P_etaeta(m,3./sqrt(2.));
}
complex<double> P_f0(complex<double> m){
	return P_pipi(m,g_f0pipi.real())+P_KK(m,7.)+P_etaeta(m,7.);
}
complex<double> P_sigmaf0(double m1){
        complex<double> m(m1,0);
	return g_f0pipi.real()*P_pipi(m,3.)/3.+7.*P_KK(m,0.46)/0.46 + 7.*P_etaeta(m,3./sqrt(2.))/3.*sqrt(2.) + consttt;
}
complex<double> D_sigma(double m){
	return M_sigma*M_sigma-m*m + P_sigma(M_sigma).real() - P_sigma(m);
}
complex<double> D_f0(double m){
	return M_f0*M_f0-m*m + P_f0(M_f0).real() - P_f0(m);
}


complex<double> propagator_sigma_new(double m,complex<double> *par){
  if(m <= 2.*Mpi)return complex<double>(0.,0.);
  complex<double> Delta = D_sigma(m)*D_f0(m) - P_sigmaf0(m)*P_sigmaf0(m);
  double rhopipi = sqrt(1.-4*Mpi*Mpi/m/m);
  complex<double> g_sigma_2pi0 = complex<double>(3.,0.)/sqrt(2);
  complex<double> g_gg_sigma = par[0];
  complex<double> g_f0_2pi0 = g_f0pipi/sqrt(2);
  complex<double> g_gg_f0 = par[2];
  return rhopipi/Delta*(g_gg_f0*g_f0_2pi0*D_sigma(m) + (g_gg_sigma*g_f0_2pi0+g_gg_f0*g_sigma_2pi0)*P_sigmaf0(m) + g_gg_sigma*g_sigma_2pi0*D_f0(m));

}


complex<double>  propagator_sigma(complex<double> s1){
 
  double gsigmapipi = 80;
  double s = sqrt(s1.real());
  double phase = 0.;
  complex<double> numerator(cos(phase),sin(phase));


  double imSigma = -3./32./3.1415*gsigmapipi*gsigmapipi*Mpi*Mpi*sqrt(1.-4.*Mpi*Mpi/s);
  double reSigma = gsigmapipi*gsigmapipi/64./3.1415/3.1415/sqrt(s*(s-4.*Mpi*Mpi));
  double fff = (s-12.*Mpi*Mpi)*sqrt(s*(s-4.*Mpi*Mpi)) + 6.*Mpi*Mpi*(s-4.*Mpi*Mpi)*log((sqrt(s*(s-4.*Mpi*Mpi)) + s)/2./Mpi/Mpi-1.);
  reSigma = reSigma*fff;
  return numerator/complex<double>(s - Msigma*Msigma -reSigma,-imSigma)*10.;

  double factor = sqrt(fabs(1.-4.*Mpi*Mpi/s))/sqrt(1.-4.*Mpi*Mpi/Msigma/Msigma);
  if(sqrt(s) <= 2.*Mpi)factor = - factor;
  double factor1 = 1.;
  double sa = Mpi*Mpi/2.;
  complex<double> mass(Msigma,0.);
  factor1 = 1.;//fabs(s - sa)/(0.93*0.93 - sa)*(0.58 + 1.66*s)*exp(-(s-0.93*0.93)/1.08);
  complex<double>  denom = mass*mass - s - mass*Wsigma*factor*complex<double>(0.,1.)*factor1;
  complex<double> denominator(denom.real(),denom.imag());
  return numerator/denominator;
}

complex<double>  propagator_f0_inter(complex<double> s1){
  double m = sqrt(s1.real());
  complex<double> parameters[] = {complex<double>(0.,0.), complex<double>(4.35778e+01, 3.73122e+00), complex<double>(0., 0.)};
  consttt = complex<double> (-6.09049e-01, 2.78206e-01);
  //return propagator_sigma_new(m,parameters);

  double s = sqrt(s1.real());
  double phase = 0.;
  complex<double> numerator(cos(phase),sin(phase));
  double factor = Msigma/sqrt(s)*sqrt(1.-4.*Mpi*Mpi/s)/sqrt(1.-4.*Mpi*Mpi/Msigma/Msigma);
  if(sqrt(s) <= 2.*Mpi)factor = 1;
  double factor1 = 1.;
  double sa = Mpi*Mpi/2.;
  complex<double> mass(1.4,0);
  factor1 = fabs(s - sa)/(0.93*0.93 - sa)*(0.58 + 1.66*s)*exp(-(s-0.93*0.93)/1.08);
  complex<double>  denom = mass*mass - s + mass*0.3*factor*complex<double>(0.,1.)*factor1;
  complex<double> denominator(denom.real(),denom.imag());
  return numerator/denominator;

}


complex<double>  propagator_f0(complex<double> s1){
  double m = sqrt(s1.real());
  complex<double> parameters[] = {complex<double>(4.47772e+01,0.00000e+00), complex<double>(0., 0.), complex<double>(0., 0.)};
  consttt = complex<double> (-6.09049e-01, 2.78206e-01);
 // return propagator_sigma_new(m,parameters);

  // flatte formula https://arxiv.org/pdf/1102.0206.pdf 
  double s = s1.real();
  double phase = 0.;
  complex<double> numerator(cos(phase),sin(phase));
  complex<double> factor_rhopipi = 2.*sqrt(fabs(s/4.-Mpi*Mpi))/sqrt(s)*complex<double>(0.,1.);
  complex<double> factor_rhoKK = 2.*sqrt(fabs(s/4.-MK*MK))/sqrt(s)*complex<double>(0.,1.);
  if(sqrt(s) <= 2.*Mpi)factor_rhopipi = -2.*sqrt(fabs(Mpi*Mpi-s/4.))/sqrt(s);
  if(sqrt(s) <= 2.*MK)factor_rhoKK = -2.*sqrt(fabs(MK*MK-s/4.))/sqrt(s);
  complex<double> mass(Mf0,0);
  complex<double>  denom = Mf0*Mf0 - s - 0.165*(factor_rhopipi + 4.21*factor_rhoKK);
  //factor_rhopipi = sqrt(fabs(s/4.-Mpi*Mpi))/sqrt(fabs(Mf0*Mf0/4.-Mpi*Mpi));
  //if(sqrt(s) <= 2.*MK)factor_rhopipi = sqrt(fabs(s/4.-MK*MK))/sqrt(fabs(Mf0*Mf0/4.-MK*MK));
  //denom = mass*mass - s - sqrt(s)*Wf0*factor_rhopipi*complex<double>(0.,1.);
  complex<double> denominator(denom.real(),denom.imag());
  return numerator/denominator;
}





complex<double>  propagator_f2(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  double factor = sqrt(1.-4.*Mpi*Mpi/s)/sqrt(1.-4.*Mpi*Mpi/Mf2/Mf2);
  if(sqrt(s) <= 2.*Mpi)factor = 1;
  complex<double> denominator(s-Mf2*Mf2,sqrt(s)*Wf2*factor);
  return numerator/denominator;
}

complex<double>  propagator_h1(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Mh1*Mh1,sqrt(s)*Wh1);
  return numerator/denominator;
}



complex<double>  propagator_pi1300(complex<double> s1){
  double s = s1.real();
  complex<double> numerator(1.0,0.0);
  complex<double> denominator(s-Mpi1300*Mpi1300,sqrt(s)*Wpi1300);
  return numerator/denominator;
}






Double_t RhoWidth(Double_t s)
{

  return Wrho/sqrt(s)*Mrho*TMath::Power(fabs(s/4.0-Mpi*Mpi)/(Mrho*Mrho/4.0-Mpi*Mpi),1.5);
  s = s*1000000.;
  Double_t MRho = 775.49;
  Double_t WRho = 149.10;

  Double_t BrKC = 0.489;
  Double_t BrK0 = 0.342;

  Double_t BrPi0Gamma = 0.0006;
  Double_t BrEtaGamma = 0.0003;

  Double_t MPhi = 1019.456;

  Double_t PhiWidth(Double_t);
  Double_t PKPKM(Double_t);
  Double_t PKLKS(Double_t);
  Double_t QPGamma(Int_t, Double_t);

  Double_t Fval = WRho/s*MRho*MRho*TMath::Power((s/4.0-Mpi*Mpi)/(MRho*MRho/4.0-Mpi*Mpi),1.5)
        + 0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrKC*PKPKM(s)/PKPKM(MPhi*MPhi)/s
                + 0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrK0*PKLKS(s)/PKLKS(MPhi*MPhi)/s
                + WRho*BrPi0Gamma*QPGamma(0,s)/QPGamma(0,MRho*MRho)
                + WRho*BrEtaGamma*QPGamma(1,s)/QPGamma(1,MRho*MRho);
  return Fval/1000.;
}

Double_t OmgWidth(Double_t s)
{

  s = s*1000000.;
  Double_t MOmg = 782.65;
  Double_t WOmg =  8.44;
  Double_t MPhi = 1019.456;
  Double_t WPhi = 4.26;

  Double_t MPi  = 139.56995;
  Double_t MPi0 = 134.9766;
  Double_t BrKC = 0.489;
  Double_t BrK0 = 0.342;
  Double_t BrEtaGamma = 0.00065;
  Double_t BrOmg3Pi   = 0.891;
  Double_t BrOmg2Pi   = 0.017;
  Double_t BrOmgPiG   = 0.087;

  Double_t PhiWidth(Double_t);
  Double_t PKPKM(Double_t);
  Double_t PKLKS(Double_t);
  Double_t QPGamma(Int_t, Double_t);
  Double_t FAS_ASPO(Double_t);

  Double_t Temp1, Temp2, Temp3, Temp4, Temp5, Temp6;
  Double_t X;

  Double_t TwoE;

  TwoE = TMath::Sqrt(s)/1000.0;

  if ( TMath::Sqrt(s) < 750.0 ) {
    Temp1 = 1.0E-6;
  } else {
    X = TMath::Sqrt(s)/1000.0 - 0.75;
    Temp1 = 1.0469 + 10.457*X + 53.009*X*X+591.66*X*X*X + 1487.7*X*X*X*X-7623.2*X*X*X*X*X+6949.5*X*X*X*X*X*X;
  }
  X = MOmg/1000.0 - 0.75;
  Temp2 = 1.0469 + 10.457*X + 53.009*X*X+591.66*X*X*X + 1487.7*X*X*X*X-7623.2*X*X*X*X*X+6949.5*X*X*X*X*X*X;

  Temp3 = TMath::Power(0.5*(1.0-MPi0*MPi0/s)*TMath::Sqrt(s),3.0);
  Temp4 = TMath::Power(0.5*(1.0-MPi0*MPi0/MOmg/MOmg)*TMath::Sqrt(MOmg*MOmg),3.0);

  Temp5 = TMath::Power((s/4.0-MPi*MPi),1.5);
  Temp6 = TMath::Power((MOmg*MOmg/4.0-MPi*MPi),1.5);

  Double_t Fval = WOmg*(BrOmg3Pi*FAS_ASPO(TwoE)/FAS_ASPO(MOmg/1000.0) +
            BrOmgPiG*QPGamma(0,s)/QPGamma(0,MOmg*MOmg)    +
            BrOmg2Pi*Temp5*MOmg*MOmg/(Temp6*s))           +
                        0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrKC*PKPKM(s)/PKPKM(MPhi*MPhi)/s +
                        0.5*PhiWidth(MPhi*MPhi)*MPhi*MPhi*BrK0*PKLKS(s)/PKLKS(MPhi*MPhi)/s +
                        WOmg*BrEtaGamma*QPGamma(1,s)/QPGamma(1,MOmg*MOmg);

  return Fval/1000.;
}

Double_t PhiWidth(Double_t s)
{
  Double_t MPhi = 1019.456;
  Double_t WPhi = 4.26;
  Double_t MPi  = 139.569995;
  Double_t MPi2;

  Double_t Br[4] = {0.489, 0.342, 0.1532, 0.01309};
  Double_t Br2Pi = 0.000073;
  Double_t BrPi0Gamma = 0.00124;

  Double_t MPhi2 = 1019.456*1019.456;

  Double_t PKPKM(Double_t);
  Double_t PKLKS(Double_t);
  Double_t QPGamma(Int_t, Double_t);

  Double_t FAS_ASPO(Double_t);

  Double_t W1, W2;
  Double_t TwoE;

  TwoE = TMath::Sqrt(s)/1000.0;
  W1   = FAS_ASPO(TwoE);
  W2   = FAS_ASPO(MPhi/1000.0);

  MPi2 = MPi*MPi;


  Double_t Fval = MPhi2*WPhi*(Br[0]*PKPKM(s)/PKPKM(MPhi2) +
                                  Br[1]*PKLKS(s)/PKLKS(MPhi2) +
                      Br2Pi*TMath::Power((s/4.0-MPi2)/(MPhi2/4.0-MPi2),1.5))/s +
                            WPhi*(Br[2]*W1/W2+BrPi0Gamma*QPGamma(0,s)/QPGamma(0,MPhi2) +
                                  Br[3]*QPGamma(1,s)/QPGamma(1,MPhi2));

  return Fval;
}

Double_t PKPKM(Double_t s)
{
  Double_t MKC = 493.677;

  if ( s < 4.0*MKC*MKC ) return 0.0;

  return TMath::Power((s/4.0-MKC*MKC),1.5);
}

Double_t PKLKS(Double_t s)
{
  Double_t MK0 = 497.672;

  if ( s < 4.0*MK0*MK0 ) return 0.0;

  return TMath::Power((s/4.0-MK0*MK0),1.5);
}

Double_t QPGamma(Int_t Mode, Double_t s)
{
  Double_t PMass;

  if ( Mode == 0 ) PMass = 134.9766;
  if ( Mode == 1 ) PMass = 547.3;

  if ( s < PMass*PMass ) return 0.0;

  Double_t Fval =  TMath::Power(0.5*(1.0-PMass*PMass/s),3)*TMath::Power(s,1.5);

  return Fval;
}

Double_t Z(Double_t s)
{
  Double_t V,e;

  if ( s < 4.0*mK*mK ) return 0.0;
  e = sqrt(s/4.);
  V = sqrt(1.-mK*mK/e/e);

  return alpha1*TMath::Pi()/V/(1.-exp(-alpha1*TMath::Pi()/V))*(1.+alpha1*alpha1/4./V/V);
}


Double_t FAS_ASPO(Double_t TwoE)
{

  Double_t Polinom(Int_t, Double_t, Double_t*);

  const Int_t K_max = 4;
  Int_t I_mark = 100;
  Int_t Mode = 1;

  bool Calc = true;

  Double_t A[4] = {-6.12394622, 25.0341405, -34.1311022, 15.5413717};
  Double_t B[4] = {5.29354148, -7.90990714, -2.26007613, 5.21453902};
  double CC[4] = {-0.5643611, 2.69953793, -4.32966739, 2.33116866};
  Double_t D[4] = {-0.0548334238,  0.31600391, -0.609523718, 0.393667808};

  Double_t A1[4] = {-4.91624401 , 19.8655606, -26.9136128 , 12.2412286};
  Double_t D1[4] = {-0.00794774189,0.0522269164,-0.114526409 , 0.0838126536};

  Double_t LM[6] = {1.1,0.875,0.75,0.62,0.52,0.46};

  Double_t POL, TEMP1, TEMP2;

  if ( Calc ) {
    TEMP2 = D[1]*LM[4]+D[2]*LM[4]*LM[4]+D[3]*LM[4]*LM[4]*LM[4];
    TEMP1 = D1[0]+D1[1]*LM[4]+D1[2]*LM[4]*LM[4]+D1[3]*LM[4]*LM[4]*LM[4];
    D[0]  = TEMP1 - TEMP2;

    TEMP2 = CC[1]*LM[3]+CC[2]*LM[3]*LM[3]+CC[3]*LM[3]*LM[3]*LM[3];
    TEMP1 = D[0]+D[1]*LM[3]+D[2]*LM[3]*LM[3]+D[3]*LM[3]*LM[3]*LM[3];
    CC[0]  = TEMP1 - TEMP2;

    TEMP2 = A1[1]*LM[2]+A1[2]*LM[2]*LM[2]+A1[3]*LM[2]*LM[2]*LM[2];
    TEMP1 = CC[0]+CC[1]*LM[2]+CC[2]*LM[2]*LM[2]+CC[3]*LM[2]*LM[2]*LM[2];
    A1[0] = TEMP1 - TEMP2;

    TEMP1 = A1[0]+A1[1]*LM[1]+A1[2]*LM[1]*LM[1]+A1[3]*LM[1]*LM[1]*LM[1];
    TEMP2 = A[1]*LM[1]+A[2]*LM[1]*LM[1]+A[3]*LM[1]*LM[1]*LM[1];
    A[0]  = TEMP1 - TEMP2;

    TEMP1 = A[0] + A[1]*LM[0]+A[2]*LM[0]*LM[0]+A[3]*LM[0]*LM[0]*LM[0];
    TEMP2 = B[1]*LM[0]+B[2]*LM[0]*LM[0]+B[3]*LM[0]*LM[0]*LM[0];
    B[0]  = TEMP1 - TEMP2;

    Calc = false;
  }
  POL = 0.0;
  if ( TwoE >= LM[0] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*B[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*B[i];
    }
  }

  if ( TwoE >= LM[1] && TwoE < LM[0] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*A[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*A[i];
    }
  }

  if ( TwoE >= LM[2] && TwoE < LM[1] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*A1[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*A1[i];
    }
  }

  if ( TwoE >= LM[3] && TwoE < LM[2] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*CC[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*CC[i];
    }
  }

  if ( TwoE >= LM[4] && TwoE < LM[3] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*D[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*D[i];
    }
  }

  if ( TwoE >= LM[5] && TwoE < LM[4] ) {
    for( Int_t i = 0; i < K_max; i++ ) {
      if ( Mode == 1 ) POL += TMath::Power(TwoE,(Double_t)i)*D1[i];
      if ( Mode == 2 ) POL += TMath::Power(TwoE,(Double_t)(i-1))*i*D1[i];
    }
  }
  if ( TwoE < LM[5] ) {
    if ( Mode == 1 ) POL = D1[0]+D1[1]*LM[5]+D1[2]*LM[5]*LM[5]+D1[3]*LM[5]*LM[5]*LM[5];
    if ( Mode == 2 ) POL = 3.0*LM[5]*LM[5]*D1[3]+2.0*LM[5]*D1[2]+D1[1];
  }


  Double_t Fval = (POL/0.393728*(0.00749/0.0361478));

  return Fval;
}

double A1Width(double *e, double *par){
    double norm[] = {0.00157213, 0.00341551, 0.00739029, 0.0162539, 0.037892, 0.0943931, 0.216952, 0.370924, 0.516185, 0.641809, 0.749868, 0.83974, 0.930229, 1.00788, 1.08098, 1.15526, 1.23872, 1.30904, 1.36955, 1.44756, 1.50457, 1.58155, 1.64559, 1.73735, 1.79019, 1.86764, 1.94434, 2.01752};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor*Ma1/e[0];
}









