// Generator 2pi2pi0
// Everything should be in GeV
#include "Cmd3Generator2pi2pi0_ke.h"
#include <sys/types.h>
#include <dirent.h>
#include<stdio.h>
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TArrow.h"
using namespace std;

double minmajoranta = 200.;//1.5e06;
double SCALE = 0;
//double RATIO = 0.;
bool TOY = true;
string if_name_AA_data = "../4pisel/histograms/data/4pi_data.dat";
string if_name_AA_data_c = "../2pip2pim/histograms/data/4pic_data.dat";

/*
parameter number & amplitude #   
0 - lambda
1 - omega pi0
2 - a1pi (rho pi)
3 - a1pi (rho pi) phase
4 - a1pi (sigma pi)
5 - a1pi (sigma pi) phase
6 - rho f0
7 - rho f0 phase
8 - rho+ rho-
9 - rho+ rho- phase
10 - rho sigma
11 - rho sigma phase
12 - ph sp
13 - ph sp phase
14 - a2 pi
15 - a2 pi phase
16 - h1 pi
17 - h1 pi phase
if_name_AA_data_c = "../2pip2pim/histograms/v1/4picmc_990_16304.root.dat";
if_name_AA_data = "/home/eakozyrev/diskD/4pi_new/4pisel/histograms/2pi2pi0_v1/mc_990_16304.root.dat";

*/

Cmd3Generator2pi2pi0_ke Part_ana;

//typedef HepMC::FourVector HepLorentzVector;

//--------------------------------------------------------------------------------------------------------------------------------------------------
//static TClassHolder<Cmd3Generator2pi2pi0_ke, Cmd3PrimaryGenerator> fP3i2pi0_ke_gen("P3i2pi0_ke_gen");
//typedef HepMC::FourVector HepLorentzVector;
//**************************************************************************************************************************************************

Cmd3Generator2pi2pi0_ke::Cmd3Generator2pi2pi0_ke(){}

double Cmd3Generator2pi2pi0_ke::factor_ampl(double enrgy, int number){
ifstream stream1(("/media/eakozyrev/diskD/4pi_new/mygenerator/data/"+results).c_str());
TString btv;
double energy;
double energy0;
double data[100];
double data0[100];
double res = -10000;
bool check = false;
int npar = number;
while(stream1.eof()==0){
        int i = -1;
        stream1 >> energy;
        if(stream1.eof()==1)break;
        while(btv!="end"){
            if(i>-1){
                data[i] = btv.Atof();
            }
            stream1 >> btv;
            if(stream1.eof()==1)break;
            i++;
        }
        btv = "as";
        if(enrgy <= energy){
            res = data[npar];
            if(check == true){
		res = data[npar]*(enrgy - energy0)/(energy - energy0) + data0[npar]*(- enrgy + energy)/(energy - energy0);
		if(number%4==2 && enrgy != energy)res = data0[npar];
	    }
            break;
        }
        check = true;
        for(int i = 0; i < 100; i++){data0[i] = data[i];}
        energy0 = energy;
  }
  if(res < -1000)res = data[npar];
  return res;
}


double Cmd3Generator2pi2pi0_ke::factor_ampl1(double enrgy, int number){

double res = 0;
double width = 0.09;
double step = 0.00005;
double en = enrgy - width/2.;
while(en < (enrgy + width/2.)){

	res = res+step*factor_ampl(en, number);
	en = en + step;
}

return res/width;
}

double Cmd3Generator2pi2pi0_ke::cross_norm(double *e, double *par){

    ifstream streamd(("/home/eakozyrev/diskD/4pi_new/mygenerator/data/new/"+scross_norm).c_str());
    int i = 0;double norm[1000];
    while(streamd.eof()==0){streamd >> norm[i]; i++; }
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

void Cmd3Generator2pi2pi0_ke::change(int num, double value){

	if(num == 1)omegapi0_ = value;
	if(num == 2){
	  if(value == -1000)a1pi_rhopi_ =a1pi_rhopi_+factor_ampl(Energy,5);
	  else a1pi_rhopi_ = value;
	}
	if(num == 3){
	  if(value == -1000)a1pi_rhopi_phase = a1pi_rhopi_phase +factor_ampl(Energy,7);
	  else a1pi_rhopi_phase = value;
	}
	if(num == 4){
	  if(value == -1000)a1pi_sigmapi_ = a1pi_sigmapi_ + factor_ampl(Energy,9);
	  else a1pi_sigmapi_ = value;
	}
	if(num == 5){
	  if(value== -1000)a1pi_sigmapi_phase = a1pi_sigmapi_phase  + factor_ampl(Energy,11);	  
	  else a1pi_sigmapi_phase = value;
	}
	if(num == 6){
	  if(value == -1000)rhof0_ = rhof0_ +  factor_ampl(Energy,13);	 
	  else rhof0_ = value;
	}
	if(num == 7){
	  if(value == -1000)rhof0phase_ = rhof0phase_ + factor_ampl(Energy,15);	 
	  else rhof0phase_ = value;
	}
	if(num == 8){
	  if(value == -1000)rhoprhom_  = rhoprhom_ + factor_ampl(Energy,17);	
	  else rhoprhom_ = value;
	}
	if(num == 9){
	  if(value == -1000)rhoprhom_phase = rhoprhom_phase +  factor_ampl(Energy,19);
	  else rhoprhom_phase = value;
	}
	if(num == 10){
	  if(value == -1000)rhosigma_ = rhosigma_ + factor_ampl(Energy,21);
	  else rhosigma_ = value;
	}
	if(num == 11){
	  if(value == -1000)rhosigmaphase_ =rhosigmaphase_ + factor_ampl(Energy,23);
	  else 	  rhosigmaphase_ = value;
	}
	if(num == 12){
	  if(value == -1000)phsp_ =   phsp_ +  factor_ampl(Energy,25);
	  else 	  phsp_ = value;

	}
	if(num == 13){
	  if(value == -1000)phspphase_ =phspphase_ + factor_ampl(Energy,27);
	  else 	  phspphase_ = value;
	}
	if(num == 14){
	  if(value == -1000)a2pi_ =  a2pi_ +  factor_ampl(Energy,29);
	  else 	  a2pi_ = value;
	}
	if(num == 15){
	  if(value == -1000)a2pi_ph =a2pi_ph + factor_ampl(Energy,31);
	  else a2pi_ph = value;
	}
	if(num == 16){
	  if(value == -1000)h1pi =  h1pi +  factor_ampl(Energy,33);
	  else 	  h1pi = value;
	}
	if(num == 17){
	  if(value == -1000)h1pi_ph = h1pi_ph+  factor_ampl(Energy,35);
	  else 	  h1pi_ph = value;
	}
	if(num == 18){
	  if(value == -1000)a1pi_rhopi_II = a1pi_rhopi_II +  factor_ampl(Energy,37);
	  else 	  a1pi_rhopi_II = value;
	}
	if(num == 19){
	  if(value == -1000)a1pi_rhopi_phase_II = a1pi_rhopi_phase_II + factor_ampl(Energy,39);
	  else a1pi_rhopi_phase_II = value;
	}
	if(num == 20){
	  if(value == -1000)rhof0_inter = rhof0_inter + factor_ampl(Energy,41);
	  else 	  rhof0_inter = value;
	}
	if(num == 21){	
	  if(value == -1000)rhof0phase_inter = rhof0phase_inter + factor_ampl(Energy,43);
	  else 	  rhof0phase_inter = value;
	}
	if(num == 22){	
	  if(value == -1000)rhof2 = rhof2 + factor_ampl(Energy,45);
	  else 	  rhof2 = value;
	}
	if(num == 23){
	  if(value == -1000)rhof2_phase = rhof2_phase +  factor_ampl(Energy,47);
	  else rhof2_phase =value;
	}
	if(num == 24){	
	  if(value == -1000)rhoprhom_II = rhoprhom_II +  factor_ampl(Energy,49);
	  else 	  rhoprhom_II = value;
	}
	if(num == 25){
	  if(value == -1000)rhoprhom_phase_II = rhoprhom_phase_II + factor_ampl(Energy,51);
	  else 	  rhoprhom_phase_II = value;
	}
}


void Cmd3Generator2pi2pi0_ke::fill_parameters(){

    isISR = 0;
    omegapi0_ = factor_ompi0();
    omegapi0phase_ = 0.;
    
    a1pi_rhopi_ = factor_a1pi_rhopi();
    a1pi_rhopi_phase = factor_a1pi_rhopi_ph();

    a1pi_sigmapi_ = factor_a1pi_sigmapi();
    a1pi_sigmapi_phase = factor_a1pi_sigmapi_ph();

    rhof0_ = factor_rho_f0();
    rhof0phase_ = factor_rho_f0_ph();

    rhosigma_ = factor_rho_sigma();
    rhosigmaphase_ = factor_rho_sigma_ph();

    phsp_ = factor_phsp();
    phspphase_ = factor_phsp_ph();

    rhoprhom_ = factor_rhop_rhom();
    rhoprhom_phase = factor_rhop_rhom_ph();

    h1pi = factor_h1_rhopi();
    h1pi_ph = factor_h1_rhopi_ph();
    
    a2pi_ = factor_a2pi();
    a2pi_ph = factor_a2pi_ph();
    a1pi_rhopi_II =  factor_ampl(Energy,36);
    a1pi_rhopi_phase_II = factor_ampl(Energy,38);

    rhof0phase_inter = factor_ampl(Energy,42);
    rhof0_inter =  factor_ampl(Energy,40);

    rhof2 = factor_ampl(Energy,44);
    rhof2_phase = factor_ampl(Energy,46);
    rhoprhom_II = factor_ampl(Energy,48);
    rhoprhom_phase_II = factor_ampl(Energy,50);
    
    double total = omegapi0_+a1pi_rhopi_+a1pi_sigmapi_+rhof0_+rhosigma_+phsp_;
    total = total + rhoprhom_ + h1pi + a2pi_ + a1pi_rhopi_II;
    
    fompi0 = omegapi0_/total;dompi0 = factor_ampl(Energy,3);
    fa1pi_rhopi = a1pi_rhopi_/total; da1pi_rhopi = factor_ampl(Energy,5);
    fa1pi_sigmapi = a1pi_sigmapi_/total; da1pi_sigmapi = factor_ampl(Energy,9);
    frhoprhom = rhoprhom_/total; drhoprhom = factor_ampl(Energy,17);
    frhof0 = rhof0_/total; drhof0 = factor_ampl(Energy,13);
    frhosigma = rhosigma_/total; drhosigma = factor_ampl(Energy,21);
    fh1pi=h1pi/total; dh1pi= factor_ampl(Energy,33);
    
    
    
}


void Cmd3Generator2pi2pi0_ke::fill0_parameters(){
    isISR = 0;
    omegapi0_ = 0;a1pi_rhopi_ = 0.;a1pi_sigmapi_ = 0.;rhof0_ = 0.;rhosigma_ = 0.;phsp_ = 0.;
    rhoprhom_ = 0;h1pi = 0.;a2pi_ = 0.;a1pi_rhopi_II = 0.;
    rhoprhom_II = 0.;rhof2 = 0.;
}
void Cmd3Generator2pi2pi0_ke::fill_omega_pi0(){
    omegapi0_ = factor_ompi0();
}
void Cmd3Generator2pi2pi0_ke::fill_a1pi_rhopi(){
    a1pi_rhopi_ = factor_a1pi_rhopi();
}
void Cmd3Generator2pi2pi0_ke::fill_a1pi_sigmapi(){
    a1pi_sigmapi_ = factor_a1pi_sigmapi();
}
void Cmd3Generator2pi2pi0_ke::fill_rhof0(){
    rhof0_ = factor_rho_f0();
}
void Cmd3Generator2pi2pi0_ke::fill_rhosigma(){
    rhosigma_ = factor_rho_sigma();
}
void Cmd3Generator2pi2pi0_ke::fill_phsp(){
   phsp_ = factor_phsp();
}
void Cmd3Generator2pi2pi0_ke::fill_rhoprhom(){
    rhoprhom_ = factor_rhop_rhom();
}
void Cmd3Generator2pi2pi0_ke::fill_h1pi(){
    h1pi = factor_h1_rhopi();
}
void Cmd3Generator2pi2pi0_ke::fill_a2pi(){
    a2pi_ = factor_a2pi();
}



void Cmd3Generator2pi2pi0_ke::print(){
	cout << "============================================================" << endl;
	cout << "Energy = " << Energy << endl;
	cout << "isISR = " << isISR << endl; 

	cout << "omegapi0_ = " << omegapi0_ << endl;

	cout << "a1pi_rhopi_ = " << a1pi_rhopi_ << endl;
	cout << "a1pi_rhopi_phase = " << a1pi_rhopi_phase << endl;

	cout << "a1pi_sigmapi_ = " << a1pi_sigmapi_ << endl;
	cout << "a1pi_sigmapi_phase = " << a1pi_sigmapi_phase << endl;

	cout << "rhoprhom_ = " << rhoprhom_ << endl;
	cout << "rhoprhom_phase = " << rhoprhom_phase << endl;
	cout << "rhoprhom_II  = " << rhoprhom_II << endl;
        cout << "rhoprhom_phase_II = " << rhoprhom_phase_II << endl;

	cout << "rhof0_ = " << rhof0_ << endl;
	cout << "rhof0phase_ = " <<rhof0phase_ << endl;

	cout << "rhof0_inter = " << rhof0_inter << endl;
	cout << "rhof0phase_inter = " <<rhof0phase_inter << endl;

	cout << "rhosigma_ = " << rhosigma_ << endl;
	cout << "rhosigmaphase_ = " << rhosigmaphase_ << endl;

	cout << "rhof2 = " << rhof2 << endl;
	cout << "rhof2_phase = " << rhof2_phase << endl;

	cout << "h1pi = " << h1pi << endl;
	cout << "h1pi_ph = " << h1pi_ph << endl;

        cout << "a1pi_rhopi_II = " << a1pi_rhopi_II << endl;
        cout << "a1pi_rhopi_phase_II = " << a1pi_rhopi_phase_II << endl;

	cout << "phsp_ = " << phsp_ << endl;
	cout << "phspphase_ = " << phspphase_ << endl;

	cout << "a2pi_ = " << a2pi_ << endl;
	cout << "a2pi_ph = " << a2pi_ph << endl;
        cout << "============================================================" << endl;
}

double Cmd3Generator2pi2pi0_ke::levicivita(int A, int B, int C, int D){
	if(A == 1 && B == 2 && C == 3 && D == 4)return 1;
	else if(A == 1 && B == 2 && C == 4 && D == 3)return -1;
	else if(A == 1 && B == 3 && C == 2 && D == 4)return 1;
	else if(A == 1 && B == 3 && C == 4 && D == 2)return -1;
	else if(A == 1 && B == 4 && C == 2 && D == 3)return 1;
	else if(A == 1 && B == 4 && C == 3 && D == 2)return -1;
	else if(A == 2 && B == 3 && C == 4 && D == 1)return 1;
	else if(A == 2 && B == 4 && C == 3 && D == 1)return -1;
	else if(A == 3 && B == 2 && C == 4 && D == 1)return 1;
	else if(A == 3 && B == 4 && C == 2 && D == 1)return -1;
	else if(A == 4 && B == 2 && C == 3 && D == 1)return 1;
	else if(A == 4 && B == 3 && C == 2 && D == 1)return -1;
	else if(A == 3 && B == 4 && C == 1 && D == 2)return 1;
	else if(A == 4 && B == 3 && C == 1 && D == 2)return -1;
	else if(A == 2 && B == 4 && C == 1 && D == 3)return 1;
	else if(A == 4 && B == 2 && C == 1 && D == 3)return -1;
	else if(A == 2 && B == 3 && C == 1 && D == 4)return 1;
	else if(A == 3 && B == 2 && C == 1 && D == 4)return -1;
	else if(A == 4 && B == 1 && C == 2 && D == 3)return 1;
	else if(A == 3 && B == 1 && C == 2 && D == 4)return -1;
	else if(A == 4 && B == 1 && C == 3 && D == 2)return 1;
	else if(A == 2 && B == 1 && C == 3 && D == 4)return -1;
	else if(A == 3 && B == 1 && C == 4 && D == 2)return 1;
	else if(A == 2 && B == 1 && C == 4 && D == 3)return -1;
	else return 0.;
}

double Cmd3Generator2pi2pi0_ke::TLorentzVectorC_index(TLorentzVectorC P, int index){
	
	if(index == 0)return P.E().real();
	else if(index == 1)return P.X().real();
	else if(index == 2)return P.Y().real();
	else if(index == 3)return P.Z().real();
	else return 0.;
	
}

double Cmd3Generator2pi2pi0_ke::svertka_leviciv(TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
	double res = 0.;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			for(int k = 0; k < 4; k++){
				for(int l = 0; l < 4; l++){
				
					res = res + levicivita(i,j,k,l)*TLorentzVectorC_index(P1,i)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
					
				}
			}
		}
	}
	return res;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::sver_vec_leviciv(TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
	double res0 = 0.;
	double res1 = 0.;
	double res2 = 0.;
	double res3 = 0.;
	for(int j = 0; j < 4; j++){
		for(int k = 0; k < 4; k++){
			for(int l = 0; l < 4; l++){
				
				res0 = res0 + levicivita(1,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
				res1 = res1 + levicivita(2,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
				res2 = res2 + levicivita(3,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);
				res3 = res3 + levicivita(4,j,k,l)*TLorentzVectorC_index(P2,j)*TLorentzVectorC_index(P3,k)*TLorentzVectorC_index(P4,l);

			}
		}
	}
	
	return TLorentzVectorC(res0, res1, res2, res3);
}



// ------------------------------- PHASE SPACE------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_ph_sp(double *e, double *par){
    double norm[] = {0.027171, 0.0258879, 0.0247892, 0.0238614, 0.0230847, 0.0223951, 0.021841, 0.0211563, 0.0207638, 0.0203469, 0.0200074, 0.0196746, 0.0192915, 0.0190267, 0.0187609, 0.0185032, 0.0182173, 0.0180387, 0.0177281, 0.0174934, 0.017325, 0.0172337, 0.0170029, 0.0168753, 0.0167488, 0.0165608, 0.0164036, 0.0162988};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_ph_sp_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC result = (P3-P4);//(0,1./sqrt(2),1./sqrt(2),0);
    //TLorentzVectorC result = result0*(P1+P2).Dot(P3);
    //result = result - result0*(P1+P2).Dot(P4);
    double ener = P.E().real();
    return result/sqrt(cross_ph_sp(&ener, &ener));


}



TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_ph_sp(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    return H_ph_sp_beforepi0sim(P,P1,P2,P3,P4) + H_ph_sp_beforepi0sim(P,P2,P1,P3,P4) - H_ph_sp_beforepi0sim(P,P1,P2,P4,P3) - H_ph_sp_beforepi0sim(P,P2,P1,P4,P3);

}


double Cmd3Generator2pi2pi0_ke::factor_phsp(){
    return factor_ampl(Energy,24);
}

double Cmd3Generator2pi2pi0_ke::factor_phsp_ph(){
    return factor_ampl(Energy,26);
}


// ------------------------- OMEGA PI0 -----------------------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_omega_pi0(double *e, double *par){
    double norm[] = {1.87049e-07, 1.56049e-06, 7.99364e-06, 3.37705e-05, 0.000135718, 0.000596526, 0.0107153, 0.100773, 0.216074, 0.309316, 0.401808, 0.4623,0.538993, 0.560054, 0.602989, 0.667697, 0.667369, 0.785723, 0.806614, 0.835757, 0.854405, 0.862221, 0.87398, 0.800103, 0.758419, 0.752406, 0.783933, 0.925019};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_omega_pi0_init(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){


    TLorentzVectorC W = P2 + P3 + P4;
    TLorentzVectorC k = P3 + P4;
    TLorentzVectorC result = (P3-P4)*P.Dot(k)*W.Dot(W) + W*P.Dot(P3-P4)*W.Dot(k) + k*P.Dot(W)*W.Dot(P3-P4);
    result = result - W*P.Dot(k)*W.Dot(P3-P4) - k*P.Dot(P3-P4)*W.Dot(W) - (P3-P4)*P.Dot(W)*W.Dot(k);
    result = result*propagator_omega(W.M()*W.M())*propagator_rho(k.M()*k.M());
    return result;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_omega_pi0_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    return - H_omega_pi0_init(P,P1,P2,P3,P4) + H_omega_pi0_init(P,P1,P4,P3,P2) + H_omega_pi0_init(P,P1,P3,P2,P4);
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_omega_pi0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC a = H_omega_pi0_beforepi0sim(P,P1,P2,P3,P4) + H_omega_pi0_beforepi0sim(P,P2,P1,P3,P4);
    double ener = P.E().real();
    return  a/sqrt(cross_omega_pi0(&ener, &ener));
}


double Cmd3Generator2pi2pi0_ke::factor_ompi0(){
    return factor_ampl(Energy,2);
}
double Cmd3Generator2pi2pi0_ke::factor_ompi0_ph(){
    return -1.;
    return factor_ampl(Energy,2);
}

//============================================================
// ------------A1 (rho pi) PI first part----------------------
//============================================================
double Cmd3Generator2pi2pi0_ke::cross_a1_rhopi_pi(double *e, double *par){
    double norm[] = {9.38941e-06, 3.37591e-05, 0.000104477, 0.000286863, 0.000767181, 0.00200169, 0.00523157, 0.0145385, 0.0411876, 0.114388, 0.2936, 0.639556, 1.20124, 2.06013, 3.09992, 4.44872, 6.03592, 7.64013, 9.14649, 10.8322, 11.9934, 12.988, 14.3327, 14.3168, 15.4297, 16.1528, 16.5091, 17.1639};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P1 + P2 + P3;
    TLorentzVectorC k = P2 + P3;
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + (P3-P2)*P.Dot(W)*W.Dot(k);
    result = result - W*P.Dot(P3-P2)*W.Dot(k);
    result = result + W*P.Dot(k)*W.Dot(P3-P2);
    result = result - k*W.Dot(P3-P2)*P.Dot(W);
    double ener = P.E().real();
    double ph = 0.;//1.*(W.M()*W.M()/Ma1/Ma1-1.);
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    return -result*propagator_rho(k.M()*k.M())*propagator_a1(W.M()*W.M())/sqrt(cross_a1_rhopi_pi(&ener, &ener))*FF; 
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_a1_rhopi_pi_beforepi0sim(P,P1,P2,P3,P4);
    result = result + H_a1_rhopi_pi_beforepi0sim(P,P2,P1,P3,P4);
    result = result - H_a1_rhopi_pi_beforepi0sim(P,P1,P2,P4,P3);
    result = result - H_a1_rhopi_pi_beforepi0sim(P,P2,P1,P4,P3);

    return result;

}

double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi(){
    return factor_ampl(Energy,4);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_ph(){
    return factor_ampl(Energy,6);
}

//=============================================
// ------A1 (rho pi) PI second part -----------
//=============================================

double Cmd3Generator2pi2pi0_ke::cross_a1_rhopi_pi_II(double *e, double *par){
    double norm[] = {1.89775e-06, 6.93526e-06, 2.02051e-05, 5.56619e-05, 0.000144712, 0.000371023, 0.000920294, 0.00242596, 0.00627131, 0.0179908, 0.0418408, 0.0831727, 0.159426, 0.26784, 0.408592, 0.63246, 0.920412, 1.38438, 2.03832, 3.08834, 4.68597, 7.24819, 12.4944, 19.9858, 32.6581, 45.0234, 56.3054, 64.2716};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_beforepi0sim_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    double ener = P.E().real();

/*
    TLorentzVectorC W = P3 + P4;
    TLorentzVectorC result(0.,0.,0.,0.);
    result = W*P.Dot(P3-P4);
    double ph = 0.;//1.7*(W.M()*W.M()/Mf0/Mf0 - 1.);
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    result = result*propagator_rho(W.Dot(W))*propagator_f0((P1+P2).Dot(P1+P2))*FF;
    return -result/sqrt(cross_rho_f0(&ener, &ener));
    
*/	
    TLorentzVectorC W = P1 + P2;
    TLorentzVectorC e = P3 - P4;
    W = P1 + P2 + P4;
    TLorentzVectorC k = P1 + P4;
    e = P4 - P1;
    TLorentzVectorC result = (W*P.Dot(P3) - P3*W.Dot(P))*(e.Dot(P2)*k.Dot(W)-e.Dot(W)*k.Dot(P2));
    result = result*propagator_pi1300(W.Dot(W))*propagator_rho(k.Dot(k));
 
    return result/sqrt(cross_a1_rhopi_pi_II(&ener, &ener));

}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_a1_rhopi_pi_beforepi0sim_II(P,P1,P2,P3,P4);
    //result = result + H_a1_rhopi_pi_beforepi0sim_II(P,P2,P1,P3,P4);
    //result = result - H_a1_rhopi_pi_beforepi0sim_II(P,P1,P2,P4,P3);
    //result = result - H_a1_rhopi_pi_beforepi0sim_II(P,P2,P1,P4,P3);
    return result;
}

double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_II(){
    return factor_ampl(Energy,36);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_rhopi_ph_II(){
    return factor_ampl(Energy,38);
}


//==================================================
// ---------A1 (sigma pi) PI -----------------------
//==================================================

double Cmd3Generator2pi2pi0_ke::cross_a1_sigmapi_pi(double *e, double *par){
    double norm[] = {2.35448e-05, 7.50863e-05, 0.000206688, 0.000509195, 0.00120981, 0.00277813, 0.00631554, 0.0142642, 0.0309426, 0.0652028, 0.135337, 0.268977, 0.499254, 0.88869, 1.45045, 2.3462, 3.49934, 4.77602, 6.22402, 7.7155, 9.27811, 10.428, 11.9285, 13.1289, 14.2201, 15.8061, 16.4864, 18.3851};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_sigmapi_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC W = P1 + P2 + P3;
	TLorentzVectorC k = P1 + P2;
	double ener = P.E().real();
	TLorentzVectorC result(0.,0.,0.,0.);
        result = result + P3*P.Dot(W)*k.Dot(W);
	result = result - k*P.Dot(W)*P3.Dot(W);
	result = result - W*P3.Dot(P)*k.Dot(W);
	result = result + W*P3.Dot(W)*k.Dot(P);
        double ph = 0.;//3.*(W.M()*W.M()/Ma1/Ma1-1.);
        complex<double>  FF = complex<double> (cos(ph),sin(ph));
	result = result*propagator_a1(W.M()*W.M())*propagator_sigma(k.M()*k.M())*FF;
	return -result/sqrt(cross_a1_sigmapi_pi(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a1_sigmapi_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = res + H_a1_sigmapi_pi_beforepi0sim(P,P1,P2,P3,P4) - H_a1_sigmapi_pi_beforepi0sim(P,P1,P2,P4,P3);
    res = res + H_a1_sigmapi_pi_beforepi0sim(P,P2,P1,P3,P4) - H_a1_sigmapi_pi_beforepi0sim(P,P2,P1,P4,P3);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_a1pi_sigmapi(){
    return factor_ampl(Energy,8);
}
double Cmd3Generator2pi2pi0_ke::factor_a1pi_sigmapi_ph(){
    return factor_ampl(Energy,10);
}

double ratio_cross(double *e, double *par){

    return Part_ana.cross_a1_rhopi_pi(e,par)/Part_ana.cross_a1_sigmapi_pi(e,par);

}

//===============================================
//------------------rho f2-----------------------
//===============================================

double Cmd3Generator2pi2pi0_ke::cross_rho_f2(double *e, double *par){
    double norm[] = {0.00203014, 0.00305518, 0.00340913, 0.00319095, 0.00271788, 0.00212643, 0.00159291, 0.00111007, 0.000769828, 0.000456646, 0.000340798, 0.000284551, 0.000250472, 0.000240615, 0.00024844, 0.000257863, 0.000258196, 0.000253757, 0.000262872, 0.000250844, 0.000238593, 0.000216461, 0.000181339, 0.000164848, 0.000141188, 0.000147802, 0.000161769, 0.000176241};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f2(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    double ener = P.E().real();
    TLorentzVectorC W = P1 + P2;
    TLorentzVectorC e = P3 - P4;

    double nomin = W.Dot(W).real();
    TLorentzVectorC result = (P2 - W*W.Dot(P2)/nomin)*(e.Dot(P1) - e.Dot(W)*P1.Dot(W)/nomin)/2.;
    result = result + (P1 - W*W.Dot(P1)/nomin)*(e.Dot(P2) - e.Dot(W)*P2.Dot(W)/nomin)/2.;
    result = result - (e-W*W.Dot(e)/nomin)*(P1.Dot(P2) - W.Dot(P1)*W.Dot(P2)/nomin)/3.;
    result = result - P*result.Dot(P)/P.Dot(P).real();
    result = result*propagator_rho((P3+P4).Dot(P3+P4))*propagator_f2(W.Dot(W));   

    return result/sqrt(cross_rho_f2(&ener, &ener));

}





//-------------------------------- A2 PI------------------------------------------
double Cmd3Generator2pi2pi0_ke::cross_a2_pi(double *e, double *par){
    double norm[] = {5.09254e-14, 1.12193e-12, 1.03475e-11, 6.11362e-11, 2.81282e-10, 1.07244e-09, 3.61526e-09, 1.11212e-08, 3.23475e-08, 8.83757e-08, 2.39315e-07, 6.36937e-07, 1.67296e-06, 4.37807e-06, 1.14225e-05, 2.95183e-05, 7.61204e-05, 0.000201503, 0.00052664, 0.00121244, 0.00256362, 0.0049281, 0.00854306, 0.0138288, 0.0221141, 0.0323682, 0.0459657, 0.0659924};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a2_pi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC W = P1 + P2 + P3;
    double ener = P.E().real();
    TLorentzVectorC result = - sver_vec_leviciv(W,P4,P1)*svertka_leviciv(W, P3-P2, P1, P4);
    TLorentzVectorC btv = (P3-P2)*W.Dot(W)*P4.Dot(P1) + P1*P4.Dot(W)*(P3-P2).Dot(W);
    btv = btv + W*W.Dot(P1)*P4.Dot(P3-P3) - P1*W.Dot(W)*P4.Dot(P3-P3);
    btv = btv - W*P4.Dot(P1)*(P3-P2).Dot(W) - (P3-P2)*W.Dot(P1)*W.Dot(P4);
    result = result + btv*(P4.Dot(P1) - W.Dot(P4)*W.Dot(P1)/W.Dot(W));
    result = result*propagator_a2(W.Dot(W))*propagator_rho((P3+P2).Dot(P3+P2));
    return result/sqrt(cross_a2_pi(&ener, &ener));

}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_a2_pi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC res = H_a2_pi_beforepi0sim(P,P1,P2,P3,P4)  + H_a2_pi_beforepi0sim(P,P2,P1,P3,P4);
    return res;
    
}

double Cmd3Generator2pi2pi0_ke::factor_a2pi(){
    return factor_ampl(Energy,28);
}
double Cmd3Generator2pi2pi0_ke::factor_a2pi_ph(){
    return factor_ampl(Energy,30);
}

//-----------------------------------------------
// ------------rho sigma-------------------------
//-----------------------------------------------

double Cmd3Generator2pi2pi0_ke::cross_rho_sigma(double *e, double *par){

    double norm[] = {0.000876075, 0.00203927, 0.00413222, 0.00768718, 0.0136843, 0.0242913, 0.0427909, 0.0782594, 0.150395, 0.282294, 0.520047, 0.805074, 1.22889, 1.65961, 2.33903, 2.99221, 3.97308, 4.4554, 5.1025, 5.70802, 6.67461, 6.67031, 7.22386, 7.26811, 7.27732, 7.9501, 7.80068, 8.19164};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_sigma_beforesim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    double ph = 1.5;//2.98311 -2.57009*ener + 0.378539*ener*ener;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    TLorentzVectorC result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    result = result*propagator_rho(W.Dot(W))*FF*propagator_sigma((P1+P2).Dot(P1+P2));
    return -result/sqrt(cross_rho_sigma(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_sigma(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = res + H_rho_sigma_beforesim(P,P1,P2,P3,P4);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_rho_sigma(){
    return factor_ampl(Energy,20);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_sigma_ph(){
    return factor_ampl(Energy,22);
}



//---------------------------------------------------------------------------
// -----------------------------------RHO F0---------------------------------
//---------------------------------------------------------------------------
double Cmd3Generator2pi2pi0_ke::cross_rho_f0(double *e, double *par){

    double norm[] = {0.000120128, 0.0002924, 0.000618387, 0.00117923, 0.00212543, 0.00378536, 0.0067353, 0.0122104, 0.0233578, 0.0438315, 0.0820602, 0.127734, 0.194391, 0.260065, 0.361597, 0.46138, 0.626127, 0.746501, 0.949543, 1.23827, 1.72136, 2.1918, 3.04271, 3.61597, 4.39168, 4.85355, 5.56185, 5.57303};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    double ph = 0.;//1.3*((P1+P2).M()-0.98)/0.08;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    result = result*propagator_rho(W.Dot(W))*propagator_f0((P1+P2).Dot(P1+P2))*FF;
    return result/sqrt(cross_rho_f0(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = H_rho_f0_before_sim(P,P1,P2,P3,P4);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0(){
    return factor_ampl(Energy,12);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_ph(){
    return factor_ampl(Energy,14);
}

//--------------------------------------------------------------------------
// -----------------------------------RHO f0_inter--------------------------
//--------------------------------------------------------------------------
double Cmd3Generator2pi2pi0_ke::cross_rho_f0_inter(double *e, double *par){

    double norm[] = {0.000120128, 0.0002924, 0.000618387, 0.00117923, 0.00212543, 0.00378536, 0.0067353, 0.0122104, 0.0233578, 0.0438315, 0.0820602, 0.127734, 0.194391, 0.260065, 0.361597, 0.46138, 0.626127, 0.746501, 0.949543, 1.23827, 1.72136, 2.1918, 3.04271, 3.61597, 4.39168, 4.85355, 5.56185, 5.57303};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_inter_before_sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P3 + P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P3-P4)*P.Dot(W) - W*P.Dot(P3-P4);
    double ph = 0.;//-2.33713 + 6.89022*ener -3.23062*ener*ener;
    complex<double>  FF = complex<double> (cos(ph),sin(ph));
    result = result*propagator_rho(W.Dot(W))*propagator_f0_inter((P1+P2).Dot(P1+P2))*FF;
    return -result/sqrt(cross_rho_f0_inter(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rho_f0_inter(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC res(0.,0.,0.,0.);
    res = H_rho_f0_inter_before_sim(P,P1,P2,P3,P4);
    return res;
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_inter(){
    return factor_ampl(Energy,12);
}

double Cmd3Generator2pi2pi0_ke::factor_rho_f0_inter_ph(){
    return factor_ampl(Energy,14);
}


//==================================================
//---------------------rho+ rho- -------------------
//==================================================

double Cmd3Generator2pi2pi0_ke::cross_rhop_rhom(double *e, double *par){

    double norm[] = {1.98223e-07, 1.43368e-06, 6.12057e-06, 2.02618e-05, 5.76425e-05, 0.000145236, 0.000351006, 0.000812736, 0.00177885, 0.00404068, 0.0088983, 0.0193122, 0.0398636, 0.0832108, 0.162197, 0.304643, 0.569157, 1.01843, 1.75296, 2.84571, 4.69834, 6.51253, 9.54061, 12.1267, 15.4348, 19.176, 22.3248, 27.298};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}


TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom_beforpi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC rhop = P1 + P3;
    TLorentzVectorC rhom = P2 + P4;
    TLorentzVectorC omega = P1 - P3;
    TLorentzVectorC eee = P2 - P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P2-P4)*(P1+P3-P2-P4).Dot(P1-P3) + (P1-P3)*(P1+P3-P2-P4).Dot(P2-P4);
    double ph1 = 0.;double ph2 = 0.;
    complex<double>  FF = complex<double> (cos(ph1),sin(ph1))*complex<double> (cos(ph2),sin(ph2));
    result = -result*propagator_rho(rhop.Dot(rhop))*propagator_rho(rhom.Dot(rhom))*FF;
    
    return result/sqrt(cross_rhop_rhom(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC result = H_rhop_rhom_beforpi0sim(P,P1,P2,P3,P4) + H_rhop_rhom_beforpi0sim(P,P2,P1,P3,P4);
    result = result - H_rhop_rhom_beforpi0sim(P,P1,P2,P4,P3) - H_rhop_rhom_beforpi0sim(P,P2,P1,P4,P3);
    return result;
}
double Cmd3Generator2pi2pi0_ke::factor_rhop_rhom(){
    return factor_ampl(Energy,16);
}
double Cmd3Generator2pi2pi0_ke::factor_rhop_rhom_ph(){
    return factor_ampl(Energy,18);
}


//==========================================
//---------------rho+rho-_II ---------------
//==========================================

double Cmd3Generator2pi2pi0_ke::cross_rhop_rhom_II(double *e, double *par){
    double norm[] = {5.09254e-14, 1.12193e-12, 1.03475e-11, 6.11362e-11, 2.81282e-10, 1.07244e-09, 3.61526e-09, 1.11212e-08, 3.23475e-08, 8.83757e-08, 2.39315e-07, 6.36937e-07, 1.67296e-06, 4.37807e-06, 1.14225e-05, 2.95183e-05, 7.61204e-05, 0.000201503, 0.00052664, 0.00121244, 0.00256362, 0.0049281, 0.00854306, 0.0138288, 0.0221141, 0.0323682, 0.0459657, 0.0659924};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_rhop_rhom_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC rhop = P1 + P3;
    TLorentzVectorC rhom = P2 + P4;
    TLorentzVectorC omega = P1 - P3;
    TLorentzVectorC eee = P2 - P4;
    double ener = P.E().real();
    TLorentzVectorC result(0.,0.,0.,0.);
    result = (P1 + P3 - P2 - P4)*(P2-P4).Dot(P1-P3);
    double ph1 = 0.; double ph2 = 0.;
    complex<double>  FF = complex<double> (cos(ph1),sin(ph1))*complex<double> (cos(ph2),sin(ph2));
    result = result*propagator_rho(rhop.Dot(rhop))*propagator_rho(rhom.Dot(rhom));
    
    return result/sqrt(cross_rhop_rhom_II(&ener, &ener));

}






//==================================================================
// ------------------h1 (rho pi) pi --------------------------------
//==================================================================

double Cmd3Generator2pi2pi0_ke::cross_h1_rhopi(double *e, double *par){
    double norm[] = {1.2186e-05, 4.41931e-05, 0.00013166, 0.000352621, 0.000941573, 0.00232855, 0.00587275, 0.0156255, 0.0420362, 0.117464, 0.298101, 0.656967, 1.31859,  2.26348, 3.4833, 4.9666, 6.33923, 7.45706, 8.57597, 9.34138, 10.635, 10.9817, 11.8845, 11.787, 12.3091, 12.6336, 12.9628, 13.0517};
    double factor;
    if(e[0] <= 0.625)factor = norm[0];
    else if(e[0] >= 1.975)factor = norm[27];
    else {
        int index = (int)((e[0] - 0.625)/0.05);
        double e1 = 0.625 + index*0.05;
        double e2 = 0.625 + (index+1.)*0.05;
        factor = norm[index]*(e2-e[0])/(e2-e1) + norm[index+1]*(-e1+e[0])/(e2-e1);
    }
    return factor;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi_beforepi0sim(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
    TLorentzVectorC W = P2 + P3 + P4;
    TLorentzVectorC k = P3 + P4;
    TLorentzVectorC result = (P3-P4)*P.Dot(W)*W.Dot(k);
    result = result - W*P.Dot(P3-P4)*W.Dot(k);
    result = result + W*P.Dot(k)*W.Dot(P3-P4);
    result = result - k*W.Dot(P3-P4)*P.Dot(W);
    double ener = P.E().real();
    return result*propagator_rho(k.Dot(k))*propagator_h1(W.Dot(W))/sqrt(cross_h1_rhopi(&ener, &ener));
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi_bef(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_h1_rhopi_beforepi0sim(P,P1,P4,P3,P2);
    result = result - H_h1_rhopi_beforepi0sim(P,P1,P2,P3,P4);
    result = result + H_h1_rhopi_beforepi0sim(P,P1,P3,P2,P4);
    return result;
}

TLorentzVectorC  Cmd3Generator2pi2pi0_ke::H_h1_rhopi(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	
    TLorentzVectorC result(0.,0.,0.,0.);
    result = result + H_h1_rhopi_bef(P,P1,P2,P3,P4);
    result = result + H_h1_rhopi_bef(P,P2,P1,P3,P4);
    return result;
}

double Cmd3Generator2pi2pi0_ke::factor_h1_rhopi(){
    return factor_ampl(Energy,32);
}
double Cmd3Generator2pi2pi0_ke::factor_h1_rhopi_ph(){
    return factor_ampl(Energy,34);
}


// ----------- Cross section function -------------------------------------
double  Cmd3Generator2pi2pi0_ke::cross2pi2pi0(double E){

    if(E < 0.7) return 0.;
    double en[100],cr[100];
    en[0] = 0.7; cr[0] = 0.0005;
    en[1] = 0.8; cr[1] = 0.1;
    en[2] = 0.8400000; cr[2] = 0.2;
    en[3] = 0.8750000; cr[3] = 0.3;
    en[4] = 0.9000000; cr[4] = 0.4;
    en[5] =  0.9250000; cr[5] = 1.184512;
    en[6] =   0.9500001; cr[6] = 2.972132;
    en[7] =  0.9750000; cr[7] =  4.066643;
    en[8] =  1.000000; cr[8] =  5.724451;
    en[9] =  1.025000; cr[9] =  8.523726;
    en[10] =  1.050000; cr[10] =  8.702085;
    en[11] =  1.075000; cr[11] =  11.97795;
    en[12] =  1.100000; cr[12] =  13.26839;
    en[13] = 1.125000; cr[13] =  15.06046;
    en[14] = 1.150000; cr[14] =  16.71859;
    en[15] = 1.175000; cr[15] =  17.34131;
    en[16] = 1.200000; cr[16] =  19.80570;
    en[17] = 1.225000; cr[17] =  20.18546;
    en[18] = 1.250000; cr[18] =  21.54455;
    en[19] = 1.275000; cr[19] =  23.38247;
    en[20] = 1.300000; cr[20] =  25.77932;
    en[21] = 1.325000; cr[21] =  26.29928;
    en[22] = 1.350000; cr[22] =  26.91924;
    en[23] = 1.375000; cr[23] =  28.40883;
    en[24] = 1.400000; cr[24] =  30.18806;
    en[25] = 1.425000; cr[25] =  31.97551;
    en[26] = 1.450000; cr[26] =  32.10345;
    en[27] = 1.475000; cr[27] =  32.26585;
    en[28] = 1.500000; cr[28] =  31.69104;
    en[29] = 1.525000; cr[29] =  29.54573;
    en[30] = 1.550000; cr[30] =  29.83364;
    en[31] = 1.575000; cr[31] =  27.92464;
    en[32] = 1.600000; cr[32] =  27.31240;
    en[33] = 1.625000; cr[33] =  25.89993;
    en[34] = 1.650000; cr[34] =  24.90640;
    en[35] = 1.675000; cr[35] =  23.58823;
    en[36] = 1.700000; cr[36] =  22.66760;
    en[37] = 1.725000; cr[37] =  20.37718;
    en[38] = 1.750000; cr[38] =  19.49763;
    en[39] = 1.775000; cr[39] =  16.94212;
    en[40] = 1.800000; cr[40] =  15.12758;
    en[41] = 1.825000; cr[41] =  12.53257;
    en[42] = 1.850000; cr[42] =  11.96361;
    en[43] = 1.875000; cr[43] =  10.39149;
    en[44] = 1.900000; cr[44] =  9.308997;
    en[45] = 1.925000; cr[45] =  9.284678;
    en[46] = 1.950000; cr[46] =  8.576995;
    en[47] = 1.975000; cr[47] =  8.684048;
    en[48] = 2.000000; cr[48] =  9.254687;
    en[49] = 2.025000; cr[49] =  9.181911;
    if(E > en[49])return cr[49];
    int i = 0;
    while(E > en[i]){i++;}
    if(i < 2){return cr[i];}
    return (E-en[i-1])/(en[i]- en[i-1])*(cr[i]-cr[i-1]) + cr[i-1];

}

double  Cmd3Generator2pi2pi0_ke::cross4pi(double E){
    if(E < 0.6125) return 0.;
    double en[100],cr[100];
    en[0] = 0.6125; cr[0] = 0.02;
    en[1] = 0.6375; cr[1] = 0.04;
    en[2] = 0.6625; cr[2] = 0.02;
    en[3] = 0.6875; cr[3] = 0.01;
    en[4] = 0.7125; cr[4] = 0.02;
    en[5] = 0.7375; cr[5] = 0.03;
    en[6] = 0.7625; cr[6] = 0.05;
    en[7] = 0.7875; cr[7] = 0.11;
    en[8] = 0.8125; cr[8] = 0.11;
    en[9] = 0.8375; cr[9] = 0.12;
    en[10] = 0.8625; cr[10] = 0.17;
    en[11] = 0.8875; cr[11] = 0.26;
    en[12] = 0.9125; cr[12] = 0.33;
    en[13] = 0.9375; cr[13] = 0.57;
    en[14] = 0.9625; cr[14] = 0.71;
    en[15] = 0.9875; cr[15] = 0.89;
    en[16] = 1.0125; cr[16] = 1.2;
    en[17] = 1.0375; cr[17] = 1.61;
    en[18] = 1.0625; cr[18] = 2.17;
    en[19] = 1.0875; cr[19] = 3.29;
    en[20] = 1.1125; cr[20] = 4.49;
    en[21] = 1.1375; cr[21] = 5.95;
    en[22] = 1.1625; cr[22] = 7.37;
    en[23] = 1.1875; cr[23] = 8.84;
    en[24] = 1.2125; cr[24] = 10.79;
    en[25] = 1.2375; cr[25] = 12.62;
    en[26] = 1.2625; cr[26] = 14.56;
    en[27] = 1.2875; cr[27] = 16.39;
    en[28] = 1.3125; cr[28] = 19.06;
    en[29] = 1.3375; cr[29] = 21.14;
    en[30] = 1.3625; cr[30] = 23.37;
    en[31] = 1.3875; cr[31] = 25.76;
    en[32] = 1.4125; cr[32] = 27.53;
    en[33] = 1.4375; cr[33] = 29.95;
    en[34] = 1.4625; cr[34] = 30.32;
    en[35] = 1.4875; cr[35] = 32.04;
    en[36] = 1.5125; cr[36] = 30.98;
    en[37] = 1.5375; cr[37] = 30.11;
    en[38] = 1.5625; cr[38] = 28.26;
    en[39] = 1.5875; cr[39] = 26.81;
    en[40] = 1.6125; cr[40] = 24.66;
    en[41] = 1.6375; cr[41] = 22.69;
    en[42] = 1.6625; cr[42] = 20.95;
    en[43] = 1.6875; cr[43] = 18.78;
    en[44] = 1.7125; cr[44] = 17.25;
    en[45] = 1.7375; cr[45] = 15.33;
    en[46] = 1.7625; cr[46] = 13.37;
    en[47] = 1.7875; cr[47] = 11.61;
    en[48] = 1.8125; cr[48] = 10.23;
    en[49] = 1.8375; cr[49] = 8.87;
    en[50] = 1.8625; cr[50] = 7.67;
    en[51] = 1.8875; cr[51] = 7.29;
    en[52] = 1.9125; cr[52] = 7.17;
    en[53] = 1.9375; cr[53] = 6.93;
    en[54] = 1.9625; cr[54] = 6.54;
    en[55] = 1.9875; cr[55] = 6.04;
    en[56] = 2.0125; cr[56] = 6.18;
    en[57] = 2.0375; cr[57] = 5.66;
    en[58] = 2.0625; cr[58] = 5.68;
    en[59] = 2.0875; cr[59] = 5.34;
    en[60] = 2.1125; cr[60] = 4.92;
    en[61] = 2.1375; cr[61] = 4.83;
    en[62] = 2.1625; cr[62] = 4.59;
    en[63] = 2.1875; cr[63] = 4.28;
    if(E > en[63])return cr[63];
    int i = 0;
    while(E > en[i]){i++;}
    if(i < 2){return cr[i];}
    return (E-en[i-1])/(en[i]- en[i-1])*(cr[i]-cr[i-1]) + cr[i-1];
}




//------------------ Generate events ------------------------------------------
double radiatorW(double s, double x){
    double L = 2.*TMath::Log(sqrt(s)/me);
    //return 2.*alpha/3.141592*(L-1.)*(1.-x+x*x/2.)/x;
    double dzeta3 = 1.2020569;
    double dzeta2 = 1.64493407;
    double beta = 2.*alpha1/3.141592*(L-1.);
    double delta2 = (9./8. - 2*dzeta2)*L*L - (45./16.-11./2.*dzeta2-3.*dzeta3)*L-6./5.*dzeta2*dzeta2-9./2.*dzeta3-6.*dzeta2*TMath::Log(2.)+3./8.*dzeta2+57./12.;
    double delta = 1. + alpha1/3.141592*(1.5*L+1./3.*3.141592*3.141592-2.) + alpha1*alpha1/3.141592/3.141592*delta2;
    //cout << delta*beta*pow(x,beta-1.) << " " << beta/2.*(2.-x) << " " << beta*beta/8.*((2.-x)*(3.*TMath::Log(1.-x)-4.*TMath::Log(x))-4.*TMath::Log(1.-x)/x-6.+x) << endl;
    return delta*beta*pow(x,beta-1.)-beta/2.*(2.-x)+beta*beta/8.*((2.-x)*(3.*TMath::Log(1.-x)-4.*TMath::Log(x))-4.*TMath::Log(1.-x)/x-6.+x);
}

double radiatorW_th(double E, double x, double costh){
    double ss = E*E;
    double s = sqrt(1. - costh*costh);
    double c = costh;
    double one = s*s;
    double two = - x*x*s*s*s*s/2./(x*x - 2.*x + 2.);
    double three = - me*me/ss*4.*((1.-2.*x)*s*s-x*x*c*c*c*c)/(x*x - 2.*x + 2.);
    double P = (one + two + three)/pow(s*s+4.*me*me/ss*c*c,2.);
    return P;
}

double radiatorW_th_integral(double E, double x){
    double L = 2.*TMath::Log(E/me);
    return L - 1.;
}

double radiatorW_delta_x_beta(double s, double x){
    double L = 2.*TMath::Log(sqrt(s)/me);
    double dzeta3 = 1.2020569;
    double dzeta2 = 1.64493407;
    double beta = 2.*alpha1/3.141592*(L-1.);
    double delta2 = (9./8. - 2*dzeta2)*L*L - (45./16.-11./2.*dzeta2-3.*dzeta3)*L-6./5.*dzeta2*dzeta2-9./2.*dzeta3-6.*dzeta2*TMath::Log(2.)+3./8.*dzeta2+57./12.;
    double delta = 1. + alpha1/3.141592*(1.5*L+1./3.*3.141592*3.141592-2.) + alpha1*alpha1/3.141592/3.141592*delta2;
    return delta*pow(x,beta) - 0.5*beta*(2*x-0.5*x*x) + beta*beta/8.*0.0612457; // actually, it is valid just for x = 0.001
}
//-------------- Visible cross section ---------------------------------
double Cmd3Generator2pi2pi0_ke::cross(double E){
    //cout << E << " " << cross2pi2pi0(E) << endl;
    if(Mode2pi2pi0 == 1)return cross2pi2pi0(E);
    return cross4pi(E);
}

double Cmd3Generator2pi2pi0_ke::visible_cross_section(double E){
    double s = E*E;
    double result = 0.;
    double x = 0.001;
    double stepx = 0.00001;
    result = cross(E)*radiatorW_delta_x_beta(s,x);
    double xA = x;
    double xB = xA;
    while(sqrt(fabs(s*(1.-x - stepx))) > 4.*Mpi && (x+stepx) < 1.){
        xA = xB;
        xB = xA + stepx;
        double prirost = (cross(sqrt(fabs(s*(1.-xA))))*radiatorW(s,xA) + cross(sqrt(fabs(s*(1.-0.5*xA-0.5*xB))))*radiatorW(s,0.5*xA + 0.5*xB)*4. + cross(sqrt(fabs(s*(1.-xB))))*radiatorW(s,xB))/6.*stepx;
        stepx = stepx/prirost*result*0.001;
    if(sqrt(fabs(s*(1.-x - stepx))) < 4.*Mpi)break;
        result = result + prirost;
        x = x + stepx;
    }
    return result;
}

// -----------  generate x ---------------------------------------
double Cmd3Generator2pi2pi0_ke::generate_photon_energy(double E){
    double radiatorW_delta_x_beta(double, double);
    double radiatorW(double, double);
    double xnew = gRandom->Rndm();
    double s = E*E;
    double Eprime = E;
    double crvis_x = 0.;
    double crvis = visible_cross_section(E);
    double x = 0.001;
    double stepx = 0.00001;
    crvis_x = cross(E)*radiatorW_delta_x_beta(s,x);
    //cout << cross(E) << endl;
    //cout << "E = " << E <<  " cross(E) = " << cross(E) << "    radiatorW_delta_x_beta(s,x) = " << radiatorW_delta_x_beta(s,x) << "   visible_cross_section(E) = " << visible_cross_section(E) << endl;
    if(xnew < crvis_x/crvis){return 0.;}
    double xA = x;
    double xB = xA;
    while(xnew >= crvis_x/crvis){
        xA = xB;
        xB = xA + stepx;
        double prirost = (cross(sqrt(fabs(s*(1.-xA))))*radiatorW(s,xA) + cross(sqrt(fabs(s*(1.-0.5*xA-0.5*xB))))*radiatorW(s,0.5*xA + 0.5*xB)*4. + cross(sqrt(fabs(s*(1.-xB))))*radiatorW(s,xB))/6.*stepx;
        stepx = stepx/prirost*crvis_x*0.001;
        crvis_x = crvis_x + prirost;
        x = x + stepx;
        if(sqrt(fabs(s*(1.-x))) < 4.*Mpi){x = x - stepx;break;}
        //cout<<"x,step,cross ="<<","<<x<<", "<<stepx<< " , " << crvis_x <<endl;
    }
    double xmax = (s-16.*Mpi*Mpi)/s;
    if(x >= xmax)x = xmax;
    return x*E/2.;
}


double Cmd3Generator2pi2pi0_ke::generate_photon_angle(double E, double x){

    double xnew = gRandom->Rndm();
    xnew = xnew*radiatorW_th_integral(E, x);
    double result = 0.;
    int nrun = 10000000;
    double step = (double)2./nrun;
    double costh = -1.;
    int i = 0;
    double xA = costh;
    double xB = xA;
    while(result < xnew){
        xA = xB;
        xB = xA + step;
        costh = costh + step;
        result = result + 2.*(radiatorW_th(E,x,xA)+4.*radiatorW_th(E,x,0.5*xA+0.5*xB)+radiatorW_th(E,x,xB))/6.*step;
        if(fabs(costh) < 0.996)step = (double)2./1000;
        i++;
        //cout << result << " " << costh << " " << step << endl;
    }
    xnew = gRandom->Rndm();
    if(xnew > 0.5){costh = -costh;}
    //cout << "costh = " << costh << endl;
    return costh;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_a1_rhopi_pi(P,P1,P2,P3,P4) + H_a1_rhopi_pi(P,P3,P2,P1,P4); 
    H = H + H_a1_rhopi_pi(P,P1,P4,P3,P2) + H_a1_rhopi_pi(P,P3,P4,P1,P2);

    return H;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_a1_sigmapi_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_a1_sigmapi_pi(P,P1,P2,P3,P4) + H_a1_sigmapi_pi(P,P3,P2,P1,P4); 
    H = H + H_a1_sigmapi_pi(P,P1,P4,P3,P2) + H_a1_sigmapi_pi(P,P3,P4,P1,P2);

    return H;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_rho_f0_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_rho_f0(P,P1,P2,P3,P4) + H_rho_f0(P,P3,P2,P1,P4); 
    H = H + H_rho_f0(P,P1,P4,P3,P2) + H_rho_f0(P,P3,P4,P1,P2);

    return H;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_rho_f0_inter_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_rho_f0_inter(P,P1,P2,P3,P4) + H_rho_f0_inter(P,P3,P2,P1,P4); 
    H = H + H_rho_f0_inter(P,P1,P4,P3,P2) + H_rho_f0_inter(P,P3,P4,P1,P2);

    return H;
}
	   
TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_rho_sigma_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_rho_sigma(P,P1,P2,P3,P4) + H_rho_sigma(P,P3,P2,P1,P4); 
    H = H + H_rho_sigma(P,P1,P4,P3,P2) + H_rho_sigma(P,P3,P4,P1,P2);

    return H;
}
	      
TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_ph_sp_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_ph_sp(P,P1,P2,P3,P4) + H_ph_sp(P,P3,P2,P1,P4); 
    H = H + H_ph_sp(P,P1,P4,P3,P2) + H_ph_sp(P,P3,P4,P1,P2);

    return H;
}
	      
TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_a2_pi_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_a2_pi(P,P1,P2,P3,P4) + H_a2_pi(P,P3,P2,P1,P4); 
    H = H + H_a2_pi(P,P1,P4,P3,P2) + H_a2_pi(P,P3,P4,P1,P2);

    return H;
}

TLorentzVectorC Cmd3Generator2pi2pi0_ke::H_a1_rhopi_pi_c_II(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){
	TLorentzVectorC H = H_a1_rhopi_pi_II(P,P1,P2,P3,P4) + H_a1_rhopi_pi_II(P,P3,P2,P1,P4); 
    H = H + H_a1_rhopi_pi_II(P,P1,P4,P3,P2) + H_a1_rhopi_pi_II(P,P3,P4,P1,P2);

    return H;
}

// -----------------------------------------------------------------------
TLorentzVectorC Cmd3Generator2pi2pi0_ke::matrix(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = H_omega_pi0(P,P1,P2,P3,P4)*omegapi0_*complex<double>(cos(omegapi0phase_),sin(omegapi0phase_));
    complex<double> fact = a1pi_rhopi_*complex<double>(cos(a1pi_rhopi_phase),sin(a1pi_rhopi_phase))/(1.+a1pi_sigmapi_);
    if(a1pi_rhopi_!=0.)H = H + (H_a1_rhopi_pi(P,P1,P2,P3,P4) + H_a1_sigmapi_pi(P,P1,P2,P3,P4)*a1pi_sigmapi_*complex<double>(cos(a1pi_sigmapi_phase),sin(a1pi_sigmapi_phase)))*fact;
    if(rhof0_!=0.)H = H + H_rho_f0(P,P1,P2,P3,P4)*rhof0_*complex<double>(cos(rhof0phase_),sin(rhof0phase_));
    if(rhof0_inter!=0.)H = H + H_rho_f0_inter(P,P1,P2,P3,P4)*rhof0_inter*complex<double>(cos(rhof0phase_inter),sin(rhof0phase_inter));
    if(rhosigma_!=0.)H = H + H_rho_sigma(P,P1,P2,P3,P4)*rhosigma_*complex<double>(cos(rhosigmaphase_),sin(rhosigmaphase_));
    if(rhoprhom_II!=0.)H = H + (H_rhop_rhom(P,P1,P2,P3,P4)*rhoprhom_*complex<double>(cos(rhoprhom_phase),sin(rhoprhom_phase))+H_rhop_rhom_II(P,P1,P2,P3,P4))/(1.+rhoprhom_)*rhoprhom_II*complex<double>(cos(rhoprhom_phase_II),sin(rhoprhom_phase_II));
    if(a2pi_!=0.)H = H + H_a2_pi(P,P1,P2,P3,P4)*a2pi_*complex<double>(cos(a2pi_ph),sin(a2pi_ph));
    if(phsp_!=0.)H = H + H_ph_sp(P,P1,P2,P3,P4)*phsp_*complex<double>(cos(phspphase_),sin(phspphase_));
    if(h1pi!=0.)H = H + H_h1_rhopi(P,P1,P2,P3,P4)*h1pi*complex<double>(cos(h1pi_ph),sin(h1pi_ph));
    if(a1pi_rhopi_II!=0.)H = H + H_a1_rhopi_pi_II(P,P1,P2,P3,P4)*a1pi_rhopi_II*complex<double>(cos(a1pi_rhopi_phase_II),sin(a1pi_rhopi_phase_II));
    if(rhof2!=0.)H = H + H_rho_f2(P,P1,P2,P3,P4)*rhof2*complex<double>(cos(rhof2_phase),sin(rhof2_phase));
    
    return H;

}


double Cmd3Generator2pi2pi0_ke::matrix_squared(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P1,P2,P3,P4); 

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();

}

double Cmd3Generator2pi2pi0_ke::matrix_squared_c(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P1,P2,P3,P4) + matrix(P,P3,P2,P1,P4); 
    H = H + matrix(P,P1,P4,P3,P2) + matrix(P,P3,P4,P1,P2);

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();
}

double Cmd3Generator2pi2pi0_ke::matrix_squared_m000(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P2,P3,P4,P1) + matrix(P,P2,P3,P4,P1); 
    H = H + matrix(P,P1,P2,P4,P3);

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();
}

double Cmd3Generator2pi2pi0_ke::matrix_squared_m0001(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

   TLorentzVectorC Pa(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1a(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2a(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3a(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4a(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));
    return matrix_squared_m000(Pa,P1a,P2a,P3a,P4a);
}

double Cmd3Generator2pi2pi0_ke::matrix_squared_mmp0(TLorentzVectorC P, TLorentzVectorC P1, TLorentzVectorC P2, TLorentzVectorC P3, TLorentzVectorC P4){

    TLorentzVectorC H = matrix(P,P4,P2,P1,P3) + matrix(P,P4,P1,P2,P3); 

    return H.X().real()*H.X().real() + H.X().imag()*H.X().imag() + H.Y().real()*H.Y().real() + H.Y().imag()*H.Y().imag();
}

double Cmd3Generator2pi2pi0_ke::matrix_squared_mmp01(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

    TLorentzVectorC Pa(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1a(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2a(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3a(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4a(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));

    return matrix_squared_mmp0(Pa,P1a,P2a,P3a,P4a);
}



double Cmd3Generator2pi2pi0_ke::matrix_squared1(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

    TLorentzVectorC Pc(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1c(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2c(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3c(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4c(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));

    return matrix_squared(Pc,P1c,P2c,P3c,P4c);
}

double Cmd3Generator2pi2pi0_ke::matrix_squared1_c(TLorentzVector P, TLorentzVector P1, TLorentzVector P2, TLorentzVector P3, TLorentzVector P4){

    TLorentzVectorC Pc(complex<double>(P.X(),0.),complex<double>(P.Y(),0.),complex<double>(P.Z(),0.),complex<double>(P.E(),0.));
    TLorentzVectorC P1c(complex<double>(P1.X(),0.),complex<double>(P1.Y(),0.),complex<double>(P1.Z(),0.),complex<double>(P1.E(),0.));
    TLorentzVectorC P2c(complex<double>(P2.X(),0.),complex<double>(P2.Y(),0.),complex<double>(P2.Z(),0.),complex<double>(P2.E(),0.));
    TLorentzVectorC P3c(complex<double>(P3.X(),0.),complex<double>(P3.Y(),0.),complex<double>(P3.Z(),0.),complex<double>(P3.E(),0.));
    TLorentzVectorC P4c(complex<double>(P4.X(),0.),complex<double>(P4.Y(),0.),complex<double>(P4.Z(),0.),complex<double>(P4.E(),0.));

    return matrix_squared_c(Pc,P1c,P2c,P3c,P4c);
}

//**************************************************************************************************************************************************
void  Cmd3Generator2pi2pi0_ke::find_majoranta0(double E){
    int i = 0;
    TGenPhaseSpace event;
    TLorentzVector P(0.,0.,0., E);
    TLorentzVector* P1 = new TLorentzVector();
    TLorentzVector* P2 = new TLorentzVector();
    TLorentzVector* P3 = new TLorentzVector();
    TLorentzVector* P4 = new TLorentzVector();
    //TH1D *u = new TH1D("u","u",10000,-3.14,2.*3000.14);
    double masses[4] = {mpiz, mpiz, Mpi, Mpi};
    event.SetDecay(P, 4, masses);
    majoranta = 0.;
    while(i < 10000){
        double wes = event.Generate();
        P1 = event.GetDecay(0);
        P2 = event.GetDecay(1);
        P3 = event.GetDecay(2);
        P4 = event.GetDecay(3);
        double batva = matrix_squared1(P,*P1,*P2,*P3,*P4)*wes;
	if(batva > majoranta){majoranta = batva;}
        i++;
    }
    P1->Delete();
    P2->Delete();
    P3->Delete();
    P4->Delete();
    //cout << "Majoranta =  "<< majoranta << endl;
    majoranta = majoranta*1.1;
}

void  Cmd3Generator2pi2pi0_ke::find_min0(double E){
    int i = 0;
    TGenPhaseSpace event;
    TLorentzVector P(0.,0.,0., E);
    TLorentzVector* P1 = new TLorentzVector();
    TLorentzVector* P2 = new TLorentzVector();
    TLorentzVector* P3 = new TLorentzVector();
    TLorentzVector* P4 = new TLorentzVector();
    //TH1D *u = new TH1D("u","u",10000,-3.14,2.*3000.14);
    double masses[4] = {mpiz, mpiz, Mpi, Mpi};
    event.SetDecay(P, 4, masses);
    minmajoranta = 100000.;
    while(i < 100000){
        //cout << i << endl;
        event.Generate();
        P1 = event.GetDecay(0);
        P2 = event.GetDecay(1);
        P3 = event.GetDecay(2);
        P4 = event.GetDecay(3);
        double batva = matrix_squared1(P,*P1,*P2,*P3,*P4);
        if(batva < minmajoranta){minmajoranta = batva;}
        //u->Fill(batva);
        i++;
    }
    //u->Draw();
    P1->Delete();
    P2->Delete();
    P3->Delete();
    P4->Delete();
    //cout << "Majoranta =  "<< majoranta << endl;
    minmajoranta = -minmajoranta*4;
    if(minmajoranta<0)minmajoranta = 0;
}
//----------------------------------------------------------------------------------------------------------------------------------------------------

double cros = 0;
TGenPhaseSpace event;
TLorentzVector* P1 = new TLorentzVector();
TLorentzVector* P2 = new TLorentzVector();
TLorentzVector* P3 = new TLorentzVector();
TLorentzVector* P4 = new TLorentzVector();
TLorentzVector P(0.,0.,0., 0.);
double masses[4] = {mpiz, mpiz, Mpi, Mpi};
double masses_c[4] = {Mpi, Mpi, Mpi, Mpi};
int niter = 20000;

double mass3pi[400000];
double mass3pi_mc[400000];
TLorentzVectorC La1pi_rho_pi[400000];
TLorentzVectorC La1pi_rho_pi_II[400000];
TLorentzVectorC La1pi_rho_pi_II_mc[400000];
TLorentzVectorC La1pi_rho_pi_II_c[400000];
TLorentzVectorC La1pi_rho_pi_II_c_mc[400000];
TLorentzVectorC Lomega_pi[400000];
TLorentzVectorC La1pi_sigma_pi[400000];
TLorentzVectorC La1pi_rho_f0[400000];
TLorentzVectorC Lrho_f0_inter[400000];
TLorentzVectorC La1pi_rho_sigma[400000];
TLorentzVectorC L_rhop_rhom[400000];
TLorentzVectorC L_ph_sp[400000];
TLorentzVectorC L_a2pi[400000];
TLorentzVectorC L_h1pi[400000];
TLorentzVectorC L_rho_f2[400000];
TLorentzVectorC L_rhop_rhom_II[400000];

TLorentzVectorC La1pi_rho_pi_mc[400000];
TLorentzVectorC Lomega_pi_mc[400000];
TLorentzVectorC La1pi_sigma_pi_mc[400000];
TLorentzVectorC La1pi_rho_f0_mc[400000];
TLorentzVectorC Lrho_f0_inter_mc[400000];
TLorentzVectorC La1pi_rho_sigma_mc[400000];
TLorentzVectorC L_rhop_rhom_mc[400000];
TLorentzVectorC L_ph_sp_mc[400000];
TLorentzVectorC L_a2pi_mc[400000];
TLorentzVectorC L_h1pi_mc[400000];
TLorentzVectorC L_rho_f2_mc[400000];
TLorentzVectorC L_rhop_rhom_II_mc[400000];

TLorentzVectorC La1pi_rho_pi_c[400000];
TLorentzVectorC La1pi_sigma_pi_c[400000];
TLorentzVectorC La1pi_rho_f0_c[400000];
TLorentzVectorC Lrho_f0_inter_c[400000];
TLorentzVectorC La1pi_rho_sigma_c[400000];
TLorentzVectorC L_ph_sp_c[400000];
TLorentzVectorC L_a2pi_c[400000];
TLorentzVectorC L_rho_f2_c[400000];

TLorentzVectorC La1pi_rho_pi_c_mc[400000];
TLorentzVectorC La1pi_sigma_pi_c_mc[400000];
TLorentzVectorC La1pi_rho_f0_c_mc[400000];
TLorentzVectorC Lrho_f0_inter_c_mc[400000];
TLorentzVectorC La1pi_rho_sigma_c_mc[400000];
TLorentzVectorC L_ph_sp_c_mc[400000];
TLorentzVectorC L_a2pi_c_mc[400000];
TLorentzVectorC L_rho_f2_c_mc[400000];


vector<double> qwe;
vector<TLorentzVector> qwe1;
vector<TLorentzVector> qwe2;
vector<TLorentzVector> qwe3;
vector<TLorentzVector> qwe4;

vector<double> qwe_c;
vector<TLorentzVector> qwe1_c;
vector<TLorentzVector> qwe2_c;
vector<TLorentzVector> qwe3_c;
vector<TLorentzVector> qwe4_c;
vector<TLorentzVector> qwe1_c_mc;
vector<TLorentzVector> qwe2_c_mc;
vector<TLorentzVector> qwe3_c_mc;
vector<TLorentzVector> qwe4_c_mc;

double Nexp = 0;
double Nphsp = 0;

double Nexp2pi2pi0 = 0;
double Nphsp2pi2pi0 = 0;

double Nbckgr = 0;
double Nmc = 0;

int N_amplitudes = 27;


double energy_cm[100000];
double energy_cm_aver = 0;

double results[100];
double Ratio,dratio;



void Cmd3Generator2pi2pi0_ke::fill_sample(){

    qwe.clear();
    qwe1.clear();
    qwe2.clear();
    qwe3.clear();
    qwe4.clear();
    P.SetE(Energy);
    int i = 0;
    event.SetDecay(P, 4, masses);
    while(i < niter){
       double weight = event.Generate();
       qwe.push_back(weight);
       P1 = event.GetDecay(0);
       P2 = event.GetDecay(1);
       P3 = event.GetDecay(2);
       P4 = event.GetDecay(3);
       qwe1.push_back(*P1);
       qwe2.push_back(*P2);
       qwe3.push_back(*P3);
       qwe4.push_back(*P4);
       i++;
    }
}

void Cmd3Generator2pi2pi0_ke::fill_sample_c(){

    qwe_c.clear();
    qwe1_c.clear();
    qwe2_c.clear();
    qwe3_c.clear();
    qwe4_c.clear();
    P.SetE(Energy);
    int i = 0;
    event.SetDecay(P, 4, masses_c);
    while(i < niter){
       double weight = event.Generate();
       qwe_c.push_back(weight);
       P1 = event.GetDecay(0);
       P2 = event.GetDecay(1);
       P3 = event.GetDecay(2);
       P4 = event.GetDecay(3);
       qwe1_c.push_back(*P1);
       qwe2_c.push_back(*P2);
       qwe3_c.push_back(*P3);
       qwe4_c.push_back(*P4);
       i++;
    }
/*
   TH1D *hE = new TH1D("hE","hE",1000,0,1.);
   double mass[] = {0.1,0.1};
   double en = 1.;
   P.SetE(en);
   event.SetDecay(P, 2, mass);
   i = 0;
   while(i < 10000){
	event.Generate();
        hE->Fill(event.GetDecay(0)->P());
        i++;
   }
   hE->Draw();
*/
}



double Cmd3Generator2pi2pi0_ke::cross_section(){

   P.SetE(Energy);
   cros = 0;
   double cross_m000 = 0.;
   double cross_mmp0 = 0.;   
   int i = 0;
   double diffcr[100];
   double diffcros = 0.;
   for(int i = 0; i < 100; i++){cross_s[i] = 0.;diffcr[i] = 0.;dcross_s[i] = 0.;}
   TLorentzVectorC CV,CV1,CV2,CV3,CV4;
   
   while(i < niter){
     if(i%1000==0)cout << i << endl;
     CV = qwe1[i]+qwe2[i]+qwe3[i]+ qwe4[i];
     CV1 = qwe1[i];
     CV2 = qwe2[i];
     CV3 = qwe3[i];
     CV4 = qwe4[i];
     double matr2 = matrix_squared1(P, qwe1[i],qwe2[i],qwe3[i], qwe4[i]);
     double matr_mmp0 = matrix_squared_mmp01(P, qwe1[i],qwe2[i],qwe3[i], qwe4[i]);
     double matr_m000 = matrix_squared_m0001(P, qwe1[i],qwe2[i],qwe3[i], qwe4[i]);
     
     cros = cros + matr2*qwe[i];
     cross_m000 = cros + matr_mmp0*qwe[i];
     cross_mmp0 = cros + matr_mmp0*qwe[i];
     
     matr2 = mod(H_omega_pi0(CV,CV1,CV2,CV3,CV4)*omegapi0_);
     cross_s[0] = cross_s[0] + matr2*qwe[i];
     
     complex<double> fact = a1pi_rhopi_/(1.+a1pi_sigmapi_);
     matr2 = mod((H_a1_rhopi_pi(CV,CV1,CV2,CV3,CV4) + H_a1_sigmapi_pi(CV,CV1,CV2,CV3,CV4)*a1pi_sigmapi_*complex<double>(cos(a1pi_sigmapi_phase), sin(a1pi_sigmapi_phase)))*fact);
     cross_s[1] = cross_s[1] + matr2*qwe[i]; 
     
     matr2 = mod(H_rho_f0(CV,CV1,CV2,CV3,CV4)*rhof0_*complex<double>(cos(rhof0phase_), sin(rhof0phase_)) + H_rho_sigma(CV,CV1,CV2,CV3,CV4)*rhosigma_*complex<double>(cos(rhosigmaphase_), sin(rhosigmaphase_)));
     cross_s[2] = cross_s[2] + matr2*qwe[i]; 
     
     matr2 = mod((H_rhop_rhom(CV,CV1,CV2,CV3,CV4)*rhoprhom_*complex<double>(cos(rhoprhom_phase),sin(rhoprhom_phase)) +H_rhop_rhom_II(CV,CV1,CV2,CV3,CV4))/(1.+rhoprhom_)*rhoprhom_II); 
     cross_s[3] = cross_s[3] + matr2*qwe[i]; 
     
     matr2 = mod(H_rho_sigma(CV,CV1,CV2,CV3,CV4)*rhosigma_); 
     cross_s[4] = cross_s[4] + matr2*qwe[i];
     
     matr2 = mod(H_h1_rhopi(CV,CV1,CV2,CV3,CV4)*h1pi);
     cross_s[5] = cross_s[5] + matr2*qwe[i];
     
     matr2 = mod(H_rho_f2(CV,CV1,CV2,CV3,CV4)*rhof2);
     cross_s[6] = cross_s[6] + matr2*qwe[i];
     i++;
   }
   i = 0;

   if(a1pi_rhopi_!=0.){
     dcross_s[1] = a1pi_rhopi_;
     change(2,-1000);
     dcross_s[1] = dcross_s[1] - a1pi_rhopi_;
     dcross_s[1] = dcross_s[1]/a1pi_rhopi_;
   }
   if(rhof0_!=0){
     dcross_s[2] = rhof0_;
     change(6,-1000);
     dcross_s[2] =  dcross_s[2] - rhof0_; 
     dcross_s[2] = dcross_s[2]/rhof0_;
   }
   if(rhoprhom_!=0){
     dcross_s[3] = rhoprhom_;
     change(8,-1000);
     dcross_s[3] = dcross_s[3] - rhoprhom_;
     dcross_s[3] = dcross_s[3]/rhoprhom_;
   }
   if(h1pi!=0){
     dcross_s[5] = h1pi;
     change(16,-1000);
     dcross_s[5] = dcross_s[5] - h1pi;
     dcross_s[5] = dcross_s[5]/h1pi;
   }
   if(rhof2!=0){
     dcross_s[6] = rhof2;
     change(22,-1000);
     dcross_s[6] = dcross_s[6] - rhof2;
     dcross_s[6] = dcross_s[6]/rhof2;
   }
   for(int l = 2; l < 25; l = l + 2){
     fill_parameters();
     change(l,-1000);
     diffcr[l] = 0.;
     //cout << "l = " << l << endl;
     while(i < niter){
       diffcr[l] = diffcr[l]+ matrix_squared1(P, qwe1[i],qwe2[i],qwe3[i], qwe4[i])*qwe[i];
       i++;
     }
     i = 0;
     //cout << "diffcr[l] = " << diffcr[l] << " cros = " << cros << endl;
     diffcros = diffcros + pow(diffcr[l] - cros,2.);
   }
   diffcros = sqrt(diffcros);
   for(int i = 0; i <= 6; i++){
     cross_s[i] = cross_s[i]/cros;
   }
   dcross_s[0] = cross_s[0]*diffcros/cros;
   for(int i = 1; i <= 6; i++){
     dcross_s[i] = sqrt(pow(dcross_s[i],2.) + pow(diffcros/cros,2.));
     dcross_s[i] = dcross_s[i]*cross_s[i];
   }

   for(int i = 0; i <= 6; i++){
     cout << " cross_s[i] = " << cross_s[i] << " +/- " << dcross_s[i]<< endl;
   }
   
   cross_s[20] = cross_m000/cros;
   cross_s[21] = cross_mmp0/cros;   
   return 0.5*cros/niter;

}



double Cmd3Generator2pi2pi0_ke::cross_section_c(){

   P.SetE(Energy);
   cros = 0;
   int i = 0;
   TRandom1 rand;
   double xnew = rand.Rndm();
   //find_majoranta0(Energy);
   double here = 0;
   while(i < niter){
        double matr2 = matrix_squared1_c(P, qwe1_c[i],qwe2_c[i],qwe3_c[i], qwe4_c[i]);
        xnew = rand.Rndm();
        if(1==1/*matr2 > xnew*majoranta/qwe1_c[i].Theta() > 0.4 && qwe1_c[i].Theta() < 3.1415 - 0.4 &&
                qwe2_c[i].Theta() > 0.4 && qwe2_c[i].Theta() < 3.1415 - 0.4 &&
                qwe3_c[i].Theta() > 0.8 && qwe3_c[i].Theta() < 3.1415 - 0.8 &&
                qwe4_c[i].Theta() > 0.8 && qwe4_c[i].Theta() < 3.1415 - 0.8*/){
            cros = cros + matr2*qwe_c[i];
            here++;
        }
        i++;
   }
   i = 0;
   return 0.25*cros/niter;
}


void Cmd3Generator2pi2pi0_ke::readd(double minenergy, double maxenergy, string from = "mc_phsp/"){
        
    cout << "--------- The sample of 2pi2pi0 events --------------- " << endl;
    cout << "minenergy = " << minenergy << endl;
    cout << "maxenergy = " << maxenergy << endl;
    energy_cm_aver = 0;
    Nphsp2pi2pi0 = 0;
    int i = 0;
    TLorentzVectorC P4p[4];
    int charge[4];
    double beam, A, B, C1;
    //TH1D *minvmc = new TH1D("minvmc","minvmc",400,0,2);
    // TH1D *minv_exp = new TH1D("minv_exp","minv_exp",400,0,2);
    
    DIR *dir = opendir(("../4pisel/histograms/"+from).c_str());
    if(dir)
    {
        struct dirent *ent;
        while((ent = readdir(dir)) != NULL)
        {
            char *dd = ent->d_name;
            if(dd[strlen(dd)-2]!='a')continue;
            std::stringstream ss;
            ss << dd[3] << dd[4] << dd[5] << dd[6] << dd[7];
            TString s = ss.str();
            beam = s.Atof();
            beam = 2.*beam/1000.;
            string names = ent->d_name;
            double wind = 0.;
            //if(minenergy > 1.875 && minenergy < 1.885)wind = 0.006;
            if(beam < minenergy-0.001 || beam > maxenergy + 0.001){continue;}
            ifstream streamNphsp(("../4pisel/histograms/" + from + names).c_str());
            cout << ("../4pisel/histograms/" + from + names).c_str() << endl;
	    while(streamNphsp.eof() == 0){
	      
            	    	streamNphsp >> beam;
            	    	beam = 2.*beam/1000.;
            	    	for(int k = 0; k < 4; k++){
            	    	    streamNphsp >> A >> B >> C1;
            	    	}
			
            	    	streamNphsp >> beam;
            	    	beam = 2.*beam/1000.;
            	    	for(int k = 0; k < 4; k++){
            	    	    streamNphsp >> A >> B >> C1;
            	    	    if(k < 2)P4p[k] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 + mpiz*mpiz));
            	    	    else P4p[k] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 +Mpi*Mpi));
            	    	}

            	    	int onep = 0;
            	    	int onen = 1;
            	    	int twop = 2;
            	    	int twon = 3;
			/*double m =  (P4p[1] + P4p[2] + P4p[3]).M();
			if(m > 0.75 && m < 0.825)continue;
 			m =  (P4p[0] + P4p[2] + P4p[3]).M();
			if(m > 0.75 && m < 0.825)continue;*/
                        
                        double betax = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Px().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
                        double betay = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Py().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
                        double betaz = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Pz().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
			//P4p[onep].Boost(-betax,-betay,-betaz);
			//P4p[onen].Boost(-betax,-betay,-betaz);
			//P4p[twop].Boost(-betax,-betay,-betaz);
			//P4p[twon].Boost(-betax,-betay,-betaz);
            	    	La1pi_rho_pi_mc[(int)Nphsp2pi2pi0] = H_a1_rhopi_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	Lomega_pi_mc[(int)Nphsp2pi2pi0] = H_omega_pi0(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	La1pi_sigma_pi_mc[(int)Nphsp2pi2pi0] = H_a1_sigmapi_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	La1pi_rho_f0_mc[(int)Nphsp2pi2pi0] = H_rho_f0(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
			Lrho_f0_inter_mc[(int)Nphsp2pi2pi0] = H_rho_f0_inter(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	La1pi_rho_sigma_mc[(int)Nphsp2pi2pi0] = H_rho_sigma(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	L_rhop_rhom_mc[(int)Nphsp2pi2pi0] = H_rhop_rhom(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	L_ph_sp_mc[(int)Nphsp2pi2pi0] = H_ph_sp(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	L_a2pi_mc[(int)Nphsp2pi2pi0] = H_a2_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	L_h1pi_mc[(int)Nphsp2pi2pi0] = H_h1_rhopi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	La1pi_rho_pi_II_mc[(int)Nphsp2pi2pi0] = H_a1_rhopi_pi_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
			L_rho_f2_mc[(int)Nphsp2pi2pi0] = H_rho_f2(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
		        L_rhop_rhom_II_mc[(int)Nphsp2pi2pi0] = H_rhop_rhom_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            	    	energy_cm_aver = energy_cm_aver + beam;
            	    	//minvmc->Fill((P4p[onep]+P4p[twon]+P4p[twop]).M(), Lomega_pi_mc[i].modX() + Lomega_pi_mc[i].modY());
            	    	Nphsp2pi2pi0++;
			if(Nphsp2pi2pi0 > 80000)break;
            	    }
            	    //if(Nphsp2pi2pi0 > 40000 && beam < 1.2)break;
        }
   }
    //minvmc->Draw();
    energy_cm_aver = energy_cm_aver/Nphsp2pi2pi0;
    cout << "energy_cm_aver in MC = " << energy_cm_aver << endl;
    cout << "Nphsp_2pi2pi0 = " << Nphsp2pi2pi0 << endl;
    Nexp2pi2pi0 = 0;
    energy_cm_aver = 0;
    
    ifstream streamNexp(if_name_AA_data.c_str());
    Nexp2pi2pi0 = 0.;
    i = 0;
    while(streamNexp.eof() == 0){
        streamNexp >> beam;
        beam = 2.*beam/1000.;
	//cout << minenergy << " " << beam << " " << maxenergy << endl;
	//	getchar();
        if(beam < maxenergy && beam > minenergy){

	    if(TOY == true)streamNexp >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A;
	  
            for(int i = 0; i < 4; i++){
              streamNexp >> A >> B >> C1;
              if(i < 2)P4p[i] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 + mpiz*mpiz));
              else P4p[i] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 + Mpi*Mpi));
            }

            int onep = 0;
            int onen = 1;
            int twop = 2;
            int twon = 3;
 	    /*double m =  (P4p[1] + P4p[2] + P4p[3]).M();
	    if(m > 0.75 && m < 0.825)continue;
 	    m =  (P4p[0] + P4p[2] + P4p[3]).M();
	    if(m > 0.75 && m < 0.825)continue;*/	    
            //cout << "4pi_data = " << (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E() << endl;
            double betax = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Px().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
            double betay = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Py().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
            double betaz = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Pz().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	    //P4p[onep].Boost(-betax,-betay,-betaz);
	    //P4p[onen].Boost(-betax,-betay,-betaz);
	    //P4p[twop].Boost(-betax,-betay,-betaz);
	    //P4p[twon].Boost(-betax,-betay,-betaz);
            La1pi_rho_pi[i] = H_a1_rhopi_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            Lomega_pi[i] = H_omega_pi0(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            La1pi_sigma_pi[i] = H_a1_sigmapi_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            La1pi_rho_f0[i] = H_rho_f0(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	    Lrho_f0_inter[i] = H_rho_f0_inter(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            La1pi_rho_sigma[i] = H_rho_sigma(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            L_rhop_rhom[i] = H_rhop_rhom(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            L_ph_sp[i] = H_ph_sp(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            L_a2pi[i] = H_a2_pi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            L_h1pi[i] = H_h1_rhopi(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            La1pi_rho_pi_II[i] = H_a1_rhopi_pi_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	    L_rho_f2[i] = H_rho_f2(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	    L_rhop_rhom_II[i] = H_rhop_rhom_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
            Nexp2pi2pi0++;
            //minv_exp->Fill((P4p[onep]+P4p[twon]+P4p[twop]).M());//, Lomega_pi_mc[i].modX() + Lomega_pi_mc[i].modY());
            energy_cm_aver = energy_cm_aver + beam;
            i++;
          }
          else{
            streamNexp >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A >> A >> B >> C1;
          }
     }
    //minvmc->SetNormFactor(minv_exp->GetEntries()); 
    //minv_exp->Draw("e");
    //minvmc->Draw("same");
    energy_cm_aver = energy_cm_aver/i;
    cout << "energy_cm_aver in EXP = " << energy_cm_aver << endl;
    cout << "Nexp_2pi2pi0 = " << Nexp2pi2pi0 << endl;
    cout << "-------------------------------------------------" << endl;
}


void Cmd3Generator2pi2pi0_ke::readd_c(double minenergy, double maxenergy, string from = "mc_phsp/"){
    cout << "----------------- The sample of 2pi+ 2pi- events ------------------ " << endl;
    cout << "minenergy = " << minenergy << endl;
    cout << "maxenergy = " << maxenergy << endl;
    energy_cm_aver = 0;
    //fourmomentum4pi2011mcafterkinfit.txt");//matrix_700_2011.dat");
    Nphsp = 0;
    int i = 0;
    TLorentzVectorC P4p[4];
    int charge[4];
    double beam, A, B, C1;
    TH1D *hener = new TH1D("hener","hener",10000,0,2);
    DIR *dir = opendir(("../2pip2pim/histograms/" + from).c_str());
    if(dir)
    {
        struct dirent *ent;
        while((ent = readdir(dir)) != NULL)
        {

            char *dd = ent->d_name;
            if(dd[strlen(dd)-2]!='a')continue;
            std::stringstream ss;
            ss << dd[7] << dd[8] << dd[9] << dd[10] << dd[11];
            TString s = ss.str();
            beam = s.Atof();
            beam = 2.*beam/1000.;
            string names = ent->d_name;
            double wind = 0.;
            if(beam < minenergy-0.001 || beam > maxenergy+0.001){continue;}
            ifstream streamNphsp(("../2pip2pim/histograms/" +from  + names).c_str());
	    cout << ("../2pip2pim/histograms/" + from + names).c_str() << endl;
            while(streamNphsp.eof() == 0){

	      streamNphsp >> A;
	      for(int k = 0; k < 4; k++){
		streamNphsp >> A >> B >> C1;
	      }
	      
	      streamNphsp >> beam;
	      beam = 2.*beam/1000.;
	      //cout << beam << " ";
	      for(int k = 0; k < 4; k++){
		streamNphsp >> A >> B >> C1;
		P4p[k] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 + Mpi*Mpi));
	      }
	      
	      if(streamNphsp.eof() == 1)break;
	      int onep = 0;
	      int twop = 1;
	      int onen = 2;
	      int twon = 3;
	      // 0 2 1 3
	      energy_cm_aver = energy_cm_aver + beam;
	      hener->Fill((P4p[onep]+P4p[twop]).M(),La1pi_rho_pi_c_mc[i].modX() + La1pi_rho_pi_c_mc[i].modY());
	      double betax = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Px().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	      double betay = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Py().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	      double betaz = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Pz().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	      //P4p[onep].Boost(-betax,-betay,-betaz);
	      //P4p[onen].Boost(-betax,-betay,-betaz);
	      //P4p[twop].Boost(-betax,-betay,-betaz);
	      //P4p[twon].Boost(-betax,-betay,-betaz);
	      
	      La1pi_rho_pi_c_mc[(int)Nphsp] = H_a1_rhopi_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      La1pi_sigma_pi_c_mc[(int)Nphsp] = H_a1_sigmapi_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      La1pi_rho_f0_c_mc[(int)Nphsp] = H_rho_f0_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      Lrho_f0_inter_c_mc[(int)Nphsp] = H_rho_f0_inter_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      La1pi_rho_sigma_c_mc[(int)Nphsp] = H_rho_sigma_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      L_ph_sp_c_mc[(int)Nphsp] = H_ph_sp_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      L_a2pi_c_mc[(int)Nphsp] = H_a2_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      La1pi_rho_pi_II_c_mc[(int)Nphsp] = H_a1_rhopi_pi_c_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      L_rho_f2_c_mc[(int)Nphsp] = H_rho_f2(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	      Nphsp++;
	      if(Nphsp > 80000)break;
	    }
        }
    }
    
    
    energy_cm_aver = energy_cm_aver/Nphsp;
    cout << "energy_cm_aver in MC = " << energy_cm_aver << endl;
    cout << "Nphsp_2pi+2pi- = " << Nphsp << endl;

    energy_cm_aver = 0;
    //TH1D *minv = new TH1D("minv","minv",400,0,2);
    ifstream streamNexp(if_name_AA_data_c.c_str()); 
    Nexp = 0.;
    i = 0;
     
    while(streamNexp.eof() == 0){
      
      if(TOY==true)streamNexp >> A >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A >> A >> B >> C1;
      
      streamNexp >> beam;
      beam = 2.*beam/1000.;
      if(beam < maxenergy && beam > minenergy){
	for(int k = 0; k < 4; k++){
	  streamNexp >> A >> B >> C1;
	  P4p[k] = TLorentzVectorC(A,B,C1,sqrt(A*A+B*B+C1*C1 + Mpi*Mpi));
	}
	
	int onep = 0;
	int twop = 1;
	int onen = 2;
	int twon = 3;
	double betax = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Px().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	double betay = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Py().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	double betaz = (P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).Pz().real()/(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon]).E().real();
	//P4p[onep].Boost(-betax,-betay,-betaz);
	//P4p[onen].Boost(-betax,-betay,-betaz);
	//P4p[twop].Boost(-betax,-betay,-betaz);
	//P4p[twon].Boost(-betax,-betay,-betaz);
	La1pi_rho_pi_c[i] = H_a1_rhopi_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	La1pi_sigma_pi_c[i] = H_a1_sigmapi_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	La1pi_rho_f0_c[i] = H_rho_f0_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	Lrho_f0_inter_c[(int)Nphsp] = H_rho_f0_inter_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	La1pi_rho_sigma_c[i] = H_rho_sigma_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	L_ph_sp_c[i] = H_ph_sp_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	L_a2pi_c[i] = H_a2_pi_c(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	La1pi_rho_pi_II_c[i] = H_a1_rhopi_pi_c_II(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	L_rho_f2_c[i] = H_rho_f2(P4p[onep]+P4p[onen]+P4p[twop]+P4p[twon],P4p[onep],P4p[onen],P4p[twop],P4p[twon]);
	energy_cm_aver = energy_cm_aver + beam;
	Nexp++;
	i++;
      }
      else{
	streamNexp >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A >> A >> B >> C1;
      }
    }
    
    energy_cm_aver = energy_cm_aver/i;
    cout << "energy_cm_aver in EXP= " << energy_cm_aver << endl;
    cout << "Nexp_2pi+2pi- = " << Nexp << endl;
    
}




double Likelyhood_4pic(double *par){

    double likelihood = 0.;
    double likelihoodphsp = 0.;

    double likelihood2pi2pi0 = 0.;
    double likelihoodphsp2pi2pi0 = 0.;

    for(int i = 0; i < Nexp; i++){
        TLorentzVectorC btv(0.,0.,0.,0.);
        complex<double> fact = par[2]*complex<double>(cos(par[3]),sin(par[3]))/(1.+par[4]);
        btv = (La1pi_rho_pi_c[i] + La1pi_sigma_pi_c[i]*par[4]*complex<double>(cos(par[5]),sin(par[5])))*fact;
        btv = btv + La1pi_rho_f0_c[i]*par[6]*complex<double>(cos(par[7]),sin(par[7]));
        btv = btv + La1pi_rho_sigma_c[i]*par[10]*complex<double>(cos(par[11]),sin(par[11]));
        btv = btv + L_ph_sp_c[i]*par[12]*complex<double>(cos(par[13]),sin(par[13]));
        btv = btv + L_a2pi_c[i]*par[14]*complex<double>(cos(par[15]),sin(par[15]));
        btv = btv + La1pi_rho_pi_II_c[i]*par[18]*complex<double>(cos(par[19]),sin(par[19]));
	btv = btv + Lrho_f0_inter_c[i]*par[20]*complex<double>(cos(par[21]),sin(par[21]));
        btv = btv + L_rho_f2_c[i]*par[22]*complex<double>(cos(par[23]),sin(par[23]));
	if((btv.modX() + btv.modY()) == 0)continue;
        likelihood = likelihood + log(btv.modX() + btv.modY());
	Likexp_c[i] =  log(btv.modX() + btv.modY());
   }

    for(int i = 0; i < Nphsp; i++){
        TLorentzVectorC btv(0.,0.,0.,0.);
        complex<double> fact = par[2]*complex<double>(cos(par[3]),sin(par[3]))/(1.+par[4]);
        btv = (La1pi_rho_pi_c_mc[i] + La1pi_sigma_pi_c_mc[i]*par[4]*complex<double>(cos(par[5]),sin(par[5])))*fact;
        btv = btv + La1pi_rho_f0_c_mc[i]*par[6]*complex<double>(cos(par[7]),sin(par[7]));
        btv = btv + La1pi_rho_sigma_c_mc[i]*par[10]*complex<double>(cos(par[11]),sin(par[11]));
        btv = btv + L_ph_sp_c_mc[i]*par[12]*complex<double>(cos(par[13]),sin(par[13]));
        btv = btv + L_a2pi_c_mc[i]*par[14]*complex<double>(cos(par[15]),sin(par[15]));
        btv = btv + La1pi_rho_pi_II_c_mc[i]*par[18]*complex<double>(cos(par[19]),sin(par[19]));
	btv = btv + Lrho_f0_inter_c_mc[i]*par[20]*complex<double>(cos(par[21]),sin(par[21]));
        btv = btv + L_rho_f2_c_mc[i]*par[22]*complex<double>(cos(par[23]),sin(par[23]));
	if((btv.modX() + btv.modY()) == 0)continue;
        likelihoodphsp = likelihoodphsp + btv.modX() + btv.modY();
	Likmc_c[i] =   log(btv.modX() + btv.modY());
    }

    for(int i = 0; i < Nexp2pi2pi0; i++){

        TLorentzVectorC btv = Lomega_pi[i]*par[1];
        complex<double> fact = par[2]*complex<double>(cos(par[3]),sin(par[3]))/(1.+par[4]);
        btv = btv + (La1pi_rho_pi[i] + La1pi_sigma_pi[i]*par[4]*complex<double>(cos(par[5]),sin(par[5])))*fact;
        btv = btv + La1pi_rho_f0[i]*par[6]*complex<double>(cos(par[7]),sin(par[7]));
        btv = btv + (L_rhop_rhom[i]*par[8]*complex<double>(cos(par[9]),sin(par[9])) +
L_rhop_rhom_II[i])*par[24]*complex<double>(cos(par[25]),sin(par[25]))/(1.+par[8]);
        btv = btv + La1pi_rho_sigma[i]*par[10]*complex<double>(cos(par[11]),sin(par[11]));
        btv = btv + L_ph_sp[i]*par[12]*complex<double>(cos(par[13]),sin(par[13]));
        btv = btv + L_a2pi[i]*par[14]*complex<double>(cos(par[15]),sin(par[15]));
        btv = btv + L_h1pi[i]*par[16]*complex<double>(cos(par[17]),sin(par[17]));        
        btv = btv + La1pi_rho_pi_II[i]*par[18]*complex<double>(cos(par[19]),sin(par[19]));
	btv = btv + Lrho_f0_inter[i]*par[20]*complex<double>(cos(par[21]),sin(par[21]));
        btv = btv + L_rho_f2[i]*par[22]*complex<double>(cos(par[23]),sin(par[23]));
	if((btv.modX() + btv.modY()) == 0)continue;
        likelihood2pi2pi0 = likelihood2pi2pi0 + log(btv.modX() + btv.modY());
	Likexp[i] =  log(btv.modX() + btv.modY());
    }

    for(int i = 0; i < Nphsp2pi2pi0; i++){
        TLorentzVectorC btv = Lomega_pi_mc[i]*par[1];
        complex<double> fact = par[2]*complex<double>(cos(par[3]),sin(par[3]))/(1.+par[4]);
        btv = btv + (La1pi_rho_pi_mc[i] + La1pi_sigma_pi_mc[i]*par[4]*complex<double>(cos(par[5]),sin(par[5])))*fact;
        btv = btv + La1pi_rho_f0_mc[i]*par[6]*complex<double>(cos(par[7]),sin(par[7]));
        btv = btv + (L_rhop_rhom_mc[i]*par[8]*complex<double>(cos(par[9]),sin(par[9])) +
L_rhop_rhom_II_mc[i])*par[24]*complex<double>(cos(par[25]),sin(par[25]))/(1.+par[8]);
        btv = btv + La1pi_rho_sigma_mc[i]*par[10]*complex<double>(cos(par[11]),sin(par[11]));
        btv = btv + L_ph_sp_mc[i]*par[12]*complex<double>(cos(par[13]),sin(par[13]));
        btv = btv + L_a2pi_mc[i]*par[14]*complex<double>(cos(par[15]),sin(par[15]));
        btv = btv + L_h1pi_mc[i]*par[16]*complex<double>(cos(par[17]),sin(par[17]));
        btv = btv + La1pi_rho_pi_II_mc[i]*par[18]*complex<double>(cos(par[19]),sin(par[19]));
	btv = btv + Lrho_f0_inter_mc[i]*par[20]*complex<double>(cos(par[21]),sin(par[21]));
        btv = btv + L_rho_f2_mc[i]*par[22]*complex<double>(cos(par[23]),sin(par[23]));
	if((btv.modX() + btv.modY()) == 0)continue;
        likelihoodphsp2pi2pi0 = likelihoodphsp2pi2pi0 + btv.modX() + btv.modY();
	Likmc[i] =   log(btv.modX() + btv.modY());
    }

    likelihoodphsp = likelihoodphsp/4.;
    likelihoodphsp2pi2pi0 = likelihoodphsp2pi2pi0/2.;
    double ress =  0.;
    if(Nphsp > 0)ress = ress - likelihood + Nexp*log(likelihoodphsp/Nphsp);
    if(Nphsp2pi2pi0 > 0)ress = ress - likelihood2pi2pi0 + Nexp2pi2pi0*log(likelihoodphsp2pi2pi0/Nphsp2pi2pi0);

    return ress;
}



double Cmd3Generator2pi2pi0_ke::LIKE(){

	double par[] = {0.,1.,a1pi_rhopi_,a1pi_rhopi_phase,a1pi_sigmapi_,a1pi_sigmapi_phase,rhof0_, rhof0phase_,rhoprhom_,rhoprhom_phase,rhosigma_,rhosigmaphase_,phsp_, phspphase_,a2pi_,a2pi_ph,h1pi,h1pi_ph,a1pi_rhopi_II, a1pi_rhopi_phase_II,rhof0_inter,rhof0phase_inter,0.,0.,0.,0.};
	return Likelyhood_4pic(par);

}



double drawlikelihood(double *s, double *par){
    double parttt[] = {par[0],par[1],par[2],par[3],s[0],par[4]};
    return Likelyhood_4pic(parttt);
}

double correctphase(double ph){
   int n = fabs((int)(ph/(2.*TMath::Pi())));
   if(ph >=0) ph = ph - n*2.*TMath::Pi();
   if(ph < 0) ph = ph + n*2.*TMath::Pi();
   if(ph > TMath::Pi())ph = ph - 2.*TMath::Pi();
   if(ph < -TMath::Pi())ph = 2.*TMath::Pi() + ph;
   return ph;
}

//TF1 *u = new TF1("u",drawlikelihood,0,1,1);u->SetParameter(0,1);u->SetTitle("L(a1=1,a2)");u->SetXTitle("a2");u->Draw();
//TF2 *u = new TF2("u",drawlikelihood,0.8,1,0,0.1,2);u->SetParameter(0,1);u->Draw("surf2z");
//TF1 *u1 = new TF1("u1",drawlikelihood,-10,10,2);u1->SetParameters(0,1);u1->Draw();
//TF1 *u1 = new TF1("u1",drawlikelihood,-15,15,5);u1->SetParameters(0,1,0,0,0);u1->Draw();

void minuitCMD3(int &nDim,double*gout,double&result,double par[],int flg){
    result=Likelyhood_4pic(par);}

int minimimlik(double minenergy, double maxenergy){

    Part_ana.readd_c(minenergy, maxenergy);
    Part_ana.readd(minenergy, maxenergy);
    if(energy_cm_aver < 0.3 ||/* Nexp < 100 || Nphsp < 100 || */Nexp2pi2pi0 < 100 || Nphsp2pi2pi0 < 100)return 0;
    Part_ana.Energy = energy_cm_aver;
    TFitter*minimizer=new TFitter(26);
    double minimum = 0;
    double resultssss[1000];
    int k = 0;TRandom rnm;
    
    while(k < 10){

    minimizer->SetFitMethod("loglikelihood");
    minimizer->SetParameter(0,"lambda ",                0, 0.0001,-10000,10000.2);
    minimizer->SetParameter(1,"ompi0 ",                     1,   0.001,-0.01,200);
    minimizer->SetParameter(2,"a1 (rho pi) ",               rnm.Rndm(),  0.001,0,20);
    minimizer->SetParameter(3,"a1 (rho pi) phase ",         rnm.Rndm()*6, 0.001,-12.,20);
    minimizer->SetParameter(4,"a1 (sigma pi) ",             rnm.Rndm(),   0.001,0,6);
    minimizer->SetParameter(5,"a1 (sigma pi) phase ",       rnm.Rndm()*6,  0.001,-10,10);
    minimizer->SetParameter(6,"rho f0",                     rnm.Rndm(),   0.001,0.,50);
    minimizer->SetParameter(7,"rho f0 phase ",              rnm.Rndm()*6,  0.001,-6.7,20);
    minimizer->SetParameter(8,"rhop rhom",                  0.,   0.0001,0.,100);
    minimizer->SetParameter(9,"rhop rhom phase ",           rnm.Rndm()*6,  0.001,-6.7,20);
    minimizer->SetParameter(10,"rho sigma",                 rnm.Rndm(),   0.0001,0.,60);
    minimizer->SetParameter(11,"rho sigma phase ",          rnm.Rndm()*6,  0.001,-10.7,20);
    minimizer->SetParameter(12,"phase sp",                  0.,   0.0001,0.,10);
    minimizer->SetParameter(13,"phase sp phase ",           0.,  0.001,0,20.);
    minimizer->SetParameter(14,"a2 pi",                     0.,   0.001,0.,10);
    minimizer->SetParameter(15,"a2 pi phase ",              0.,  0.001,-6.7,20.);
    minimizer->SetParameter(16,"h1 pi",                     0.,   0.001,0.,10);
    minimizer->SetParameter(17,"h1 pi phase ",              rnm.Rndm()*6,  0.001,-10.,22);
    minimizer->SetParameter(18,"a1 (rho pi) pi_II",         0.,  0.001,0,300);
    minimizer->SetParameter(19,"a1 (rho pi) pi_II phase ",  0.,  0.001,-8.,3.2*2.);
    minimizer->SetParameter(20,"rho f0_inter",              0.,   0.001,-10,100);
    minimizer->SetParameter(21,"rho f0_inter phase ",       0.,  0.001,-6.7,20);
    minimizer->SetParameter(22,"rho f2",	            0.,   0.001,0.,100);
    minimizer->SetParameter(23,"rho f2 phase ", 	    0.,  0.001,-6.7,20);
    minimizer->SetParameter(24,"rhop rhom II ",             0.,   0.001,0.,100);
    minimizer->SetParameter(25,"rhop rhom II phase ",       0.,  0.001,-6.7,20);

    
    minimizer->FixParameter(0);
    minimizer->FixParameter(1);
    //minimizer->FixParameter(2);
    //minimizer->FixParameter(3);
    //minimizer->FixParameter(4);
    //minimizer->FixParameter(5);
    //minimizer->FixParameter(6);
    //minimizer->FixParameter(7);
    minimizer->FixParameter(8);
    minimizer->FixParameter(9);
    //minimizer->FixParameter(10);
    //minimizer->FixParameter(11);
    minimizer->FixParameter(12);
    minimizer->FixParameter(13);
    minimizer->FixParameter(14);
    minimizer->FixParameter(15);
    //minimizer->FixParameter(16);
    //minimizer->FixParameter(17);
    minimizer->FixParameter(18);
    minimizer->FixParameter(19);
    minimizer->FixParameter(20);
    minimizer->FixParameter(21);


    if(minenergy < 1.4){
        minimizer->FixParameter(16);
        minimizer->FixParameter(17);
        minimizer->FixParameter(8);
        minimizer->FixParameter(9);
        minimizer->SetParameter(24,"rhop rhom II ",0.,0.001,0.,100);
        minimizer->FixParameter(24);
        minimizer->FixParameter(25);
    }
    if(minenergy < 1.4){
        minimizer->FixParameter(22);
        minimizer->FixParameter(23);
    }
    

    minimizer->SetFCN(minuitCMD3);
    cout << "Nexp = " <<Nexp << endl;

    //minimizer->ExecuteCommand("Migrad",0,0);
    double arglist[100];
    arglist[0] = 0;
    // set print level
    //minimizer->ExecuteCommand("SET PRINT",arglist,2);
    arglist[0] = 400000; // number of function calls
    arglist[1] = 0.0001; // tolerance
    //minimizer->ExecuteCommand("Simplex",arglist,2);
    minimizer->ExecuteCommand("Migrad",arglist,2);
    //minimizer->ExecuteCommand("IMProve",arglist,2);
   
    int npars = minimizer->GetNumberTotalParameters();
    double lparforlik[100];
    for(int i = 0; i < npars; i++){
      lparforlik[i] = minimizer->GetParameter(i);
    }
    for(int i = 0; i < N_amplitudes*2.+1; i = i +2){
        resultssss[i] = minimizer->GetParameter(i/2.);
        if((i/2)%2==1 && i/2. > 1)resultssss[i] = correctphase(resultssss[i]);
    }
    for(int i = 1; i < N_amplitudes*2.+2; i = i + 2){
        resultssss[i] = minimizer->GetParError((i-1)/2.);
    }
    
    cout << " ================================================ " << endl;
    cout << " Likelyhood_4pic(resultssss) = " << Likelyhood_4pic(lparforlik) << endl;
    cout << " ================================================ " << endl;
    // Likelyhood_4pic(resultssss)
    if(Likelyhood_4pic(lparforlik) < minimum){
      if(minimizer->GetParError(2)/minimizer->GetParameter(2) > 1)continue;
      results[0] = Likelyhood_4pic(lparforlik);
      minimum = Likelyhood_4pic(lparforlik);
      for(int i = 1; i < N_amplitudes*2.+2; i = i + 1){results[i] = resultssss[i];}
    }
    k++;
    }
    
/*
    Part_ana.fill_sample();Part_ana.fill_sample_c(); 
    for(int i = 0; i < N_amplitudes; i++){Part_ana.change(i,minimizer->GetParameter(i));}
    Part_ana.print();

    TMatrix left(1,50);TMatrix right(50,1);TMatrix res(1,1);
    for(int i = 0; i < 50; i++){left(0,i) = 0; right(i,0) = 0;}
    TMatrix cov(50,50);
    for(int i = 0; i < 50; i++){for(int j = 0; j < 50; j++){cov(i,j) = 0.;}}
    int nper = 0;

    for(int i =0; i < N_amplitudes; i++){
	if(minimizer->GetParError(i) > 0){
		Part_ana.change(i,minimizer->GetParameter(i)+0.001);
		left(0,nper) = Part_ana.cross_section_c()/Part_ana.cross_section();
		Part_ana.change(i,minimizer->GetParameter(i));
		Ratio = Part_ana.cross_section_c()/Part_ana.cross_section();
		left(0,nper) = left(0,nper) - Ratio;
		left(0,nper) = left(0,nper)/0.001;
		nper++;
	}
    }
    for(int i =0; i < nper; i++){right(i,0) = left(0,i);}    

    for(int i = 0; i < nper; i++){
	for(int j = 0; j < nper; j++){
		cov(i,j) = minimizer->GetCovarianceMatrixElement(i,j);
	}
    }

    res = left*cov*right;
    dratio = sqrt(res(0,0));
*/

    return 1;
}






double energy1[1000], energy2[1000];
double energy1_btv[1000], energy2_btv[1000];
int nbins = 0;
int dNevents = 4000;

void cover_data(double enA, double enB, int nevents){
    
    energy1_btv[0] = 0.95;
    double beam, A, B, C1;
    ifstream streamNexp("../4pisel/histograms/data/4pi_data.dat");
    TH1D *henergy = new TH1D("henergy","henergy",10000,0,3);
    while(streamNexp.eof()==0){
        
        streamNexp >> beam;
        streamNexp >> A >> A >> B >> C1 >> A >> A >> B >> C1 >> A >> A >> B >> C1;
        henergy->Fill(beam*2./1000.);
        
    }
    henergy->Draw();
    int current = 0;
    nbins = 0;
    for(int i = 1; i < henergy->GetNbinsX(); i++){
	bool informat = true;
        if(nbins > 0 && fabs(henergy->GetBinCenter(i) - energy2_btv[nbins-1]) > 0.05 && fabs(henergy->GetBinCenter(i) - 1.15) < 0.2 )informat = false;
        if(current < nevents && informat){
            current = current + henergy->GetBinContent(i);
        }
        else{
            current = 0;
            energy2_btv[nbins] = henergy->GetBinCenter(i);
            if(nbins == 0)energy2_btv[nbins]= 0.999;
	    if(nbins > 0)energy1_btv[nbins] = energy2_btv[nbins-1];
            nbins++;

        }	
    }
    
    int nbins_btv = nbins;
    nbins = 0;
    for(int i = 0; i < nbins_btv; i++){

        if(energy2_btv[i] < enA || energy1_btv[i] > enB)continue;
        energy1[nbins] = energy1_btv[i];
        energy2[nbins] = energy2_btv[i];
        nbins++;

    }
    cout << "nbins = " << nbins << endl;
    for(int i = 0; i < nbins; i++){
    
	cout << i << " " << energy1[i] << " " << energy2[i] << endl;	    
    	    
    }
    
}

void script_minim(string nname){

    ifstream stream(nname.c_str());
    double startt; double finishh;
    //cover_data(startt,finishh,nrun);
    ofstream streamr("data/results.dat",ios_base::app);
    ofstream streamm("data/Nevents.dat",ios_base::app);
    ofstream streamratio("data/ratio.dat",ios_base::app);

    while(stream.eof() == 0){
       stream >> startt >> startt >> finishh;
       if(stream.eof() == 1)break;
       minimimlik(startt,finishh);
       if(/*Nexp < 100 || Nphsp < 100 || */Nexp2pi2pi0 < 100 || Nphsp2pi2pi0 < 100){continue;}
       streamr << energy_cm_aver << " ";
       for(int k = 0; k < N_amplitudes*2.; k++){
            if(fabs(results[k]) < 0.000001) results[k] = 0;
            streamr << results[k] << " ";
       }
       streamr << " end" << endl;
       double ratioooo = Nexp2pi2pi0/(Nexp+0.1);
       streamm << energy_cm_aver << " " << startt << " " << finishh <<  " " << Nexp << " " << Nphsp << " " << Nexp2pi2pi0 << " " << Nphsp2pi2pi0 << " " << ratioooo << endl;
       streamratio << energy_cm_aver  << " " << Ratio << " " << dratio << endl;
     
    }
}


void script_fill_cross(){
    double estart = 0.6;
    double eend = 2.;
    double estep = 0.05;
    int nrun = (eend-estart)/estep+1; 
    ofstream streamr2("data/omega_a1pi_f0pi.dat");
    for(int i = 0; i < nrun; i++){
       cout << i << endl;
       //streamr2 << estart+i*estep + estep/2. << " ";
       Part_ana.Energy = estart+i*estep + estep/2.;
       Part_ana.fill0_parameters();
       Part_ana.rhof2 = 1.;

       //Part_ana.fill_parameters();
       Part_ana.print();
       Part_ana.fill_sample();
       streamr2 << Part_ana.cross_section() << " ";
    }

}

void draw(){

    double e[1000];
    double cr[1000];
    ifstream streamr("data/omega_a1pi_f0pi.dat");
    int i = 0;
    while(streamr.eof()==0){

        streamr >> e[i] >> cr[i];
        i++;

    }
    TCanvas *s = new TCanvas();
    TH1F *fr  = s->DrawFrame(0.4,0.,2.07,51.);
    fr->SetXTitle("E_{c.m.}, MeV");
    fr->SetYTitle("#sigma (e^{+}e^{-} #rightarrow 2#pi^{#pm}2#pi^{0}) ,nb");

    TGraphErrors *Cross2   = new TGraphErrors(i-1,e,cr,0,0);
    Cross2->SetMarkerColor(4);
    Cross2->SetMarkerStyle(20);
    Cross2->SetLineColor(4);
    Cross2->SetLineWidth(2.);
    Cross2->Draw("P");
}

double cos_angle(TLorentzVector p1, TLorentzVector p2){
	
	return (p1.E()*p2.E() - p1.Dot(p2))/p1.P()/p2.P();	
	
}

double cos_angle_cm(TLorentzVector p1, TLorentzVector p2){
	TLorentzVector p = p1 + p2;
	TLorentzVector p3 = p1;
	p3.Boost(-p.Px()/p.E(),-p.Py()/p.E(),-p.Pz()/p.E());
	return cos_angle(p3,p);	
}

double Cmd3Generator2pi2pi0_ke::generate(int niter = 100000)
{
  find_majoranta0(Energy);
  
  TRandom1 rand;
  double xnew = rand.Rndm();

  TGenPhaseSpace event;
  TLorentzVector* P1 = new TLorentzVector();
  TLorentzVector* P2 = new TLorentzVector();
  TLorentzVector* P3 = new TLorentzVector();
  TLorentzVector* P4 = new TLorentzVector();
  TLorentzVector* Ph = new TLorentzVector();
  double masses[4] = {mpiz, mpiz, Mpi, Mpi};

  TH1D *minv = new TH1D("minv","minv",200,0,2);
  TH1D *minvcn = new TH1D("minvcn","minvcn",200,0,2);
  TH1D *minvcc = new TH1D("minvcc","minvcc",200,0,2);
  TH1D *minvnn = new TH1D("minvnn","minvnn",200,0,2);
  TH1D *hcos_pippimmc = new TH1D("hcos_pippimmc","hcos_pippimmc",100,-1,1.);
  hcos_pippimmc->SetXTitle("cos (#pi#pm)_{ +- c.m.}");
  hcos_pippimmc->SetYTitle("yields");

  TH1D *hEph = new TH1D("hEph","hEph",1000,-0.05,2);
  TH1D *hEth = new TH1D("hEth","hEth",1000,-0.05,3.2);

  double cross = 0;
  ofstream matrixstream("data/matrix_toy.dat");
  for(int iii = 0; iii < niter; iii++){
  	if(iii%100 == 0)cout << iii << " "<< endl;
  	int proba = 0;
  	double cr123234546;
  	TLorentzVector P;
  	while(proba < 5000){
      		double ephoton = 0.;//generate_photon_energy(Energy);
      		double costhetaphoton = generate_photon_angle(Energy, 2.*ephoton/Energy);
      		double sinthetaphoton = sqrt(1. - costhetaphoton*costhetaphoton);
      		Ph = new TLorentzVector(ephoton*sinthetaphoton*cos(xnew*2.*3.14), ephoton*sinthetaphoton*sin(xnew*2.*3.14), ephoton*costhetaphoton,ephoton);
      		P.SetPxPyPzE(-Ph->X(),-Ph->Y(),-Ph->Z(),Energy - ephoton);
      		event.SetDecay(P, 4, masses);
      		cr123234546 = event.Generate();
      		P1 = event.GetDecay(0);
      		P2 = event.GetDecay(1);
      		P3 = event.GetDecay(2);
      		P4 = event.GetDecay(3);
		TLorentzVectorC Pc(P.X(),P.Y(),P.Z(),P.E());
		TLorentzVectorC P1c(P1->X(),P1->Y(),P1->Z(),P1->E());
		TLorentzVectorC P2c(P2->X(),P2->Y(),P2->Z(),P2->E());
		TLorentzVectorC P3c(P3->X(),P3->Y(),P3->Z(),P3->E());
		TLorentzVectorC P4c(P4->X(),P4->Y(),P4->Z(),P4->E());
		TLorentzVectorC btvv = matrix(Pc,P1c,P2c,P3c,P4c);
		btvv = btvv + matrix(Pc,P3c,P2c,P1c,P4c);
		btvv = btvv + matrix(Pc,P1c,P4c,P3c,P2c);
		btvv = btvv + matrix(Pc,P3c,P4c,P1c,P2c);
                cout << btvv.X() << " ";
                cout << btvv.Y() << " ";
                cout << btvv.Z() << " ";
                cout << btvv.E() << endl;

      		xnew = rand.Rndm();
      		if(proba > 40000.) cout << proba << " HEEEEEEEEEEEEEEEEEELP! ";
      		if(matrix_squared1(P,*P1,*P2,*P3,*P4)*cr123234546 < xnew*majoranta){proba++;continue;}
      		break;
	}
  	TLorentzVectorC zero(P.X(),P.Y(),P.Z(),P.E());
  	TLorentzVectorC one(P1->X(),P1->Y(),P1->Z(),P1->E());
  	TLorentzVectorC two(P2->X(),P2->Y(),P2->Z(),P2->E());
  	TLorentzVectorC three(P3->X(),P3->Y(),P3->Z(),P3->E());
  	TLorentzVectorC four(P4->X(),P4->Y(),P4->Z(),P4->E());
/*
  matrix << Energy*1000. << " " << one.X().real() << " " << H_omega_pi0(zero,one,two,three,four).X().imag() << " ";
  matrix << H_omega_pi0(zero,one,two,three,four).Y().real() << " " << H_omega_pi0(zero,one,two,three,four).Y().imag() << " ";
  // matrix << H_omega_pi0(zero,one,two,three,four).Z().real() << " "<< H_omega_pi0(zero,one,two,three,four).Z().imag() << " ";
  // matrix << H_omega_pi0(zero,one,two,three,four).E().real() << " " << H_omega_pi0(zero,one,two,three,four).E().imag() << " ";
  matrix << H_a1_rhopi_pi(zero,one,two,three,four).X().real() << " " << H_a1_rhopi_pi(zero,one,two,three,four).X().imag() << " ";
  matrix << H_a1_rhopi_pi(zero,one,two,three,four).Y().real() << " " << H_a1_rhopi_pi(zero,one,two,three,four).Y().imag() << " ";
  // matrix << H_a1_rhopi_pi(zero,one,two,three,four).Z().real() << " " << H_a1_rhopi_pi(zero,one,two,three,four).Z().imag() << " ";
  // matrix << H_a1_rhopi_pi(zero,one,two,three,four).E().real() << " " << H_a1_rhopi_pi(zero,one,two,three,four).E().imag() << " ";
   matrix << H_a1_sigmapi_pi(zero,one,two,three,four).X().real() << " " << H_a1_sigmapi_pi(zero,one,two,three,four).X().imag() << " ";
   matrix << H_a1_sigmapi_pi(zero,one,two,three,four).Y().real() << " " << H_a1_sigmapi_pi(zero,one,two,three,four).Y().imag() << " ";
  // matrix << H_a1_sigmapi_pi(zero,one,two,three,four).Z().real() << " " << H_a1_sigmapi_pi(zero,one,two,three,four).Z().imag() << " ";
  // matrix << H_a1_sigmapi_pi(zero,one,two,three,four).E().real() << " " << H_a1_sigmapi_pi(zero,one,two,three,four).E().imag() << endl;
  matrix << endl;*/
  	minv->Fill((*P1+*P3+*P4+*P2).M());
 	//minv->Fill((*P2+*P3+*P4).M());
  	//if(fabs((*P1+*P3+*P4).M()-0.787)<0.02){continue;}
  	//if(fabs((*P2+*P3+*P4).M()-0.787)<0.02){continue;}
  	minvcn->Fill((*P1+*P3).M());minvcn->Fill((*P1+*P4).M());
  	minvcn->Fill((*P2+*P3).M());minvcn->Fill((*P2+*P4).M());
  	minvcc->Fill((*P4+*P3).M());
  	minvnn->Fill((*P1+*P2).M());
	hcos_pippimmc->Fill(cos_angle_cm(*P3,*P4));
	//hcos_pippimmc->Fill(cos_angle_cm(*P4,*P3));
        hEph->Fill(Ph->E());
        if(Ph->E() > 0.001)hEth->Fill(Ph->Theta());
  	cross = cross + cr123234546;
  }
  cross = cross/niter;
  TCanvas *s = new TCanvas();
  s->Divide(2,2);
  s->cd(1);
  minv->Draw();
  s->cd(2);
  minvcn->Draw();
  s->cd(3);
  minvcc->Draw();
  s->cd(4);
  minvnn->Draw();

  TCanvas *s1 = new TCanvas();
  s1->Divide(2,2);
  s1->cd(1);
  hcos_pippimmc->Draw();
  s1->cd(2);
  hEph->Draw();
  s1->cd(3);
  hEth->Draw();
  return cross;
}




double goodness(double enA, double enB){

  double goodn = 0.;double nmc = 0;
  Part_ana.readd(enA,enB,"2pi2pi0_v1/");

  double en = (enA+enB)/2.;
  Part_ana.Energy  = en;
  Part_ana.fill_parameters();
  double parrs[]={0.,Part_ana.omegapi0_,Part_ana.a1pi_rhopi_,Part_ana.a1pi_rhopi_phase,Part_ana.a1pi_sigmapi_,Part_ana.a1pi_sigmapi_phase, Part_ana.rhof0_, Part_ana.rhof0phase_, Part_ana.rhoprhom_, Part_ana.rhoprhom_phase, Part_ana.rhosigma_, Part_ana.rhosigmaphase_, Part_ana.phsp_, Part_ana.phspphase_, Part_ana.a2pi_, Part_ana.a2pi_ph, Part_ana.h1pi, Part_ana.h1pi_ph, Part_ana.a1pi_rhopi_II, Part_ana.a1pi_rhopi_phase_II, Part_ana.rhof0_inter, Part_ana.rhof0phase_inter, Part_ana.rhof2, Part_ana.rhof2_phase, Part_ana.rhoprhom_II, Part_ana.rhoprhom_phase_II};

  Likelyhood_4pic(parrs);
  double likexp = 0,likmc = 0;
  for(int i = 0; i < Nexp2pi2pi0; i++){
    likexp = likexp + Likexp[i];
  }

  TH1D *hmc = new TH1D("hmc","hmc",100,likexp-500,likexp+500);
  hmc->SetYTitle("yields");
  hmc->SetXTitle("#Sigma^{Nexp}_{i}(log(|M|_{i}^{2}))");
  for(int i = 0; i < Nphsp2pi2pi0; i++){
    if(i%(int)Nexp2pi2pi0==1 && likmc > 0 && i > 0){
      hmc->Fill(likmc);
      if(likexp > likmc)goodn++;
      nmc++;
      likmc = 0.;
    }  
    likmc = likmc + Likmc[i];
  }

  hmc->SetTitle(Form("E_{c.m.} = %g -- %g GeV, goodness = %g",enA,enB,((int)(goodn/nmc*100.))/100.));
  hmc->Draw();
  TArrow *f = new TArrow(likexp,2,likexp,0,0.05,">");
  f->SetLineColor(2);
  f->Draw();
  
  cout << "likexp = " << likexp << endl;
				   
  return goodn/nmc;
}


  void drawgood(){
    double good[100], en[100],den[100];
    //ifstream stream("run_minimizer/all");
    double ena,enb;
    ena = 0.95;
    int i = 0;
    while(ena < 2.){
      //stream >> ena >> ena >> enb;
      // if(stream.eof()==1)break;
      //cout << "GGGGGGGGGGOOOOOOOOOOOOOOOOOOOOOOOOOODDDDDDDDDDDD " << goodness(ena,enb) << endl;
      ena = 0.95 + 0.08*i;
      enb = ena + 0.08;
      en[i] = 0.5*ena + 0.5*enb;
      den[i] = 0.5*(enb - ena);
      good[i] = goodness(ena,enb);
      i++;
    }

    TCanvas *s = new TCanvas();
    TH1F *frd  = s->DrawFrame(en[0]-1,-10.,en[i-1]+1,10.);
    frd->SetXTitle("");
    frd->SetYTitle("");
    
   TGraphErrors *Cross  = new TGraphErrors(i-2,en,good,den,0);
   Cross->SetMarkerColor(2);
   Cross->SetMarkerStyle(20);
   Cross->SetLineColor(2);
   Cross->SetLineWidth(2.);
   Cross->SetTitle("");
   Cross->Draw("P");

  }

