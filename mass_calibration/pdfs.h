#include <iostream>
#include <TLegend.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1D.h"
#include "TH2D.h"
#include <THStack.h>
#include "TProfile.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFractionFitter.h"
#include <string>
#include <vector>
#include <math.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TString.h>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TEfficiency.h"
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <iostream>
#include <valarray>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooBernstein.h>
#include <RooBreitWigner.h>
#include <RooVoigtian.h>
#include <RooFFTConvPdf.h>
#include <RooNumConvPdf.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooMsgService.h> 
#include "DSCB.h"

using namespace std;

class reso_pdfs {
	public:

		RooWorkspace* w;
		//// common params
		RooRealVar* alpha1;
		RooRealVar* alpha2;
		RooRealVar* n1;
		//RooRealVar* n2; NOT CURRENTLY USED IN FIT
		RooRealVar* frac_gau;
		RooRealVar* frac_gau2;
		RooRealVar* gau_reso_scale;
		RooRealVar* gau_reso_scale2;
		RooRealVar* gau_reso_scale3;
		//////

		~reso_pdfs(){
			w->Delete();
		}

		reso_pdfs(){
			w = new RooWorkspace("dpworkspace", "");
			
			//// common paramaters
			alpha1 = new RooRealVar("alpha1"  	, "alpha1 " 	, 1.,0.,2.5);
			alpha2 = new RooRealVar("alpha2"  	, "alpha2 " 	, 1.,0.,2.5);
			n1 = new RooRealVar("n1"  	, "n1 " 	, 1, 0.01, 5);//0.5,5);
			//n2 = RooRealVar("n2"  	, "n2 " 	, 1,0.5,5); N2 IS NOT CURRENTLY USED IN THE FIT
			frac_gau = new RooRealVar("frac_gau", "frac_gau", 0.5,0,1);
			frac_gau2 = new RooRealVar("frac_gau2", "frac_gau2", 0.5,0,1);
			gau_reso_scale = new RooRealVar ("gau_reso_scale", "gau_reso_scale", 0.5, 0.1, 3.5);//1.5);
			gau_reso_scale2 = new RooRealVar ("gau_reso_scale2", "gau_reso_scale2", 0.5, 0.1, 3.5);//1.5);
			gau_reso_scale3 = new RooRealVar ("gau_reso_scale3", "gau_reso_scale3", 0.5, 0.1, 3.5);
			//////
			

			///// define all pdfs

			//ETA MESON//

			//Define observable and parameters
			RooRealVar* m2mu_eta = new RooRealVar("m2mu_eta", "m2mu_eta", 0, 15);
			RooRealVar* w_eta = new RooRealVar("w_eta", "w_eta", 0.00000131, 0.000001, 0.000002); //Only used in generating breit-wigner for convolving
			
			//Create signal model
			RooAddPdf* signalModel_eta = make_signal_model("eta", m2mu_eta, 0.55, 0.535, 0.565, w_eta);

			RooRealVar sig_fraction_eta("sig_fraction_eta", "",0.2,0.02,0.8);

			//Create background model
			RooRealVar par1_eta("par1_eta", "par1_eta", 0.2, 0, 10);
			RooRealVar par2_eta("par2_eta", "par2_eta", 1.5, 0, 10);
			RooRealVar par3_eta("par3_eta", "par3_eta", 2.0, 0, 10);
			RooRealVar par4_eta("par4_eta", "par4_eta", 2.0, 0, 10);
			RooRealVar par5_eta("par5_eta", "par5_eta", 2.0, 0, 10);
			RooRealVar par6_eta("par6_eta", "par6_eta", 2.0, 0, 10);
			RooRealVar par7_eta("par7_eta", "par7_eta", 0.2, 0, 10);
			RooRealVar par8_eta("par8_eta", "par8_eta", 1.5, 0, 10);
			RooRealVar par9_eta("par9_eta", "par9_eta", 2.0, 0, 10);
			RooArgList parlist_eta(par1_eta, par2_eta, par3_eta, par4_eta, par5_eta);//, par6_eta, par7_eta, par8_eta);//, par9_eta);
			RooBernstein bkgModel_eta("bkgModel_eta", "bkgModel_eta", *m2mu_eta, parlist_eta);

			//Combine signal and bkg
			RooAddPdf model_eta("model_eta", "model_eta", RooArgList(*signalModel_eta, bkgModel_eta), RooArgList(sig_fraction_eta));

			w->import(model_eta);

			//OMEGA+PHI MESONS//
			
			//Define observable and parameters
			RooRealVar* m2mu_omegaphi = new RooRealVar("m2mu_omegaphi", "m2mu_omegaphi", 0,15);
			RooRealVar* w_omega = new RooRealVar("w_omega", "w_omega", 0.0005, 0.0005, 0.0005);//0.00849, 0.00849, 0.00849);
			RooRealVar* w_phi = new RooRealVar("w_phi", "w_phi", 0.004249, 0.0041, 0.0043); 
			RooRealVar* m_omega = new RooRealVar("m_omega", "m_omega", 0); //0 mean for convolving

			//Create signal models
			RooAddPdf* signalModel_omega = make_signal_model("omega", m2mu_omegaphi, 0.782, 0.77, 0.795, w_omega);
			RooAddPdf* signalModel_phi = make_signal_model("phi", m2mu_omegaphi, 1.02, 0.95, 1.1, w_phi);

			RooRealVar sig_fraction_omega("sig_fraction_omega", "",0.2,0.02,0.8);
			RooRealVar sig_fraction_phi("sig_fraction_phi", "",0.2,0.02,0.8);

			//For Manual Convolving of Omega
			RooBreitWigner* omega_BW = new RooBreitWigner("omega_bw", "omega_bw", *m2mu_omegaphi, *m_omega, *w_omega);
			//m2mu_omegaphi->setBins(5500, "cache"); //Changes from default binning, doesn't finish running without this statement
			//RooFFTConvPdf* omega_conv = new RooFFTConvPdf("omega_conv", "omega_conv", *m2mu_omegaphi, *omega_BW, *signalModel_omega);//Fourier Transform convolving, loops without changing binning
			RooNumConvPdf* omega_conv = new RooNumConvPdf("omega_conv", "omega_conv", *m2mu_omegaphi, *omega_BW, *signalModel_omega);//Brute force convolving, takes really long time
			

			//Creates background models
			RooRealVar par1_omegaphi("par1_omegaphi", "par1_omegaphi", 0.2, 0, 15);
			RooRealVar par2_omegaphi("par2_omegaphi", "par2_omegaphi", 1.5, 0, 15);
			RooRealVar par3_omegaphi("par3_omegaphi", "par3_omegaphi", 2.0, 0, 15);
			RooRealVar par4_omegaphi("par4_omegaphi", "par4_omegaphi", 2.0, 0, 15);
			RooRealVar par5_omegaphi("par5_omegaphi", "par5_omegaphi", 2.0, 0, 15);
			RooRealVar par6_omegaphi("par6_omegaphi", "par6_omegaphi", 2.0, 0, 15);
			RooRealVar par7_omegaphi("par7_omegaphi", "par7_omegaphi", 0.2, 0, 15);
			RooRealVar par8_omegaphi("par8_omegaphi", "par8_omegaphi", 1.5, 0, 15);
			RooRealVar par9_omegaphi("par9_omegaphi", "par9_omegaphi", 2.0, 0, 15);
			RooArgList parlist_omegaphi(par1_omegaphi, par2_omegaphi, par3_omegaphi, par4_omegaphi, par5_omegaphi);//, par6_omegaphi, par7_omegaphi, par8_omegaphi);//, par9_omegaphi);
			RooBernstein bkgModel_omegaphi("bkgModel_omegaphi", "bkgModel_omegaphi", *m2mu_omegaphi, parlist_omegaphi);

			//Combine sig and bkg
			//RooAddPdf model_omegaphi("model_omegaphi", "model_omegaphi", RooArgList(*signalModel_omega, *signalModel_phi, bkgModel_omegaphi), RooArgList(sig_fraction_omega, sig_fraction_phi));
			RooAddPdf model_omegaphi("model_omegaphi", "model_omegaphi", RooArgList(*omega_conv, *signalModel_phi, bkgModel_omegaphi), RooArgList(sig_fraction_omega, sig_fraction_phi));//Convolved omega

			w->import(model_omegaphi);
			
			//J/Psi+Psi2s MESONS//

			//Define obervable and parameters
			RooRealVar* m2mu_jpsipsi = new RooRealVar("m2mu_jpsipsi", "m2mu_jpsipsi", 0, 15);
			RooRealVar* w_jpsi = new RooRealVar("w_jpsi", "w_jpsi", 0.0000929, 0.000091, 0.000094);//only used in convolving
			RooRealVar* w_psi2s = new RooRealVar("w_psi2s", "w_psi2s", 0.000294, 0.000293, 0.000294);//only used in convolving
			
			//Creates signal models
			RooAddPdf* signalModel_jpsi = make_signal_model("jpsi", m2mu_jpsipsi, 3.1, 3.0, 3.2, w_jpsi);
			RooAddPdf* signalModel_psi2s = make_signal_model("psi2s", m2mu_jpsipsi, 3.65, 3.6, 3.8, w_psi2s);

			RooRealVar sig_fraction_jpsi("sig_fraction_jpsi", "",0.3,0.1,0.9);
			RooRealVar sig_fraction_psi2s("sig_fraction_psi2s", "",0.1,0.02,0.5);

			//Creates background models
			RooRealVar par1_jpsipsi("par1_jpsipsi", "par1_jpsipsi", 0.2, 0, 10);
			RooRealVar par2_jpsipsi("par2_jpsipsi", "par2_jpsipsi", 1.5, 0, 10);
			RooRealVar par3_jpsipsi("par3_jpsipsi", "par3_jpsipsi", 2.0, 0, 10);
			RooRealVar par4_jpsipsi("par4_jpsipsi", "par4_jpsipsi", 2.0, 0, 10);
			RooRealVar par5_jpsipsi("par5_jpsipsi", "par5_jpsipsi", 2.0, 0, 10);
			RooRealVar par6_jpsipsi("par6_jpsipsi", "par6_jpsipsi", 2.0, 0, 10);
			RooRealVar par7_jpsipsi("par7_jpsipsi", "par7_jpsipsi", 0.2, 0, 10);
			RooRealVar par8_jpsipsi("par8_jpsipsi", "par8_jpsipsi", 1.5, 0, 10);
			RooRealVar par9_jpsipsi("par9_jpsipsi", "par9_jpsipsi", 2.0, 0, 10);
			RooArgList parlist_jpsipsi(par1_jpsipsi, par2_jpsipsi, par3_jpsipsi, par4_jpsipsi, par5_jpsipsi);//, par6_jpsipsi, par7_jpsipsi, par8_jpsipsi);//, par9_jpsipsi);
			RooBernstein bkgModel_jpsipsi("bkgModel_jpsipsi", "bkgModel_jpsipsi", *m2mu_jpsipsi, parlist_jpsipsi);
			
			//Combine sig and bkg
			RooAddPdf model_jpsipsi("model_jpsipsi", "model_jpsipsi", RooArgList(*signalModel_jpsi, *signalModel_psi2s, bkgModel_jpsipsi), RooArgList(sig_fraction_jpsi,sig_fraction_psi2s));

			w->import(model_jpsipsi);
			
			//UPSILON1S+2S+3S MESONS//

			//Define observable and parameters
			RooRealVar* m2mu_upsilon = new RooRealVar("m2mu_upsilon", "m2mu_upsilon", 0, 15);
			RooRealVar* w_upsilon1s = new RooRealVar("w_upsilon1s", "w_upsilon1s", 0.0, 0.0, 0.1);
			RooRealVar* w_upsilon2s = new RooRealVar("w_upsilon2s", "w_upsilon2s", 0.0, 0.0, 0.1);
			RooRealVar* w_upsilon3s = new RooRealVar("w_upsilon3s", "w_upsilon3s", 0.0, 0.0, 0.1);
			
			//Creates signal models
			RooAddPdf* signalModel_upsilon1s = make_signal_model("upsilon1s", m2mu_upsilon, 9.46, 9.36, 9.56, w_upsilon1s);
			RooAddPdf* signalModel_upsilon2s = make_signal_model("upsilon2s", m2mu_upsilon, 10.05, 9.95, 10.1, w_upsilon2s);
			RooAddPdf* signalModel_upsilon3s = make_signal_model("upsilon3s", m2mu_upsilon, 10.3, 10.25, 10.4, w_upsilon3s);

			RooRealVar sig_fraction_upsilon1s("sig_fraction_upsilon1s", "",0.2,0.02,0.8);
			RooRealVar sig_fraction_upsilon2s("sig_fraction_upsilon2s", "",0.2,0.02,0.8);
			RooRealVar sig_fraction_upsilon3s("sig_fraction_upsilon3s", "",0.2,0.02,0.8);

			//Creates background models
			RooRealVar par1_upsilon("par1_upsilon", "par1_upsilon", 0.2, 0, 10);
			RooRealVar par2_upsilon("par2_upsilon", "par2_upsilon", 1.5, 0, 10);
			RooRealVar par3_upsilon("par3_upsilon", "par3_upsilon", 2.0, 0, 10);
			RooRealVar par4_upsilon("par4_upsilon", "par4_upsilon", 2.0, 0, 10);
			RooRealVar par5_upsilon("par5_upsilon", "par5_upsilon", 2.0, 0, 10);
			RooRealVar par6_upsilon("par6_upsilon", "par6_upsilon", 2.0, 0, 10);
			RooRealVar par7_upsilon("par7_upsilon", "par7_upsilon", 0.2, 0, 10);
			RooRealVar par8_upsilon("par8_upsilon", "par8_upsilon", 1.5, 0, 10);
			RooRealVar par9_upsilon("par9_upsilon", "par9_upsilon", 2.0, 0, 10);
			RooArgList parlist_upsilon(par1_upsilon, par2_upsilon, par3_upsilon, par4_upsilon, par5_upsilon);//, par6_upsilon, par7_upsilon, par8_upsilon);//, par9_upsilon);
			RooBernstein bkgModel_upsilon("bkgModel_upsilon", "bkgModel_upsilon", *m2mu_upsilon, parlist_upsilon);

			//Combine sig and bkg
			RooAddPdf model_upsilon("model_upsilon", "model_upsilon", RooArgList(*signalModel_upsilon1s, *signalModel_upsilon2s, *signalModel_upsilon3s, bkgModel_upsilon), RooArgList(sig_fraction_upsilon1s, sig_fraction_upsilon2s, sig_fraction_upsilon3s));

			w->import(model_upsilon);

		};
		/*//SIGNALS MADE FROM GAUSSIAN + CB
		RooAddPdf* make_signal_model(string name, RooRealVar* mass, double m_init, double m_min, double m_max, RooRealVar* width){
			RooRealVar* M = new RooRealVar(("M_"+name).c_str(), ("M_"+name).c_str(), m_init, m_min, m_max);

			RooRealVar* res_rel = new RooRealVar(("res_rel_"+name).c_str(), ("res_rel_"+name).c_str(), 0.01, 0.005, 0.1);

			RooFormulaVar* res_CB = new RooFormulaVar(("res_CB_"+name).c_str(), ("M_"+name+"*res_rel_"+name).c_str() , RooArgList(*M, *res_rel));
			RooDoubleCB* signalModel_CB = new RooDoubleCB(("signalModel_CB_"+name).c_str(), ("signalModel_CB_"+name).c_str(), *mass, *M,*res_CB,*alpha1,*n1,*alpha2,*n1);

			RooFormulaVar* res_gau = new RooFormulaVar(("res_gau_"+name).c_str(), ("gau_reso_scale*M_"+name+"*res_rel_"+name).c_str(), RooArgList(*gau_reso_scale, *M,*res_rel));
			RooGaussian* signalModel_gau = new RooGaussian(("signalModel_gau_"+name).c_str(), ("signalModel_gau_"+name).c_str(), *mass, *M,*res_gau);

			RooAddPdf* signalModel = new RooAddPdf(("signalModel_"+name).c_str(), ("signalModel_"+name).c_str(), RooArgList(*signalModel_CB, *signalModel_gau), RooArgList(*frac_gau));
			return signalModel;
		};*/

		//SIGNALS MADE FROM 3 GAUSSIANS
		RooAddPdf* make_signal_model(string name, RooRealVar* mass, double m_init, double m_min, double m_max, RooRealVar* width){
			RooRealVar* M1 = new RooRealVar(("M1_"+name).c_str(), ("M1_"+name).c_str(), m_init-0.01, m_min-0.05, m_max+0.05);
			RooRealVar* M2 = new RooRealVar(("M2_"+name).c_str(), ("M2_"+name).c_str(), m_init+0.01, m_min-0.05, m_max+0.05);
			RooRealVar* M3 = new RooRealVar(("M3_"+name).c_str(), ("M3_"+name).c_str(), m_init, m_min-0.05, m_max+0.05);

			RooRealVar* res_rel1 = new RooRealVar(("res_rel1_"+name).c_str(), ("res_rel1_"+name).c_str(), 0.01, 0.005, 0.2);
			RooRealVar* res_rel2 = new RooRealVar(("res_rel2_"+name).c_str(), ("res_rel2_"+name).c_str(), 0.01, 0.005, 0.2);
			RooRealVar* res_rel3 = new RooRealVar(("res_rel3_"+name).c_str(), ("res_rel3_"+name).c_str(), 0.01, 0.005, 0.2);

			RooFormulaVar* res_gau1 = new RooFormulaVar(("res_gau1_"+name).c_str(), ("gau_reso_scale*M1_"+name+"*res_rel1_"+name).c_str(), RooArgList(*gau_reso_scale, *M1,*res_rel1));
			RooGaussian* signalModel_gau1 = new RooGaussian(("signalModel_gau1_"+name).c_str(), ("signalModel_gau1_"+name).c_str(), *mass, *M1,*res_gau1);
		
			RooFormulaVar* res_gau2 = new RooFormulaVar(("res_gau2_"+name).c_str(), ("gau_reso_scale2*M2_"+name+"*res_rel2_"+name).c_str(), RooArgList(*gau_reso_scale2, *M2,*res_rel2));
			RooGaussian* signalModel_gau2 = new RooGaussian(("signalModel_gau2_"+name).c_str(), ("signalModel_gau2_"+name).c_str(), *mass, *M2,*res_gau2);

			RooFormulaVar* res_gau3 = new RooFormulaVar(("res_gau3_"+name).c_str(), ("gau_reso_scale3*M2_"+name+"*res_rel3_"+name).c_str(), RooArgList(*gau_reso_scale3, *M3,*res_rel3));
			RooGaussian* signalModel_gau3 = new RooGaussian(("signalModel_gau3_"+name).c_str(), ("signalModel_gau3_"+name).c_str(), *mass, *M3,*res_gau3);

			RooAddPdf* signalModel = new RooAddPdf(("signalModel_"+name).c_str(), ("signalModel_"+name).c_str(), RooArgList(*signalModel_gau1, *signalModel_gau2, *signalModel_gau3), RooArgList(*frac_gau));
			//RooAddPdf* signalModel2 = new RooAddPdf(("signalModel2_"+name).c_str(), ("signalModel2_"+name).c_str(), RooArgList(, *signalModel), RooArgList(*frac_gau2));
			return signalModel;
		};
		
		/*//SIGNAL MADE FROM VOIGTIAN + CB
		RooAddPdf* make_signal_model(string name, RooRealVar* mass, double m_init, double m_min, double m_max, RooRealVar* width){
			RooRealVar* M = new RooRealVar(("M_"+name).c_str()      , ("M_"+name).c_str()          , m_init, m_min, m_max);
			RooRealVar* res_rel = new RooRealVar(("res_rel_"+name).c_str(), ("res_rel_"+name).c_str(), 0.01, 0.005, 0.1);
			RooFormulaVar* res_CB = new RooFormulaVar(("res_CB_"+name).c_str(), ("M_"+name+"*res_rel_"+name).c_str() , RooArgList(*M, *res_rel));
			RooDoubleCB* signalModel_CB = new RooDoubleCB(("signalModel_CB_"+name).c_str(), ("signalModel_CB_"+name).c_str(), *mass, *M,*res_CB,*alpha1,*n1,*alpha2,*n1);
			RooFormulaVar* res_gau = new RooFormulaVar(("res_gau_"+name).c_str()  	, ("gau_reso_scale*M_"+name+"*res_rel_"+name).c_str() , RooArgList(*gau_reso_scale, *M,*res_rel));
			//RooGaussian* signalModel_gau = new RooGaussian(("signalModel_gau_"+name).c_str(), ("signalModel_gau_"+name).c_str(), *mass, *M,*res_gau);
			RooVoigtian* signalModel_voi = new RooVoigtian(("signalModel_voi_"+name).c_str(), ("signalModel_voi_"+name).c_str(), *mass, *M, *width, *res_gau);
			RooAddPdf* signalModel = new RooAddPdf(("signalModel_"+name).c_str(), ("signalModel_"+name).c_str(), RooArgList(*signalModel_CB, *signalModel_voi), RooArgList(*frac_gau));
			return signalModel;
		};
		*/
		/*//SIGNAL MADE FROM CB+GAUS CONVOLVED WITH BW, ENDS UP IN LOOP
		RooAddPdf* make_signal_model(string name, RooRealVar* mass, double m_init, double m_min, double m_max, RooRealVar* width){
			RooRealVar* M = new RooRealVar(("M_"+name).c_str()      , ("M_"+name).c_str()          , m_init, m_min, m_max);
			RooRealVar* res_rel = new RooRealVar(("res_rel_"+name).c_str(), ("res_rel_"+name).c_str(), 0.01, 0.005, 0.1);
			RooFormulaVar* res_CB = new RooFormulaVar(("res_CB_"+name).c_str(), ("M_"+name+"*res_rel_"+name).c_str() , RooArgList(*M, *res_rel));
			RooDoubleCB* signalModel_CB = new RooDoubleCB(("signalModel_CB_"+name).c_str(), ("signalModel_CB_"+name).c_str(), *mass, *M,*res_CB,*alpha1,*n1,*alpha2,*n1);
			RooFormulaVar* res_gau = new RooFormulaVar(("res_gau_"+name).c_str()  	, ("gau_reso_scale*M_"+name+"*res_rel_"+name).c_str() , RooArgList(*gau_reso_scale, *M,*res_rel));
			RooGaussian* signalModel_gau = new RooGaussian(("signalModel_gau_"+name).c_str(), ("signalModel_gau_"+name).c_str(), *mass, *M,*res_gau);
			//RooVoigtian* signalModel_voi = new RooVoigtian(("signalModel_voi_"+name).c_str(), ("signalModel_voi_"+name).c_str(), *mass, *M, *width, *res_gau);
			RooAddPdf* signalModel_CB_gau = new RooAddPdf(("signalModel_"+name).c_str(), ("signalModel_"+name).c_str(), RooArgList(*signalModel_CB, *signalModel_gau), RooArgList(*frac_gau));
			RooBreitWigner* signalModel_BW = new RooBreitWigner(("signalModel_BW_"+name).c_str(), ("signalModel_BW_"+name).c_str(), *mass, *M, *width);
			RooFFTConvPdf* signalModel = new RooFFTConvPdf(("signalModel_conv_"+name).c_str(), ("signalModel_conv_"+name).c_str(), *mass, *signalModel_BW, *signalModel_CB_gau);
			return signalModel_CB_gau;
		};*/
};


