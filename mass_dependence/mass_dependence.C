#include "TFileCollection.h"
#include "TChain.h"
#include "TFile.h"
#include <TTreeReader.h>
#include "TH1D.h"
#include "TH2D.h"
#include <TTreeReaderValue.h>
#include "TLorentzVector.h"
#include <iostream>
#include <TLegend.h>
#include <TFrame.h>
#include <sstream>
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "TH1F.h"
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
#include <cctype>
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
#include <RooHist.h>
#include <RooPolynomial.h>
#include <RooBernstein.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooMsgService.h> 
#include <RooAddPdf.h> 
#include <RooAddition.h> 
#include <RooMinuit.h> 
#include <RooFitResult.h> 

void mass_dependence(string treepath = "/mnt/hadoop/scratch/gandreas/ScoutingRunC/ScoutingCaloMuon/crab_20181212_141702/181212_131713/0000/scout_1.root", string reso = "upsilon1"){

	///////////////////////////////
	/////VAR AND MASS ANALYSIS/////
	///////////////////////////////

	//JPSI MASS FOR SCALING
	double real_mass;
	if (reso == "phi") real_mass = 1.019461;
	if (reso == "jpsi") real_mass = 3.0969;
	if (reso == "upsilon1") real_mass = 9.46030;

	///FOR SIDEBAND SUBTRACTED PLOT
	//DIMUON
	//PT
	int bins = 13;
	TH1F* sigpt = new TH1F("sigpt", "signal region", bins, 0., 30.);
	TH1F* sbpt = new TH1F("sbpt", "sideband region", bins, 0., 30.);
	TH1F* sigsubpt;

	//PHI
	TH1F* sigphi = new TH1F("sigphi", "signal region", bins, -3.15, 3.15);
	TH1F* sbphi = new TH1F("sbphi", "sideband region", bins, -3.15, 3.15);
	TH1F* sigsubphi;

	//ETA
	TH1F* sigeta = new TH1F("sigeta", "signal region", bins, -3.15, 3.15);
	TH1F* sbeta = new TH1F("sbeta", "sideband region", bins, -3.15, 3.15);
	TH1F* sigsubeta;

	//MUON 1
	//PT
	TH1F* m1sigpt = new TH1F("m1sigpt", "signal region", bins, 0., 30.);
	TH1F* m1sbpt = new TH1F("m1sbpt", "sideband region", bins, 0., 30.);
	TH1F* m1sigsubpt;

	//PHI
	TH1F* m1sigphi = new TH1F("m1sigphi", "signal region", bins, -3.15, 3.15);
	TH1F* m1sbphi = new TH1F("m1sbphi", "sideband region", bins, -3.15, 3.15);
	TH1F* m1sigsubphi;

	//ETA
	TH1F* m1sigeta = new TH1F("m1sigeta", "signal region", bins, -3.15, 3.15);
	TH1F* m1sbeta = new TH1F("m1sbeta", "sideband region", bins, -3.15, 3.15);
	TH1F* m1sigsubeta;

	//MUON 2
	//PT
	TH1F* m2sigpt = new TH1F("m2sigpt", "signal region", bins, 0., 30.);
	TH1F* m2sbpt = new TH1F("m2sbpt", "sideband region", bins, 0., 30.);
	TH1F* m2sigsubpt;

	//PHI
	TH1F* m2sigphi = new TH1F("m2sigphi", "signal region", bins, -3.15, 3.15);
	TH1F* m2sbphi = new TH1F("m2sbphi", "sideband region", bins, -3.15, 3.15);
	TH1F* m2sigsubphi;

	//ETA
	TH1F* m2sigeta = new TH1F("m2sigeta", "signal region", bins, -3.15, 3.15);
	TH1F* m2sbeta = new TH1F("m2sbeta", "sideband region", bins, -3.15, 3.15);
	TH1F* m2sigsubeta;

	//LOAD SIG SUBBED HISTOGRAMS
    	TFile* subtracted = TFile::Open("/home/jrowell/CMSSW_9_2_1/src/DarkPhotonAnalysisV2/DimuonAnalysis/macros/subs.root");
	TH1F* pthisto = (TH1F*) subtracted->Get("sigsubpt");
	TH1F* phihisto = (TH1F*) subtracted->Get("sigsubphi");
	TH1F* etahisto = (TH1F*) subtracted->Get("sigsubeta");
	TH1F* m1pthisto = (TH1F*) subtracted->Get("m1sigsubpt");
	TH1F* m1phihisto = (TH1F*) subtracted->Get("m1sigsubphi");
	TH1F* m1etahisto = (TH1F*) subtracted->Get("m1sigsubeta");
	TH1F* m2pthisto = (TH1F*) subtracted->Get("m2sigsubpt");
	TH1F* m2phihisto = (TH1F*) subtracted->Get("m2sigsubphi");
	TH1F* m2etahisto = (TH1F*) subtracted->Get("m2sigsubeta");
	int nbins = pthisto->GetXaxis()->GetNbins();

	//FILL WITH THE VARIABLE RANGE FOR EACH MASS PICTURE
	std::vector <double> pt_bin_avg;
	std::vector <double> phi_bin_avg;
	std::vector <double> eta_bin_avg;
	std::vector <double> m1pt_bin_avg;
	std::vector <double> m1phi_bin_avg;
	std::vector <double> m1eta_bin_avg;
	std::vector <double> m2pt_bin_avg;
	std::vector <double> m2phi_bin_avg;
	std::vector <double> m2eta_bin_avg;

	//COMPLETE MASS PICTURE
	TH1F* mass_total = new TH1F("mass total", "mass total", 1000, 0., 11.);
	TH1F* mass_pt_scaled = new TH1F("pt scaled", "pt scaled", 1000, 0., 11.);
	TH1F* mass_eta_scaled = new TH1F("eta scaled", "eta scaled", 1000, 0., 11.);

	//SAVED SCALES FROM PT AND ETA
	std::vector <double> pt_scales_saved, eta_scales_saved;
	//FOR JPSI RANGE
	if (reso == "jpsi"){
		pt_scales_saved.push_back(1.0);
		pt_scales_saved.push_back(0.997836);
		pt_scales_saved.push_back(1.00052);
		pt_scales_saved.push_back(1.00027);
		pt_scales_saved.push_back(1.0);
		pt_scales_saved.push_back(1.0001);
		pt_scales_saved.push_back(0.999923);
		pt_scales_saved.push_back(0.999977);
		pt_scales_saved.push_back(0.99996);
		pt_scales_saved.push_back(1.00047);
		pt_scales_saved.push_back(0.999653);
		pt_scales_saved.push_back(1.00018);

		eta_scales_saved.push_back(0.809813);
		eta_scales_saved.push_back(1.00321);
		eta_scales_saved.push_back(1.00252);
		eta_scales_saved.push_back(1.00163);
		eta_scales_saved.push_back(1.00171);
		eta_scales_saved.push_back(1.00368);
		eta_scales_saved.push_back(1.00375);
		eta_scales_saved.push_back(1.00184);
		eta_scales_saved.push_back(1.00183);
		eta_scales_saved.push_back(1.00235);
		eta_scales_saved.push_back(1.00246);
		eta_scales_saved.push_back(0.98278);
	}
	//FOR UPSILON1S
	if (reso == "upsilon1"){
		pt_scales_saved.push_back(1.00077);
		pt_scales_saved.push_back(0.999892);
		pt_scales_saved.push_back(1.00026);
		pt_scales_saved.push_back(0.99945);
		pt_scales_saved.push_back(1.00059);
		pt_scales_saved.push_back(1.00072);
		pt_scales_saved.push_back(1.00016);
		pt_scales_saved.push_back(1.00077);
		pt_scales_saved.push_back(0.998164);
		pt_scales_saved.push_back(1.00322);
		pt_scales_saved.push_back(0.999801);
		pt_scales_saved.push_back(0.996914);

		eta_scales_saved.push_back(1.00192);
		eta_scales_saved.push_back(1.001);
		eta_scales_saved.push_back(1.00253);
		eta_scales_saved.push_back(1.00102);
		eta_scales_saved.push_back(1.00233);
		eta_scales_saved.push_back(1.00219);
		eta_scales_saved.push_back(1.00158);
		eta_scales_saved.push_back(1.00186);
		eta_scales_saved.push_back(1.00049);
		eta_scales_saved.push_back(1.00217);
		eta_scales_saved.push_back(0.999558);
		eta_scales_saved.push_back(1.00077);

	}
	//FOR PHI
	if (reso == "phi"){
		pt_scales_saved.push_back(1.0);
		pt_scales_saved.push_back(1.0);
		pt_scales_saved.push_back(1.00188);
		pt_scales_saved.push_back(0.999463);
		pt_scales_saved.push_back(1.00118);
		pt_scales_saved.push_back(0.999054);
		pt_scales_saved.push_back(0.999658);
		pt_scales_saved.push_back(1.00572);
		pt_scales_saved.push_back(1.00214);
		pt_scales_saved.push_back(0.998208);
		pt_scales_saved.push_back(1.00128);
		pt_scales_saved.push_back(1.00797);

		eta_scales_saved.push_back(1.0);
		eta_scales_saved.push_back(1.00617);
		eta_scales_saved.push_back(1.0024);
		eta_scales_saved.push_back(1.00174);
		eta_scales_saved.push_back(1.00103);
		eta_scales_saved.push_back(1.00306);
		eta_scales_saved.push_back(1.00366);
		eta_scales_saved.push_back(1.00233);
		eta_scales_saved.push_back(1.00279);
		eta_scales_saved.push_back(1.00214);
		eta_scales_saved.push_back(1.00378);
		eta_scales_saved.push_back(1.0);
	}
	double eta_scale_factor;
	double pt_scale_factor;

	//STORES THE 12 MASS PICTURES FOR EACH VARIABLE
	std::vector <TH1F*> pt_mass_pics;
	TH1F* pt1 = new TH1F("pt1", "pt1", 1000, 0., 11.);
	TH1F* pt2 = new TH1F("pt2", "pt2", 1000, 0., 11.);
	TH1F* pt3 = new TH1F("pt3", "pt3", 1000, 0., 11.);
	TH1F* pt4 = new TH1F("pt4", "pt4", 1000, 0., 11.);
	TH1F* pt5 = new TH1F("pt5", "pt5", 1000, 0., 11.);
	TH1F* pt6 = new TH1F("pt6", "pt6", 1000, 0., 11.);
	TH1F* pt7 = new TH1F("pt7", "pt7", 1000, 0., 11.);
	TH1F* pt8 = new TH1F("pt8", "pt8", 1000, 0., 11.);
	TH1F* pt9 = new TH1F("pt9", "pt9", 1000, 0., 11.);
	TH1F* pt10 = new TH1F("pt10", "pt10", 1000, 0., 11.);
	TH1F* pt11 = new TH1F("pt11", "pt11", 1000, 0., 11.);
	TH1F* pt12 = new TH1F("pt12", "pt12", 1000, 0., 11.);
	pt_mass_pics.push_back(pt1);
	pt_mass_pics.push_back(pt2);
	pt_mass_pics.push_back(pt3);
	pt_mass_pics.push_back(pt4);
	pt_mass_pics.push_back(pt5);
	pt_mass_pics.push_back(pt6);
	pt_mass_pics.push_back(pt7);
	pt_mass_pics.push_back(pt8);
	pt_mass_pics.push_back(pt9);
	pt_mass_pics.push_back(pt10);
	pt_mass_pics.push_back(pt11);
	pt_mass_pics.push_back(pt12);

	std::vector <TH1F*> phi_mass_pics;
	TH1F* phi1 = new TH1F("phi1", "phi1", 1000, 0., 11.);
	TH1F* phi2 = new TH1F("phi2", "phi2", 1000, 0., 11.);
	TH1F* phi3 = new TH1F("phi3", "phi3", 1000, 0., 11.);
	TH1F* phi4 = new TH1F("phi4", "phi4", 1000, 0., 11.);
	TH1F* phi5 = new TH1F("phi5", "phi5", 1000, 0., 11.);
	TH1F* phi6 = new TH1F("phi6", "phi6", 1000, 0., 11.);
	TH1F* phi7 = new TH1F("phi7", "phi7", 1000, 0., 11.);
	TH1F* phi8 = new TH1F("phi8", "phi8", 1000, 0., 11.);
	TH1F* phi9 = new TH1F("phi9", "phi9", 1000, 0., 11.);
	TH1F* phi10 = new TH1F("phi10", "phi10", 1000, 0., 11.);
	TH1F* phi11 = new TH1F("phi11", "phi11", 1000, 0., 11.);
	TH1F* phi12 = new TH1F("phi12", "phi12", 1000, 0., 11.);
	phi_mass_pics.push_back(phi1);
	phi_mass_pics.push_back(phi2);
	phi_mass_pics.push_back(phi3);
	phi_mass_pics.push_back(phi4);
	phi_mass_pics.push_back(phi5);
	phi_mass_pics.push_back(phi6);
	phi_mass_pics.push_back(phi7);
	phi_mass_pics.push_back(phi8);
	phi_mass_pics.push_back(phi9);
	phi_mass_pics.push_back(phi10);
	phi_mass_pics.push_back(phi11);
	phi_mass_pics.push_back(phi12);

	std::vector <TH1F*> eta_mass_pics;
	TH1F* eta1 = new TH1F("eta1", "eta1", 1000, 0., 11.);
	TH1F* eta2 = new TH1F("eta2", "eta2", 1000, 0., 11.);
	TH1F* eta3 = new TH1F("eta3", "eta3", 1000, 0., 11.);
	TH1F* eta4 = new TH1F("eta4", "eta4", 1000, 0., 11.);
	TH1F* eta5 = new TH1F("eta5", "eta5", 1000, 0., 11.);
	TH1F* eta6 = new TH1F("eta6", "eta6", 1000, 0., 11.);
	TH1F* eta7 = new TH1F("eta7", "eta7", 1000, 0., 11.);
	TH1F* eta8 = new TH1F("eta8", "eta8", 1000, 0., 11.);
	TH1F* eta9 = new TH1F("eta9", "eta9", 1000, 0., 11.);
	TH1F* eta10 = new TH1F("eta10", "eta10", 1000, 0., 11.);
	TH1F* eta11 = new TH1F("eta11", "eta11", 1000, 0., 11.);
	TH1F* eta12 = new TH1F("eta12", "eta12", 1000, 0., 11.);
	eta_mass_pics.push_back(eta1);
	eta_mass_pics.push_back(eta2);
	eta_mass_pics.push_back(eta3);
	eta_mass_pics.push_back(eta4);
	eta_mass_pics.push_back(eta5);
	eta_mass_pics.push_back(eta6);
	eta_mass_pics.push_back(eta7);
	eta_mass_pics.push_back(eta8);
	eta_mass_pics.push_back(eta9);
	eta_mass_pics.push_back(eta10);
	eta_mass_pics.push_back(eta11);
	eta_mass_pics.push_back(eta12);

	std::vector <TH1F*> m2pt_mass_pics;
	TH1F* m2pt1 = new TH1F("m2pt1", "m2pt1", 1000, 0., 11.);
	TH1F* m2pt2 = new TH1F("m2pt2", "m2pt2", 1000, 0., 11.);
	TH1F* m2pt3 = new TH1F("m2pt3", "m2pt3", 1000, 0., 11.);
	TH1F* m2pt4 = new TH1F("m2pt4", "m2pt4", 1000, 0., 11.);
	TH1F* m2pt5 = new TH1F("m2pt5", "m2pt5", 1000, 0., 11.);
	TH1F* m2pt6 = new TH1F("m2pt6", "m2pt6", 1000, 0., 11.);
	TH1F* m2pt7 = new TH1F("m2pt7", "m2pt7", 1000, 0., 11.);
	TH1F* m2pt8 = new TH1F("m2pt8", "m2pt8", 1000, 0., 11.);
	TH1F* m2pt9 = new TH1F("m2pt9", "m2pt9", 1000, 0., 11.);
	TH1F* m2pt10 = new TH1F("m2pt10", "m2pt10", 1000, 0., 11.);
	TH1F* m2pt11 = new TH1F("m2pt11", "m2pt11", 1000, 0., 11.);
	TH1F* m2pt12 = new TH1F("m2pt12", "m2pt12", 1000, 0., 11.);
	m2pt_mass_pics.push_back(m2pt1);
	m2pt_mass_pics.push_back(m2pt2);
	m2pt_mass_pics.push_back(m2pt3);
	m2pt_mass_pics.push_back(m2pt4);
	m2pt_mass_pics.push_back(m2pt5);
	m2pt_mass_pics.push_back(m2pt6);
	m2pt_mass_pics.push_back(m2pt7);
	m2pt_mass_pics.push_back(m2pt8);
	m2pt_mass_pics.push_back(m2pt9);
	m2pt_mass_pics.push_back(m2pt10);
	m2pt_mass_pics.push_back(m2pt11);
	m2pt_mass_pics.push_back(m2pt12);

	std::vector <TH1F*> m1pt_mass_pics;
	TH1F* m1pt1 = new TH1F("m1pt1", "m1pt1", 1000, 0., 11.);
	TH1F* m1pt2 = new TH1F("m1pt2", "m1pt2", 1000, 0., 11.);
	TH1F* m1pt3 = new TH1F("m1pt3", "m1pt3", 1000, 0., 11.);
	TH1F* m1pt4 = new TH1F("m1pt4", "m1pt4", 1000, 0., 11.);
	TH1F* m1pt5 = new TH1F("m1pt5", "m1pt5", 1000, 0., 11.);
	TH1F* m1pt6 = new TH1F("m1pt6", "m1pt6", 1000, 0., 11.);
	TH1F* m1pt7 = new TH1F("m1pt7", "m1pt7", 1000, 0., 11.);
	TH1F* m1pt8 = new TH1F("m1pt8", "m1pt8", 1000, 0., 11.);
	TH1F* m1pt9 = new TH1F("m1pt9", "m1pt9", 1000, 0., 11.);
	TH1F* m1pt10 = new TH1F("m1pt10", "m1pt10", 1000, 0., 11.);
	TH1F* m1pt11 = new TH1F("m1pt11", "m1pt11", 1000, 0., 11.);
	TH1F* m1pt12 = new TH1F("m1pt12", "m1pt12", 1000, 0., 11.);
	m1pt_mass_pics.push_back(m1pt1);
	m1pt_mass_pics.push_back(m1pt2);
	m1pt_mass_pics.push_back(m1pt3);
	m1pt_mass_pics.push_back(m1pt4);
	m1pt_mass_pics.push_back(m1pt5);
	m1pt_mass_pics.push_back(m1pt6);
	m1pt_mass_pics.push_back(m1pt7);
	m1pt_mass_pics.push_back(m1pt8);
	m1pt_mass_pics.push_back(m1pt9);
	m1pt_mass_pics.push_back(m1pt10);
	m1pt_mass_pics.push_back(m1pt11);
	m1pt_mass_pics.push_back(m1pt12);	

	std::vector <TH1F*> m1eta_mass_pics;
	TH1F* m1eta1 = new TH1F("m1eta1", "m1eta1", 1000, 0., 11.);
	TH1F* m1eta2 = new TH1F("m1eta2", "m1eta2", 1000, 0., 11.);
	TH1F* m1eta3 = new TH1F("m1eta3", "m1eta3", 1000, 0., 11.);
	TH1F* m1eta4 = new TH1F("m1eta4", "m1eta4", 1000, 0., 11.);
	TH1F* m1eta5 = new TH1F("m1eta5", "m1eta5", 1000, 0., 11.);
	TH1F* m1eta6 = new TH1F("m1eta6", "m1eta6", 1000, 0., 11.);
	TH1F* m1eta7 = new TH1F("m1eta7", "m1eta7", 1000, 0., 11.);
	TH1F* m1eta8 = new TH1F("m1eta8", "m1eta8", 1000, 0., 11.);
	TH1F* m1eta9 = new TH1F("m1eta9", "m1eta9", 1000, 0., 11.);
	TH1F* m1eta10 = new TH1F("m1eta10", "m1eta10", 1000, 0., 11.);
	TH1F* m1eta11 = new TH1F("m1eta11", "m1eta11", 1000, 0., 11.);
	TH1F* m1eta12 = new TH1F("m1eta12", "m1eta12", 1000, 0., 11.);
	m1eta_mass_pics.push_back(m1eta1);
	m1eta_mass_pics.push_back(m1eta2);
	m1eta_mass_pics.push_back(m1eta3);
	m1eta_mass_pics.push_back(m1eta4);
	m1eta_mass_pics.push_back(m1eta5);
	m1eta_mass_pics.push_back(m1eta6);
	m1eta_mass_pics.push_back(m1eta7);
	m1eta_mass_pics.push_back(m1eta8);
	m1eta_mass_pics.push_back(m1eta9);
	m1eta_mass_pics.push_back(m1eta10);
	m1eta_mass_pics.push_back(m1eta11);
	m1eta_mass_pics.push_back(m1eta12);

	std::vector <TH1F*> m2eta_mass_pics;
	TH1F* m2eta1 = new TH1F("m2eta1", "m2eta1", 1000, 0., 11.);
	TH1F* m2eta2 = new TH1F("m2eta2", "m2eta2", 1000, 0., 11.);
	TH1F* m2eta3 = new TH1F("m2eta3", "m2eta3", 1000, 0., 11.);
	TH1F* m2eta4 = new TH1F("m2eta4", "m2eta4", 1000, 0., 11.);
	TH1F* m2eta5 = new TH1F("m2eta5", "m2eta5", 1000, 0., 11.);
	TH1F* m2eta6 = new TH1F("m2eta6", "m2eta6", 1000, 0., 11.);
	TH1F* m2eta7 = new TH1F("m2eta7", "m2eta7", 1000, 0., 11.);
	TH1F* m2eta8 = new TH1F("m2eta8", "m2eta8", 1000, 0., 11.);
	TH1F* m2eta9 = new TH1F("m2eta9", "m2eta9", 1000, 0., 11.);
	TH1F* m2eta10 = new TH1F("m2eta10", "m2eta10", 1000, 0., 11.);
	TH1F* m2eta11 = new TH1F("m2eta11", "m2eta11", 1000, 0., 11.);
	TH1F* m2eta12 = new TH1F("m2eta12", "m2eta12", 1000, 0., 11.);
	m2eta_mass_pics.push_back(m2eta1);
	m2eta_mass_pics.push_back(m2eta2);
	m2eta_mass_pics.push_back(m2eta3);
	m2eta_mass_pics.push_back(m2eta4);
	m2eta_mass_pics.push_back(m2eta5);
	m2eta_mass_pics.push_back(m2eta6);
	m2eta_mass_pics.push_back(m2eta7);
	m2eta_mass_pics.push_back(m2eta8);
	m2eta_mass_pics.push_back(m2eta9);
	m2eta_mass_pics.push_back(m2eta10);
	m2eta_mass_pics.push_back(m2eta11);
	m2eta_mass_pics.push_back(m2eta12);

	//READS FROM THE TREE 
	cout<<treepath<<endl;
	TChain* chain = new TChain("mmtree/tree");
	chain->Add((TString)treepath);
	TTreeReader reader(chain);
	TTreeReaderValue<std::vector<float> >          	mpt      (reader, "muonpt"     );
	TTreeReaderValue<std::vector<float> >          	meta     (reader, "muoneta"    );
	TTreeReaderValue<std::vector<float> >          	mphi     (reader, "muonphi"    );
	TTreeReaderValue<std::vector<float> >          	mcharge  (reader, "muoncharge" );
	TTreeReaderValue<std::vector<char>  >          	mid      (reader, "muonid"     );
	TTreeReaderValue<std::vector<float> >          	chi2     (reader, "chi2"       );
	TTreeReaderValue<std::vector<float> >          	dxy      (reader, "dxy"        );
	TTreeReaderValue<std::vector<float> >          	dz       (reader, "dz"         );
	TTreeReaderValue<std::vector<float> >          	cpiso    (reader, "cpiso"      );
	TTreeReaderValue<std::vector<float> >          	chiso    (reader, "chiso"      );
	TTreeReaderValue<std::vector<float> >          	phiso    (reader, "phiso"      );
	TTreeReaderValue<std::vector<float> >          	nhiso    (reader, "nhiso"      );
	TTreeReaderValue<std::vector<float> >          	puiso    (reader, "puiso"      );
	TTreeReaderValue<std::vector<float> >          	tkiso    (reader, "tkiso"      );
	TTreeReaderValue<std::vector<bool> >           	l1Result (reader, "l1Result"   );
	TTreeReaderValue<std::vector<unsigned char> >  	nmhits   (reader, "nMuonHits"  );
	TTreeReaderValue<std::vector<unsigned char> >  	nphits   (reader, "nPixelHits" );
	TTreeReaderValue<std::vector<unsigned char> >  	ntklayers(reader, "nTkLayers"  );
	TTreeReaderValue<unsigned char>                	hlt      (reader, "trig"       );
	TTreeReaderValue<unsigned>                     	nverts   (reader, "nvtx"       );
	TTreeReaderValue<std::vector<float> >          	vtxX     (reader, "vtxX"       );
	TTreeReaderValue<std::vector<float> >          	vtxY     (reader, "vtxY"       );
	TTreeReaderValue<std::vector<float> >          	vtxZ     (reader, "vtxZ"       );
	TTreeReaderValue<double>                       	rho      (reader, "rho"        );
	TTreeReaderValue<unsigned int>                 	run      (reader, "run"        );
	TTreeReaderValue<unsigned int>		       	lumSec   (reader, "lumSec"     );
     
	//FOR MUON VARS
	float mmpt = 0;
	float mmphi = 0;
	float mmeta = 0;
	float  m1pt    = 0.0;        
	float  m1eta   = 0.0;        
	float  m1phi   = 0.0;        
	float  m1iso   = 0.0;        
	float  m2pt    = 0.0;        
	float  m2eta   = 0.0;        
	float  m2phi   = 0.0;        
	float  m2iso   = 0.0;        
	float  mass    = 0.0;        
	float  mass4   = 0.0;        
	char   m1id    = 0;
	char   m2id    = 0;
	char   m3id    = 0;
	char   m4id    = 0;
	float  m1ch    = 0.; 
	float  m2ch    = 0.;
	float  m3ch    = 0.; 
	float  m4ch    = 0.; 
	vector<bool> l1Bools;
	unsigned nvtx  = 0;
	unsigned Run   = 0;
	unsigned LumSec   = 0;

	while(reader.Next()) {
		if (((*hlt) & 2) == 0) continue;   
		bool passIso=false;
		bool passIsoLoose=false;
		double ea = 0.45;
		std::vector<unsigned> goodmuons;
		for (std::size_t i = 0; i < mpt->size(); i++) {
			if ((*nphits)[i] == 0)     continue;
			if ((*ntklayers)[i] <= 5)  continue;
			if ((*chi2)[i] > 10.)      continue;
			double iso = (*cpiso)[i] + (*nhiso)[i] + (*phiso)[i] - ea*(*rho);
			goodmuons.push_back(i);
		}

		if (goodmuons.size() < 2) continue;
		unsigned idx1 = goodmuons[0];
		unsigned idx2 = goodmuons[1];
		if((*tkiso)[goodmuons[0]]<0.02 && (*tkiso)[goodmuons[1]]<0.02 ) passIso=true;
		if((*tkiso)[goodmuons[0]]<0.15 && (*tkiso)[goodmuons[1]]<0.15 ) passIsoLoose=true;
		if ((*mpt)[goodmuons[0]] < (*mpt)[goodmuons[1]]) {
			idx1 = goodmuons[1];
			idx2 = goodmuons[0];
		}
		
		//BUILDS LORENTZ VECTORS FOR DIMUON AND MUONS 1,2
		TLorentzVector mm;
		TLorentzVector m1;
		TLorentzVector m2;
		m1.SetPtEtaPhiM((*mpt)[idx1], (*meta)[idx1], (*mphi)[idx1], 0.1057);
		m2.SetPtEtaPhiM((*mpt)[idx2], (*meta)[idx2], (*mphi)[idx2], 0.1057);
		mm += m1;
		mm += m2;

		//CREATES VARIABLES FOR M1, M2
		m1pt   = m1.Pt();
		m1eta  = m1.Eta();
		m1phi  = m1.Phi();
		m1iso  = (*cpiso)[idx1] + (*phiso)[idx1] + (*nhiso)[idx1] - ea*(*rho);
		m1id   = (*mid)[idx1];
		m1ch   = (*mcharge)[idx1];

		m2pt   = m2.Pt();
		m2eta  = m2.Eta();
		m2phi  = m2.Phi();
		m2iso  = (*cpiso)[idx2] + (*phiso)[idx2] + (*nhiso)[idx2] - ea*(*rho);
		m2id   = (*mid)[idx2];
		m2ch   = (*mcharge)[idx2];

		//STORE DIMUON VARS
		mmpt = mm.Pt();
		mmphi = mm.Phi();
		mmeta= mm.Eta();

		Run=*run;
		LumSec=*lumSec;

		mass   = mm.M();

		//FILL SIGNAL AND SIDEBAND DIMUON PT PLOTS
		double sig_lo;
		double sig_hi;
		double sb_l_lo;
		double sb_l_hi;
		double sb_r_lo;
		double sb_r_hi;
		if (reso == "jpsi"){
			sig_lo = 3.0;//2.95;
			sig_hi = 3.19;//3.25;
			sb_l_lo = 2.71;//2.6;
			sb_l_hi = 2.9;
			sb_r_lo = 3.29;//3.3;
			sb_r_hi = 3.48;//3.6;
		}
		if (reso == "upsilon1"){
			sig_lo = 9.26;
			sig_hi = 9.66;
			sb_l_lo = 8.5;
			sb_l_hi = 8.96;
			sb_r_lo = 10.54;
			sb_r_hi = 11.0;
		}
		if (reso == "phi"){
			sig_lo = 0.98;
			sig_hi = 1.06;
			sb_l_lo = 0.88;
			sb_l_hi = 0.96;
			sb_r_lo = 1.08;
			sb_r_hi = 1.16;
		}

		if (mass>=sig_lo && mass<=sig_hi) {
			sigpt->Fill(mmpt);
			m1sigpt->Fill(m1pt);
			m2sigpt->Fill(m2pt);
			sigphi->Fill(mmphi);
			m1sigphi->Fill(m1phi);
			m2sigphi->Fill(m2phi);
			sigeta->Fill(mmeta);
			m1sigeta->Fill(m1eta);
			m2sigeta->Fill(m2eta);
		}
		if ( (mass>=sb_l_lo && mass<=sb_l_hi) || (mass>=sb_r_lo && mass<=sb_r_hi) ) {
			sbpt->Fill(mmpt);
			m1sbpt->Fill(m1pt);
			m2sbpt->Fill(m2pt);
			sbphi->Fill(mmphi);
			m1sbphi->Fill(m1phi);
			m2sbphi->Fill(m2phi);
			sbeta->Fill(mmeta);
			m1sbeta->Fill(m1eta);
			m2sbeta->Fill(m2eta);
		}

		//CREATE MASS PICTURE FOR EACH PT BIN
		double bin_lo, bin_hi, bin_min, bin_max;
		//MMETA
		for (int i = 1; i <= nbins-1; i++){
			bin_lo = etahisto->GetBinCenter(i);
			bin_hi = etahisto->GetBinCenter(i+1);

			bin_min = etahisto->GetBinCenter(1);
			bin_max = etahisto->GetBinCenter(13);

			if (mmeta>=bin_lo && mmeta<bin_hi){
				eta_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" Eta").c_str());
				eta_mass_pics[i-1]->SetStats(false);
				eta_mass_pics[i-1]->Fill(mass);
				eta_scale_factor = eta_scales_saved[i-1];
				break;
			} else if(mmeta<bin_min || mmeta >= bin_max) {
				eta_scale_factor = 1.0;
				break;
			}
		}
		//MMPT
		for (int i = 1; i < nbins; i++){
			//Get bin centers for bounds of the mass pictures
			bin_lo = pthisto->GetBinCenter(i);
			bin_hi = pthisto->GetBinCenter(i+1);

			bin_min = pthisto->GetBinCenter(1);
			bin_max = pthisto->GetBinCenter(13);

			//If mmpt is within the range of the first mass pic, fill with mass
			if (mmpt>=bin_lo && mmpt<bin_hi){
				pt_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" GeV pt").c_str());
				pt_mass_pics[i-1]->SetStats(false);
				pt_mass_pics[i-1]->Fill(mass*eta_scale_factor);//FILLS PT BINS WITH MASSES CORRECTED FOR ETA, THIS IS THEN USED TO GENERATE NEW PT SCALES
				pt_scale_factor = pt_scales_saved[i-1];
				break;
			} else if(mmpt<bin_min || mmpt >= bin_max) {
				pt_scale_factor = 1.0;
				break;
			}
		}
	
		mass_total->Fill(mass);
		mass_eta_scaled->Fill(mass*eta_scale_factor);
		mass_pt_scaled->Fill(mass*eta_scale_factor*pt_scale_factor); //NOW THIS IS SCALED FOR BOTH DIMUON ETA AND DIMUON PT

		//M1ETA
		for (int i = 1; i < nbins; i++){
			bin_lo = m1etahisto->GetBinCenter(i);
			bin_hi = m1etahisto->GetBinCenter(i+1);

			m1eta_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (m1eta>=bin_lo && m1eta<bin_hi){
				m1eta_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" m1eta").c_str());
				m1eta_mass_pics[i-1]->SetStats(false);
				m1eta_mass_pics[i-1]->Fill(mass);
			}
		}
		//MMPHI
		for (int i = 1; i <= nbins-1; i++){
			bin_lo = phihisto->GetBinCenter(i);
			bin_hi = phihisto->GetBinCenter(i+1);
			phi_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (mmphi>=bin_lo && mmphi<bin_hi){
				phi_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" Phi").c_str());
				phi_mass_pics[i-1]->SetStats(false);
				phi_mass_pics[i-1]->Fill(mass);
			}
		}
		//M1PT
		for (int i = 1; i <= nbins-1; i++){
			bin_lo = m1pthisto->GetBinCenter(i);
			bin_hi = m1pthisto->GetBinCenter(i+1);
			m1pt_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (m1pt>=bin_lo && m1pt<bin_hi){
				m1pt_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" m1pt").c_str());
				m1pt_mass_pics[i-1]->SetStats(false);
				m1pt_mass_pics[i-1]->Fill(mass);
			}
		}
		//M2PT
		for (int i = 1; i <= nbins-1; i++){
			bin_lo = m2pthisto->GetBinCenter(i);
			bin_hi = m2pthisto->GetBinCenter(i+1);
			m2pt_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (m2pt>=bin_lo && m2pt<bin_hi){
				m2pt_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" m2pt").c_str());
				m2pt_mass_pics[i-1]->SetStats(false);
				m2pt_mass_pics[i-1]->Fill(mass);
			}
		}
		//M2ETA
		for (int i = 1; i <= nbins-1; i++){
			bin_lo = m2etahisto->GetBinCenter(i);
			bin_hi = m2etahisto->GetBinCenter(i+1);
			m2eta_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (m2eta>=bin_lo && m2eta<bin_hi){
				m2eta_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" m2eta").c_str());
				m2eta_mass_pics[i-1]->SetStats(false);
				m2eta_mass_pics[i-1]->Fill(mass);
			}
		}
	}

	//WRITE SIGNAL AND SIDEBAND PLOTS
	TFile *subs = new TFile("subs.root", "RECREATE");
	//DIMUON
	sigsubpt = (TH1F*) sigpt->Clone("sigsubpt");
	sigsubpt->Add(sbpt, -.5);
	sigsubpt->SetTitle("signal subtracted pt");

	sigsubphi = (TH1F*) sigphi->Clone("sigsubphi");
	sigsubphi->Add(sbphi, -.5);
	sigsubphi->SetTitle("signal subtracted phi");

	sigsubeta = (TH1F*) sigeta->Clone("sigsubeta");
	sigsubeta->Add(sbeta, -.5);
	sigsubeta->SetTitle("signal subtracted eta");

	//MUON 1
	m1sigsubpt = (TH1F*) m1sigpt->Clone("m1sigsubpt");
	m1sigsubpt->Add(m1sbpt, -.5);
	m1sigsubpt->SetTitle("signal subtracted pt muon 1");

	m1sigsubphi = (TH1F*) m1sigphi->Clone("m1sigsubphi");
	m1sigsubphi->Add(m1sbphi, -.5);
	m1sigsubphi->SetTitle("signal subtracted phi muon 1");

	m1sigsubeta = (TH1F*) m1sigeta->Clone("m1sigsubeta");
	m1sigsubeta->Add(m1sbeta, -.5);
	m1sigsubeta->SetTitle("signal subtracted eta muon 1");
	
	//MUON 2
	m2sigsubpt = (TH1F*) m2sigpt->Clone("m2sigsubpt");
	m2sigsubpt->Add(m2sbpt, -.5);
	m2sigsubpt->SetTitle("signal subtracted pt muon 2");

	m2sigsubphi = (TH1F*) m2sigphi->Clone("m2sigsubphi");
	m2sigsubphi->Add(m2sbphi, -.5);
	m2sigsubphi->SetTitle("signal subtracted phi muon 2");

	m2sigsubeta = (TH1F*) sigeta->Clone("m2sigsubeta");
	m2sigsubeta->Add(m2sbeta, -.5);
	m2sigsubeta->SetTitle("signal subtracted eta muon 2");
	
	//PLOT THE SIG, SB, and SUB
	//DIMUON
	TCanvas c1("c1", "c1");
	sigsubpt->SetStats(false);
	sigsubpt->SetLineColor(kBlue);
	sigsubpt->GetXaxis()->SetTitle("GeV");
	sigsubpt->GetYaxis()->SetTitle("Events");
	sigsubpt->Draw("SAME");
	c1.SaveAs((reso+"/sig_subs/sig_sub_pt.png").c_str());

	TCanvas c2("c2", "c2");
	sigsubphi->SetStats(false);
	sigsubphi->SetLineColor(kBlue);
	sigsubphi->GetXaxis()->SetTitle("Phi");
	sigsubphi->GetYaxis()->SetTitle("Events");
	sigsubphi->Draw("SAME");
	c2.SaveAs((reso+"/sig_subs/sig_sub_phi.png").c_str());

	TCanvas c3("c3", "c3");
	sigsubeta->SetStats(false);
	sigsubeta->SetLineColor(kBlue);
	sigsubeta->GetXaxis()->SetTitle("Eta");
	sigsubeta->GetYaxis()->SetTitle("Events");
	sigsubeta->Draw("SAME");
	c3.SaveAs((reso+"/sig_subs/sig_sub_eta.png").c_str());

	//MUON 1
	TCanvas c4("c4", "c4");
	m1sigsubpt->SetStats(false);
	m1sigsubpt->SetLineColor(kBlue);
	m1sigsubpt->GetXaxis()->SetTitle("GeV");
	m1sigsubpt->GetYaxis()->SetTitle("Events");
	m1sigsubpt->Draw("SAME");
	c4.SaveAs((reso+"/sig_subs/m1sig_sub_pt.png").c_str());

	TCanvas c5("c5", "c5");
	m1sigsubphi->SetStats(false);
	m1sigsubphi->SetLineColor(kBlue);
	m1sigsubphi->GetXaxis()->SetTitle("Phi");
	m1sigsubphi->GetYaxis()->SetTitle("Events");
	m1sigsubphi->Draw("SAME");
	c5.SaveAs((reso+"/sig_subs/m1sig_sub_phi.png").c_str());

	TCanvas c6("c6", "c6");
	m1sigsubeta->SetStats(false);
	m1sigsubeta->SetLineColor(kBlue);
	m1sigsubeta->GetXaxis()->SetTitle("Eta");
	m1sigsubeta->GetYaxis()->SetTitle("Events");
	m1sigsubeta->Draw("SAME");
	c6.SaveAs((reso+"/sig_subs/m1sig_sub_eta.png").c_str());

	//MUON 2
	TCanvas c7("c7", "c7");
	m2sigsubpt->SetStats(false);
	m2sigsubpt->SetLineColor(kBlue);
	m2sigsubpt->GetXaxis()->SetTitle("GeV");
	m2sigsubpt->GetYaxis()->SetTitle("Events");
	m2sigsubpt->Draw("SAME");
	c7.SaveAs((reso+"/sig_subs/m2sig_sub_pt.png").c_str());

	TCanvas c8("c8", "c8");
	m2sigsubphi->SetStats(false);
	m2sigsubphi->SetLineColor(kBlue);
	m2sigsubphi->GetXaxis()->SetTitle("Phi");
	m2sigsubphi->GetYaxis()->SetTitle("Events");
	m2sigsubphi->Draw("SAME");
	c8.SaveAs((reso+"/sig_subs/m2sig_sub_phi.png").c_str());

	TCanvas c9("c9", "c9");
	m2sigsubeta->SetStats(false);
	m2sigsubeta->SetLineColor(kBlue);
	m2sigsubeta->GetXaxis()->SetTitle("Eta");
	m2sigsubeta->GetYaxis()->SetTitle("Events");
	m2sigsubeta->Draw("SAME");
	c9.SaveAs((reso+"/sig_subs/m2sig_sub_eta.png").c_str());

	//SIGNAL SUBTRACTED PT HISTOGRAM
	sigsubpt->Write();
	sigsubphi->Write();
	sigsubeta->Write();
	m1sigsubpt->Write();
	m1sigsubphi->Write();
	m1sigsubeta->Write();
	m2sigsubpt->Write();
	m2sigsubphi->Write();
	m2sigsubeta->Write();
	subs->Close();

	//WRITE THE MASS PICTURES FROM THE SIDEBAND PLOTS
	std::vector<double> ptmeans;
	std::vector<double> pterr;
	std::vector<double> pt0;
	std::vector<double> ptmeans_scaled;
	std::vector<double> pt_scales;

	std::vector<double> m1ptmeans;
	std::vector<double> m1pterr;
	std::vector<double> m1pt0;

	std::vector<double> m2ptmeans;
	std::vector<double> m2pterr;
	std::vector<double> m2pt0;

	std::vector<double> phimeans;
	std::vector<double> phierr;
	std::vector<double> phi0;

	std::vector<double> etameans;
	std::vector<double> etaerr;
	std::vector<double> eta0;
	std::vector<double> eta_scales;
	std::vector<double> eta_scaled;

	std::vector<double> m1etameans;
	std::vector<double> m1etaerr;
	std::vector<double> m1eta0;

	std::vector<double> m2etameans;
	std::vector<double> m2etaerr;
	std::vector<double> m2eta0;

	float left_bound;
	float right_bound;
	
	if(reso == "jpsi"){
		left_bound = 3.0;
		right_bound = 3.19;
	}
	if(reso == "upsilon1"){
		left_bound = 9.26;
		right_bound = 9.66;
	}
	if(reso == "phi"){
		left_bound = 0.98;
		right_bound = 1.06;
	}

	int i = 1;
	for (const auto &hist : pt_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");
		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);
		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		fit_func->SetParameter(1, real_mass);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_ptbin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		double scale = real_mass/(fit_func->GetParameter(1));
		pt_scales.push_back(scale);
		ptmeans.push_back(fit_func->GetParameter(1));
		ptmeans_scaled.push_back(fit_func->GetParameter(1)*scale);
		pterr.push_back(fit_func->GetParError(1));
		pt0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : eta_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");
		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);
		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		fit_func->SetParameter(1, real_mass);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_etabin"+std::to_string(i)+".png").c_str());

		double scale = real_mass/fit_func->GetParameter(1);
		eta_scales.push_back(scale);

		cout << fit_func->GetParameter(1) << endl;
		etameans.push_back(fit_func->GetParameter(1));
		eta_scaled.push_back(fit_func->GetParameter(1)*scale);
		etaerr.push_back(fit_func->GetParError(1));
		eta0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : phi_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");

		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_phibin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		phimeans.push_back(fit_func->GetParameter(1));
		phierr.push_back(fit_func->GetParError(1));
		phi0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : m1eta_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");

		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_m1etabin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		m1etameans.push_back(fit_func->GetParameter(1));
		m1etaerr.push_back(fit_func->GetParError(1));
		m1eta0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : m1pt_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");

		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_m1ptbin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		m1ptmeans.push_back(fit_func->GetParameter(1));
		m1pterr.push_back(fit_func->GetParError(1));
		m1pt0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : m2pt_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");

		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_m2ptbin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		m2ptmeans.push_back(fit_func->GetParameter(1));
		m2pterr.push_back(fit_func->GetParError(1));
		m2pt0.push_back(0);
		i++;
	}
	i = 1;
	for (const auto &hist : m2eta_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");

		TF1 *fit_func = new TF1("fit_func", "gaus", left_bound, right_bound);
		hist->Fit("fit_func", "", "", left_bound, right_bound);

		hist->GetXaxis()->SetRangeUser(left_bound, right_bound);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs((reso+"/mass_pictures/mass_picture_m2etabin"+std::to_string(i)+".png").c_str());

		cout << fit_func->GetParameter(1) << endl;
		m2etameans.push_back(fit_func->GetParameter(1));
		m2etaerr.push_back(fit_func->GetParError(1));
		m2eta0.push_back(0);
		i++;
	}

	//CREATES BIN AVERAGE VALUES FOR PLOTTING AGAINST THE DISTRIBUTION MEANS
	for(int i = 1; i < nbins; i++){
		pt_bin_avg.push_back((pthisto->GetBinCenter(i)+pthisto->GetBinCenter(i+1))/2);
		eta_bin_avg.push_back((etahisto->GetBinCenter(i)+etahisto->GetBinCenter(i+1))/2);
	}

	//CREATE PLOTS FOR MEAN OF PEAK WITH CHANGES
	TLine* line = new TLine(-5, real_mass, 30, real_mass);

	TGraphErrors* pt_gr = new TGraphErrors(12, &pt_bin_avg[0], &ptmeans[0], &pt0[0], &pterr[0]);
	TGraphErrors* ptmeans_scaled_gr = new TGraphErrors(12, &pt_bin_avg[0], &ptmeans_scaled[0], &pt0[0], &pterr[0]);

	TGraphErrors* eta_scaled_gr = new TGraphErrors(12, &eta_bin_avg[0], &eta_scaled[0], &eta0[0], &etaerr[0]);
	TGraphErrors* eta_gr = new TGraphErrors(12, &eta_bin_avg[0], &etameans[0], &eta0[0], &etaerr[0]);

	TGraphErrors* m1pt_gr = new TGraphErrors(12, &m1pt_bin_avg[0], &m1ptmeans[0], &m1pt0[0], &m1pterr[0]);
	TGraphErrors* m2pt_gr = new TGraphErrors(12, &m2pt_bin_avg[0], &m2ptmeans[0], &m2pt0[0], &m2pterr[0]);
	TGraphErrors* phi_gr = new TGraphErrors(12, &phi_bin_avg[0], &phimeans[0], &phi0[0], &phierr[0]);
	TGraphErrors* m1eta_gr = new TGraphErrors(12, &m1eta_bin_avg[0], &m1etameans[0], &m1eta0[0], &m1etaerr[0]);
	TGraphErrors* m2eta_gr = new TGraphErrors(12, &m2eta_bin_avg[0], &m2etameans[0], &m2eta0[0], &m2etaerr[0]);


	//MMPT
	TCanvas ca("ca", "ca");
	TLegend* leg1 = new TLegend();
	pt_gr->SetTitle("Jpsi peak mean vs pt bin average; pt bin range average; GeV");
	pt_gr->SetMarkerColor(4);
	pt_gr->SetMarkerStyle(21);
	//pt_gr->GetYaxis()->SetRangeUser(3.08, 3.11);
	pt_gr->Draw("AP");
	ptmeans_scaled_gr->SetMarkerColor(kRed);
	ptmeans_scaled_gr->SetMarkerStyle(21);
	ptmeans_scaled_gr->Draw("P");
	line->Draw("P");
	leg1->AddEntry(pt_gr, "uncorrected");
	leg1->AddEntry(ptmeans_scaled_gr, "corrected");
	leg1->Draw();
	ca.SaveAs((reso+"/mass_means/pt_means.png").c_str());
	//MMETA
	TCanvas cf("cf", "cf");
	TLegend* leg2 = new TLegend();
	eta_gr->SetTitle("JPsi peak mean vs eta bin average; eta bin range average; GeV");
	eta_gr->SetMarkerColor(4);
	eta_gr->SetMarkerStyle(21);
	//eta_gr->GetYaxis()->SetRangeUser(3.08, 3.11);
	eta_gr->Draw("AP");
	eta_scaled_gr->SetMarkerColor(kRed);
	eta_scaled_gr->SetMarkerStyle(21);
	eta_scaled_gr->Draw("P");
	line->Draw("P");
	leg2->AddEntry(eta_gr, "uncorrected");
	leg2->AddEntry(eta_scaled_gr, "correct");
	leg2->Draw();
	cf.SaveAs((reso+"/mass_means/eta_corrected_means.png").c_str());
	//MMPHI
	TCanvas cb("cb", "cb");
	phi_gr->SetTitle("JPsi peak mean vs phi bin average; phi bin range average; GeV");
	phi_gr->SetMarkerColor(4);
	phi_gr->SetMarkerStyle(21);
	//phi_gr->GetYaxis()->SetRangeUser(3.08, 3.09);
	phi_gr->Draw("AP");
	cb.SaveAs((reso+"/mass_means/phi_means.png").c_str());
	//M1ETA
	TCanvas cc("cc", "cc");
	m1eta_gr->SetTitle("JPsi peak mean vs m1eta bin average; m1eta bin range average; GeV");
	m1eta_gr->SetMarkerColor(4);
	m1eta_gr->SetMarkerStyle(21);
	//m1eta_gr->GetYaxis()->SetRangeUser(3.084, 3.092);
	m1eta_gr->Draw("AP");
	cc.SaveAs((reso+"/mass_means/m1eta_means.png").c_str());
	//M1PT
	TCanvas cd("cd", "cd");
	TLegend* leg = new TLegend();
	m1pt_gr->SetTitle("JPsi peak mean vs m1pt bin average; m1pt bin range average; GeV");
	m1pt_gr->SetMarkerColor(4);
	m1pt_gr->SetMarkerStyle(21);
	//m1pt_gr->GetYaxis()->SetRangeUser(3.078, 3.095);
	m1pt_gr->Draw("AP");
	cd.SaveAs((reso+"/mass_means/m1pt_means.png").c_str());
	//M2PT
	TCanvas ce("ce", "ce");
	m2pt_gr->SetTitle("JPsi peak mean vs m2pt bin average; m2pt bin range average; GeV");
	m2pt_gr->SetMarkerColor(4);
	m2pt_gr->SetMarkerStyle(21);
	//m2pt_gr->GetYaxis()->SetRangeUser(3.05, 3.13);
	m2pt_gr->Draw("AP");
	ce.SaveAs((reso+"/mass_means/m2pt_means.png").c_str());
	//M2ETA
	TCanvas cg("cg", "cg");
	m2eta_gr->SetTitle("JPsi peak mean vs m2eta bin average; m2eta bin range average; GeV");
	m2eta_gr->SetMarkerColor(4);
	m2eta_gr->SetMarkerStyle(21);
	//m2eta_gr->GetYaxis()->SetRangeUser(3.08, 3.095);
	m2eta_gr->Draw("AP");
	cg.SaveAs((reso+"/mass_means/m2eta_means.png").c_str());

	//PLOT COMPLETE AND SCALED MASS PICTURES
	TCanvas cmass1("cmass1", "cmass1");
	//mass_pt_scaled->SetStats(false);
	mass_pt_scaled->GetXaxis()->SetTitle("GeV");
	mass_pt_scaled->GetYaxis()->SetTitle("Events");
	mass_pt_scaled->SetTitle("Mass from mmpt bins, scaled by mmeta");
	TF1 *gausexpo1 = new TF1("gausexpo1", "gaus", left_bound, right_bound);
	gausexpo1->SetParameter(1, real_mass);
	mass_pt_scaled->Fit("gausexpo1", "", "", left_bound, right_bound);
	mass_pt_scaled->GetXaxis()->SetRangeUser(left_bound, right_bound);
	mass_pt_scaled->Draw("SAME");
	gausexpo1->Draw("SAME");
	cmass1.SaveAs((reso+"/mass_pictures/THE PT SCALED MASS.png").c_str());

	TCanvas cmass3("cmass3", "cmass3");
	//mass_eta_scaled->SetStats(false);
	mass_eta_scaled->GetXaxis()->SetTitle("GeV");
	mass_eta_scaled->GetYaxis()->SetTitle("Events");
	mass_eta_scaled->SetTitle("Mass from mmeta bins, scaled by mmpt");
	TF1 *gausexpo3 = new TF1("gausexpo3", "gaus", left_bound, right_bound);
	gausexpo3->SetParameter(1, real_mass);
	mass_eta_scaled->Fit("gausexpo3", "", "", left_bound, right_bound);
	mass_eta_scaled->GetXaxis()->SetRangeUser(left_bound, right_bound);
	mass_eta_scaled->Draw("SAME");
	gausexpo3->Draw("SAME");
	cmass3.SaveAs((reso+"/mass_pictures/THE ETA SCALED MASS.png").c_str());

	TCanvas cmass2("cmass2", "cmass2");
	mass_total->GetXaxis()->SetTitle("GeV");
	mass_total->GetYaxis()->SetTitle("Events");
	mass_total->SetTitle("Total Unscaled Mass");
	TF1 *gausexpo2 = new TF1("gausexpo2", "gaus", left_bound, right_bound);
	gausexpo2->SetParameter(1, real_mass);
	mass_total->Fit("gausexpo2", "", "", left_bound, right_bound);
	mass_total->GetXaxis()->SetRangeUser(left_bound, right_bound);
	mass_total->Draw("SAME");
	gausexpo2->Draw("SAME");
	cmass2.SaveAs((reso+"/mass_pictures/UNSCALE MASS.png").c_str());
	
	//TO UPDATE SCALE FACTORS
	//for (const auto scale : pt_scales){ cout << scale << endl; }
	//for (const auto scale : eta_scales){ cout << scale << endl; }
	
	double pt_width = (double)(gausexpo1->GetParameter(2))*2.355;
	double eta_width = (double)(gausexpo3->GetParameter(2))*2.355;
	double total_width = (double)(gausexpo2->GetParameter(2))*2.355;
	double pt_width_change = 1000*sqrt(abs(pow(pt_width, 2)-pow(total_width, 2)));
	double eta_width_change = 1000*sqrt(abs(pow(eta_width, 2)-pow(total_width, 2)));
	string pt_change = "increase";
	string eta_change = "increase";
	if (pt_width < total_width) pt_change = "decrease";
	if (eta_width < total_width) eta_change = "decrease"; 

	cout << endl << "results for " + reso << endl << endl;
	cout << "nominal mass: " << real_mass << endl << endl;
	cout << "mass results" << endl;
	cout << "uncorrected mass: " << gausexpo2->GetParameter(1) << " +-" << gausexpo2->GetParError(1) << endl;
	cout << "eta corrected mass: " << gausexpo3->GetParameter(1) << " +-" << gausexpo3->GetParError(1) << endl;
	cout << "pt and eta corrected mass: " << gausexpo1->GetParameter(1) << " +-" << gausexpo1->GetParError(1) << endl << endl;
	cout << "width results" << endl;
	cout << "uncorrected width: " << total_width << " +-" << gausexpo2->GetParError(2) << endl;
	cout << "eta corrected width: " << eta_width << " +-" << gausexpo3->GetParError(2) << endl;
	cout << "change in width: " << eta_width_change << " MeV " << eta_change << endl;
	cout << "pt and eta corrected width: " << pt_width << " +-" << gausexpo1->GetParError(2) << endl;
	cout << "change in width: " << pt_width_change << " MeV " << pt_change << endl << endl;
}
