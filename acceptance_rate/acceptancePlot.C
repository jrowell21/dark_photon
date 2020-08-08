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

void acceptancePlot(){
	
	std::vector<string> files;
	files.push_back("m0p5");
	files.push_back("m1");
	files.push_back("m2");
	files.push_back("m3");
	files.push_back("m4");
	files.push_back("m5");
	files.push_back("m6");
	files.push_back("m7");
	files.push_back("m8");
	files.push_back("m9");
	files.push_back("m10");

	TFile* file;
	TDirectoryFile* dir;

	TH1D* mass1;
	TH1D* mass2;

	double entries1;
	double entries2;
	double acceptance;
	double error;

	std::vector<double> acceptance_rates;
	std::vector<double> acceptance_errors;
	std::vector<double> acceptance_os;

	for (const auto fileName:files){
		file = TFile::Open(("DarkPhotonAcceptance"+fileName+".root").c_str());
		dir = (TDirectoryFile*) file->Get("TestHepMCEvt");
		mass1 = (TH1D*) dir->Get("Hist2muMass");
		mass2 = (TH1D*) dir->Get("Hist2muMassFid");

		entries1 = mass1->GetEntries();
		entries2 = mass2->GetEntries();

		acceptance = entries2/entries1;
		error = acceptance*sqrt( pow( sqrt(entries1) / entries1, 2) + pow(sqrt(entries2) / entries2, 2));

		acceptance_rates.push_back(acceptance);
		acceptance_errors.push_back(error);
		acceptance_os.push_back(0.0);
	}

	std::vector<double> mass;
	mass.push_back(0.5);
	mass.push_back(1.0);
	mass.push_back(2.0);
	mass.push_back(3.0);
	mass.push_back(4.0);
	mass.push_back(5.0);
	mass.push_back(6.0);
	mass.push_back(7.0);
	mass.push_back(8.0);
	mass.push_back(9.0);
	mass.push_back(10.0);
	
	int n = mass.size();
	TGraphErrors* gr = new TGraphErrors(n, &mass[0], &acceptance_rates[0], &acceptance_os[0], &acceptance_errors[0]);

	TF1* fit_func = new TF1("fit_func", "pol3", 0, 11);

	gr->Fit("fit_func", "", "", 0, 11);
	
	TCanvas c("c", "c");
	gr->SetTitle("Acceptance Rate vs Mass; GeV; Acceptance Rate");
	gr->SetMarkerColor(4);
	gr->SetMarkerStyle(21);
	gr->Draw("AP");
	fit_func->Draw("SAME");
	c.SaveAs("Acceptance.png");

	TFile* function_out = new TFile("acceptance_fit_function.root", "RECREATE");
	fit_func->Write();
	function_out->Close();

}
