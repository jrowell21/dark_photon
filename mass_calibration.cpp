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
#include "pdfs.h"

using namespace std;

void mass_calibration(TString year = "2017"){

  //gROOT->ProcessLine(".L ../RooFit-pdfs/CB.cc+");

  	TFile* file;
  	if (year.CompareTo("2017")) file=TFile::Open("/mnt/hadoop/scratch/gandreas/DP_histos/2017_ScoutingRunD_4fb.root");
  	else if (year.CompareTo("2018")) file=TFile::Open("/mnt/hadoop/scratch/gandreas/DP_histos/2018_ScoutingRunC_6.6fb.root");
	
	std::map<string, RooDataHist*> hist_map; //map resonance name to RooDataHist
	std::map<string, RooAddPdf*> pdf_map; //map resonance name to pdf

	std::map<string, vector<float>> hist_ranges;
	hist_ranges["Eta"] = {{0.4,0.65}};
	hist_ranges["OmegaPhi"] = {{0.65,1.3}};
	hist_ranges["JPsiPsi"] =  {{2.6, 4.}};
	hist_ranges["Upsilon"] =  {{8.5 ,11}};

	std::vector <float> mass, massError, rel_width, rel_widthError, abs_width, abs_widthError, sig_integral, sb_integral;
	std::vector <string> resonances;
	resonances.push_back("eta");
	resonances.push_back("omega");
	resonances.push_back("phi");
	resonances.push_back("jpsi");
	resonances.push_back("psi2s");
	resonances.push_back("upsilon1s");
	resonances.push_back("upsilon2s");
	resonances.push_back("upsilon3s");
	
	//Creates PDFs inside the workspace
	reso_pdfs pdfs;
	RooWorkspace *w = pdfs.w;

   //LOOP OVER MASS INDICES AND MAKE THE CARDS/WORKSPACES
	TH1F* histo=(TH1F*)file->Get("massforLimitFull");

	RooArgSet chi2s;
	std::vector<Double_t> chis;
	std::vector<Double_t> probs;
	double sub_ratio, sig, sub, tot_int1;
	for (const auto &name_range : hist_ranges){
		string name = name_range.first;

	  	float xmin = name_range.second[0];
	  	float xmax = name_range.second[1];

	  	string name_lower = name;
	  	std::for_each(name_lower.begin(), name_lower.end(), [](char & c) {c = ::tolower(c);});

		//Gets the observable from the workspace and sets the fit range
	  	RooRealVar* m2mu = w->var(("m2mu_"+name_lower).c_str());
	  	m2mu->setRange(xmin, xmax);

		//Creates a histogram from the observable and the data and gets the fit model from workspace
	  	RooDataHist* rdh = new RooDataHist(("rdh_"+name_lower).c_str(), ("rdh_"+name_lower).c_str(), RooArgSet(*m2mu), histo);
	  	RooAddPdf* model = (RooAddPdf*) w->pdf(("model_"+name_lower).c_str());
	  	hist_map[name_lower] = rdh;
	  	pdf_map[name_lower] = model;

		//Chi2 values for minimizing, also probability
	  	RooAbsReal* chi2 = model->createChi2(*rdh, RooFit::Range(("fitRange_"+name_lower).c_str()));
		Double_t p = TMath::Prob(chi2->getVal(), rdh->numEntries()-1);
	  	chi2s.add(*chi2);
		chis.push_back(chi2->getVal());
		probs.push_back(p);
	 }
	/*//Saves the chi2s and probs
	Int_t probs_size = probs.size();
	TGraph* chiprob = new TGraph(probs_size, &chis[0], &probs[0]);
	TFile* chiFile = new TFile(("chi2s8order_"+(string)year+".root").c_str(), "RECREATE");
	chiprob->Write();
	chiFile->Close();*/
	
	//Minimizes results of fit
	RooAddition totchi2 ("totchi2", "totchi2", chi2s);
	RooMinuit m(totchi2);
	m.migrad();
	m.hesse();

	//Saves the results of fit into file
	RooFitResult* Ares = new RooFitResult();
	Ares = m.save();
	cout<<"$$$$"<<Ares->covQual()<<endl;
	TFile *fitResults = new TFile(("fullRangeFitResults"+(string)year+".root").c_str(), "RECREATE");
	Ares->Write();
	fitResults->Close();

	for (const auto &name_range : hist_ranges){
		string name = name_range.first;
	  	string name_lower = name;
	  	std::for_each(name_lower.begin(), name_lower.end(), [](char & c) {c = ::tolower(c);});
		
		//Creates a frame for observable
	  	RooPlot *frame = w->var(("m2mu_"+name_lower).c_str())->frame(RooFit::Range(("fitRange_"+name_lower).c_str()));

		//Plots the histogram data
		hist_map[name_lower]->plotOn(frame);

		//Iterator to plot all signals over the range
		TIterator* pdfIter = pdf_map[name_lower]->pdfList().createIterator();
		RooAbsPdf* signal = (RooAbsPdf*)pdfIter->Next();
		while (signal){
			pdf_map[name_lower]->plotOn(frame, RooFit::NormRange(("fitRange_"+name_lower).c_str()), RooFit::Components(signal->GetName()), RooFit::LineColor(kGreen),RooFit::Name("s")); //Plots the signal
			signal = (RooAbsPdf*)pdfIter->Next();
		}

		//Plots the background and signal+background
		pdf_map[name_lower]->plotOn(frame, RooFit::NormRange(("fitRange_"+name_lower).c_str()), RooFit::Components(("bkgModel_"+name_lower).c_str()), RooFit::LineColor(kRed),RooFit::Name("b")); //Plots the background
		pdf_map[name_lower]->plotOn(frame, RooFit::NormRange(("fitRange_"+name_lower).c_str()),RooFit::Name("s+b")); //Plots the signal + background
		//pdf_map[name_lower]->plotOn(frame, RooFit::NormRange(("fitRange_"+name_lower).c_str()), RooFit::Components(("signalModel_gau_"+name_lower).c_str()), RooFit::LineColor(kOrange));//Trying to see the crystal ball and gaus individually, but need resonance names not hist_ranges
		//pdf_map[name_lower]->plotOn(frame, RooFit::NormRange(("fitRange_"+name_lower).c_str()), RooFit::Components(("signalModel_CB_"+name_lower).c_str()), RooFit::LineColor(kPink));//Need resonance names not hist_ranges
		
		//Resid Plot Top and Bottom
		RooPlot *frame_top = w->var(("m2mu_"+name_lower).c_str())->frame(RooFit::Range(("fitRange_"+name_lower).c_str()));
		hist_map[name_lower]->plotOn(frame_top);
		pdf_map[name_lower]->plotOn(frame_top, RooFit::NormRange(("fitRange_"+name_lower).c_str()));
		
		RooHist* hresid = frame_top->residHist();
		RooPlot* frame_bottom = w->var(("m2mu_"+name_lower).c_str())->frame(RooFit::Range(("fitRange_"+name_lower).c_str()));
		frame_bottom->addPlotable(hresid, "P");

		//Creates a legend
		TLegend* leg = new TLegend();
		leg->SetHeader((name).c_str(),"C");
		leg->AddEntry(frame->findObject("s"),(name + " Signal").c_str(), "l");
		leg->AddEntry(frame->findObject("b"), "Background","l");
		leg->AddEntry(frame->findObject("s+b"),"Signal + Background","l");

		//Plots the signal+background
		TCanvas c_resid("c_all", "c_all", 800, 500);
		c_resid.Divide(0,2);
		c_resid.cd(1);
		frame_top->Draw();
		frame_top->SetTitle((name+(string)year).c_str());
		frame_top->GetXaxis()->SetTitle("dimuon mass [GeV]");
		
		//Plots the residuals
		c_resid.cd(2);
		frame_bottom->Draw();
		frame_bottom->SetTitle("Residuals");
		frame_bottom->GetXaxis()->SetTitle("dimuon mass [GeV]");
		c_resid.SaveAs(("plots"+(string)year+"/"+name_lower+(string)year+"resids.png").c_str());
		
		//Plots the fit
		TCanvas c_all("c_all", "c_all", 800, 500);
		leg->Draw();
		frame->Draw();
		frame->SetTitle((name+(string)year).c_str());
		frame->GetXaxis()->SetTitle("dimuon mass [GeV]");
		c_all.SaveAs(("plots"+(string)year+"/"+name_lower+(string)year+".png").c_str());
	} 
	//Adds mass, rel/abs widths with errors to respective vectors
	for (auto reso : resonances){
		mass.push_back(w->var(("M1_"+reso).c_str())->getVal());
		massError.push_back(w->var(("M1_"+reso).c_str())->getError());
		rel_width.push_back(w->var(("res_rel1_"+reso).c_str())->getVal());
		rel_widthError.push_back(w->var(("res_rel1_"+reso).c_str())->getError());
		double abs_w = (w->var(("M1_"+reso).c_str())->getVal())*(w->var(("res_rel1_"+reso).c_str())->getVal());
		abs_width.push_back(abs_w);
		abs_widthError.push_back(abs_w*sqrt(((w->var(("M1_"+reso).c_str())->getError())/(w->var(("M1_"+reso).c_str())->getVal()))+((w->var(("res_rel1_"+reso).c_str())->getError())/(w->var(("res_rel1_"+reso).c_str())->getVal()))));
		//uncertainty added in quadrature ^^
	}

	//Graphs+Errors for rel/abs widths
	const Int_t n = resonances.size();
	TGraphErrors* gr_rel = new TGraphErrors(n,&mass[0],&rel_width[0],&massError[0],&rel_widthError[0]);
	TGraphErrors* gr_abs = new TGraphErrors(n,&mass[0],&abs_width[0],&massError[0],&abs_widthError[0]);

	//Plots absolute width vs mass
	TCanvas* c1 = new TCanvas("c1","c1",200,10,700,500);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor(21);
	c1->GetFrame()->SetBorderSize(12);
	gr_abs->SetName("abs_resos");
	gr_abs->SetTitle("Absolute resonance width vs mass; Mass GeV; Absolute Res. Width");
	gr_abs->SetMarkerColor(4);
	gr_abs->SetMarkerStyle(21);
	gr_abs->GetYaxis()->SetRangeUser(0,0.2);//(0,0.05)
	gr_abs->Draw("AP");
	c1->SaveAs(("plots"+(string)year+"/abs_res_graph"+(string)year+".png").c_str());
	
	//Plots relative width vs mass
	TCanvas* c2 = new TCanvas("c2", "c2", 200, 10, 700, 500);
	c2->SetGrid();
	c2->GetFrame()->SetFillColor(21);
	c2->GetFrame()->SetBorderSize(12);
	gr_rel->SetName("rel_resos");
	gr_rel->SetTitle("Relative resonance width vs mass; Mass GeV; Relative Res. Width");
	gr_rel->SetMarkerColor(4);
	gr_rel->SetMarkerStyle(21);
	gr_rel->GetYaxis()->SetRangeUser(0,0.05);//(0,0.05)
	gr_rel->Draw("AP");
	c2->SaveAs(("plots"+(string)year+"/rel_res_graph"+(string)year+".png").c_str());

	//Writes the mass and width data to file
	TFile *outf = new TFile(("mass_resolutions"+(string)year+".root").c_str(), "RECREATE");
	gr_rel->Write();
	gr_abs->Write();
	outf->Close();

	//Save workspace to disk
	w->writeToFile(("mass_cal_workspace"+(string)year+".root").c_str());
}
