//#include "sumwgt.h"
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

void mctreereader(string treepath = "/mnt/hadoop/scratch/gandreas/NanoAOD/505/BuToJpsiK_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen+RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2+MINIAODSIM/MCmatch_flat_Bukmm/0-50.root"){

	///////////////////////////////
	/////VAR AND MASS ANALYSIS/////
	///////////////////////////////
	//FOR SIDEBAND SUBTRACTED PLOT
	//DIMUON
	//PT
	int bins = 13;
	TH1F* sigpt = new TH1F("sigpt", "signal region", bins, 0., 30.);
	TH1F* sbpt = new TH1F("sbpt", "sideband region", bins, 0., 30.);
	TH1F* sigsubpt;

	//MUON 1
	//ETA
	TH1F* m1sigeta = new TH1F("m1sigeta", "signal region", bins, -3.15, 3.15);
	TH1F* m1sbeta = new TH1F("m1sbeta", "sideband region", bins, -3.15, 3.15);
	TH1F* m1sigsubeta;

	//JPSI MASS FOR SCALING
	double jpsi = 3.0969;

	//LOAD PREVIOUSLY CREATED FILE WITH SIG SUB PLOTS OF VARS
	TFile* sub_vars = TFile::Open("mcsubs.root");
	TH1F* pthisto = (TH1F*) sub_vars->Get("sigsubpt");
	TH1F* m1etahisto = (TH1F*) sub_vars->Get("m1sigsubeta");

	//FILL WITH THE VARIABLE RANGE FOR EACH MASS PICTURE
	std::vector <double> pt_bin_avg;
	std::vector <double> m1eta_bin_avg;

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

	TChain* chain = new TChain("Candidates");
    	chain->Add((TString)treepath);
    	TTreeReader reader(chain);

	TTreeReaderValue<float> mass (reader, "mm_kin_mass");
	TTreeReaderValue<float> mmpt (reader, "mm_kin_pt");
	TTreeReaderValue<float> m1eta (reader, "mm_kin_mu1eta");

	while(reader.Next()){

		//FILL SIGNAL AND SIDEBAND DIMUON PT PLOTS
		double sig_lo = 2.95;
		double sig_hi = 3.2;
		double sb_l_lo = 2.6;
		double sb_l_hi = 2.85;
		double sb_r_lo = 3.3;
		double sb_r_hi = 3.55;
		if (*mass>=sig_lo && *mass<=sig_hi) {
			sigpt->Fill(*mmpt);
			m1sigeta->Fill(*m1eta);
		}
		if ( (*mass>=sb_l_lo && *mass<=sb_l_hi) || (*mass>=sb_r_lo && *mass<=sb_r_hi) ) {
			sbpt->Fill(*mmpt);
			m1sbeta->Fill(*m1eta);
		}

		//CREATE MASS PICTURE FOR EACH PT BIN
		double bin_lo;
		double bin_hi;
		for (int i = 1; i < bins; i++){
			bin_lo = pthisto->GetBinCenter(i);
			bin_hi = pthisto->GetBinCenter(i+1);
			pt_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (*mmpt>=bin_lo && *mmpt<bin_hi){
				pt_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" GeV pt").c_str());
				pt_mass_pics[i-1]->SetStats(false);
				pt_mass_pics[i-1]->Fill(*mass); 
				//if (eta_scales_saved[i-1] != 0.0) {
					//mass_pt_scaled->Fill(mass*eta_scales_saved[i-1]);
				//}
			}
		}
		for (int i = 1; i < bins; i++){
			bin_lo = m1etahisto->GetBinCenter(i);
			bin_hi = m1etahisto->GetBinCenter(i+1);
			m1eta_bin_avg.push_back((bin_lo+bin_hi)/2);
			if (*m1eta>=bin_lo && *m1eta<bin_hi){
				m1eta_mass_pics[i-1]->SetTitle(("Mass from "+std::to_string(bin_lo)+"-"+std::to_string(bin_hi)+" m1eta").c_str());
				m1eta_mass_pics[i-1]->SetStats(false);
				m1eta_mass_pics[i-1]->Fill(*mass);
				//if (pt_scales_saved[i-1] != 0.0) {
					//mass_eta_scaled->Fill(mass*pt_scales_saved[i-1]);
				//}
			}
		}



	}//While loop end



	//WRITE SIGNAL AND SIDEBAND PLOTS
	TFile *subs = new TFile("mcsubs.root", "RECREATE");
	
	//DIMUON
	sigsubpt = (TH1F*) sigpt->Clone("sigsubpt");
	sigsubpt->Add(sbpt, -.5);
	sigsubpt->SetTitle("signal subtracted pt");
	sigsubpt->Write();
	//MUON1
	m1sigsubeta = (TH1F*) m1sigeta->Clone("m1sigsubeta");
	m1sigsubeta->Add(m1sbeta, -.5);
	m1sigsubeta->SetTitle("signal subtracted eta muon 1");
	m1sigsubeta->Write();

	//PLOT THE SIG, SB, and SUB
	TCanvas c1("c1", "c1");
	sigsubpt->SetStats(false);
	sigsubpt->SetLineColor(kBlue);
	sigsubpt->GetXaxis()->SetTitle("GeV");
	sigsubpt->GetYaxis()->SetTitle("Events");
	sigsubpt->Draw("SAME");
	c1.SaveAs("sig_sub_pt_mc.png");

	TCanvas c6("c6", "c6");
	m1sigsubeta->SetStats(false);
	m1sigsubeta->SetLineColor(kBlue);
	m1sigsubeta->GetXaxis()->SetTitle("Eta");
	m1sigsubeta->GetYaxis()->SetTitle("Events");
	m1sigsubeta->Draw("SAME");
	c6.SaveAs("m1sig_sub_eta_mc.png");

	//WRITE THE MASS PICTURES FROM THE SIDEBAND PLOTS
	std::vector<double> ptmeans;
	std::vector<double> ptmeans_scaled;
	std::vector<double> pterr;
	std::vector<double> pt0;
	std::vector<double> pt_scales;
	int i = 0;
	for (const auto &hist : pt_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");
		hist->GetXaxis()->SetRangeUser(2.7, 3.5);
		TF1 *fit_func = new TF1("fit_func", "gaus", 2.7, 3.5);
		hist->Fit("fit_func", "", "", 2.7, 3.5);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs(("mass_pictures/mass_picture_ptbin"+std::to_string(i+1)+".png").c_str());

		double scale = jpsi/(fit_func->GetParameter(1));
		pt_scales.push_back(scale);

		ptmeans.push_back(fit_func->GetParameter(1));
		ptmeans_scaled.push_back(fit_func->GetParameter(1)*scale);

		pterr.push_back(fit_func->GetParError(1));
		pt0.push_back(0);
		i++;
	}
	std::vector<double> m1etameans;
	std::vector<double> m1eta_scaled;
	std::vector<double> m1etaerr;
	std::vector<double> m1eta0;
	std::vector<double> m1eta_scales;
	i = 0;
	for (const auto &hist : m1eta_mass_pics){
		TCanvas c("c", "c");
		hist->GetXaxis()->SetTitle("GeV");
		hist->GetYaxis()->SetTitle("Events");
		hist->GetXaxis()->SetRangeUser(2.7, 3.5);
		TF1 *fit_func = new TF1("fit_func", "gaus", 2.7, 3.5);
		hist->Fit("fit_func", "", "", 2.7, 3.5);

		hist->Draw("SAME");
		fit_func->Draw("SAME");

		c.SaveAs(("mass_pictures/mass_picture_m1etabin"+std::to_string(i+1)+".png").c_str());

		double scale = jpsi/fit_func->GetParameter(1);
		m1eta_scales.push_back(scale);

		m1etameans.push_back(fit_func->GetParameter(1));
		m1eta_scaled.push_back(fit_func->GetParameter(1)*scale);

		m1etaerr.push_back(fit_func->GetParError(1));
		m1eta0.push_back(0);
		i++;
	}

	//CREATE PLOTS FOR MEAN OF PEAK WITH CHANGES
	TLine* line = new TLine(-5, 3.0969, 30, 3.0969);

	TGraphErrors* pt_gr = new TGraphErrors(12, &pt_bin_avg[0], &ptmeans[0], &pt0[0], &pterr[0]);
	TGraphErrors* ptmeans_scaled_gr = new TGraphErrors(12, &pt_bin_avg[0], &ptmeans_scaled[0], &pt0[0], &pterr[0]);
	TGraphErrors* m1eta_gr = new TGraphErrors(12, &m1eta_bin_avg[0], &m1etameans[0], &m1eta0[0], &m1etaerr[0]);
	TGraphErrors* m1eta_scaled_gr = new TGraphErrors(12, &m1eta_bin_avg[0], &m1eta_scaled[0], &m1eta0[0], &m1etaerr[0]);

	TCanvas c2("c2", "c2");
	TLegend* leg1 = new TLegend();
	pt_gr->SetTitle("JPsi peak mean vs pt bin average; pt bin range average; GeV");
	pt_gr->SetMarkerColor(4);
	pt_gr->SetMarkerStyle(21);
	pt_gr->GetYaxis()->SetRangeUser(3.09, 3.1);
	pt_gr->Draw("AP");
	ptmeans_scaled_gr->SetMarkerColor(kRed);
	ptmeans_scaled_gr->SetMarkerStyle(21);
	ptmeans_scaled_gr->Draw("P");
	line->Draw("P");
	leg1->AddEntry(pt_gr, "uncorrected");
	leg1->AddEntry(ptmeans_scaled_gr, "corrected");
	leg1->Draw();
	c2.SaveAs("mass_means/pt_means.png");

	TCanvas c7("c7", "c7");
	TLegend* leg2 = new TLegend();
	m1eta_gr->SetTitle("JPsi peak mean vs m1eta bin average; m1eta bin range average; GeV");
	m1eta_gr->SetMarkerColor(4);
	m1eta_gr->SetMarkerStyle(21);
	m1eta_gr->GetYaxis()->SetRangeUser(3.09, 3.1);
	m1eta_gr->Draw("AP");
	m1eta_scaled_gr->SetMarkerColor(kRed);
	m1eta_scaled_gr->SetMarkerStyle(21);
	m1eta_scaled_gr->Draw("P");
	line->Draw("P");
	leg2->AddEntry(m1eta_gr, "uncorrected");
	leg2->AddEntry(m1eta_scaled_gr, "correct");
	leg2->Draw();
	c7.SaveAs("mass_means/m1eta_corrected_means.png");

}
