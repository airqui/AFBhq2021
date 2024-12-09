
#include <TPaveStats.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <TFitResult.h>
#include <TF1.h>
#include <TSpectrum.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystemFile.h"
#include "../../../style/Style.C"
#include "../../../style/Labels.C"
#include "../hidden_cross_sections.h"
#include "histograms.h"

void plotsReco_1D( int cuts = 0, float lum = 900)
{

  folder = TString::Format("../results_FCC/selection_cuts%i", cuts);

  /*  //-------------------------------
  // 1d histograms to be plotted.. they have to exist in the selection_XXX.root file
  std::vector<TString> histonames = {"h_mjj"};//, "h_mj1_mj2", "h_y23", "h_d23", "h_thrust", "h_major_thrust", "h_minor_thrust","h_sphericity"};
  //labels

  TString histo1d_titles[] = {
    "",
    "M_{j_{1}j_{2}} [GeV]"};*/
  
  std::vector<TString> histonames = {"h_mjj", "h_mj1_mj2", "h_y23", "h_d23", "h_thrust", "h_major_thrust", "h_minor_thrust","h_sphericity"};
  //labels
  TString histo1d_titles[] = {
      "",
      "M_{j_{1}j_{2}} [GeV]",
      "M_{j_{1}}+M_{j_{2}} [GeV]",
      "y23",
      "d23 [GeV]",
      "T-principle",
      "T-major",
      "T-minor",
      "sphericity"};


  std::vector<std::vector<TH1F *>> h1_bkg;

  for (int isample = 0; isample < sizeof(samples2) / sizeof(TString); isample++)
    {

      std::vector<TH1F *> h1_bkg_temp = GetHisto1D(samples2[isample], lum, histonames);
      h1_bkg.push_back(h1_bkg_temp);
    }

  for (int i = 0; i < h1_bkg.size(); i++)
    {
      cout << samples2[i] << ": " << h1_bkg.at(i).at(1)->Integral() << " ";
    }
  cout << endl;

  SetQQbarStyle();
  TGaxis::SetMaxDigits(3);

  for (int k = 1; k < sizeof(histo1d_titles) / sizeof(TString); k++)
    {


      float xmin = 0.4, ymin = 0.7, xmax = 0.8, ymax = 0.9;

      TLegend *leg = new TLegend(xmin, 0.5, xmax, 0.9); //(0.4,0.3,0.5,0.6);
      leg->SetTextSize(0.035);

      TCanvas *canvas1 = new TCanvas(TString::Format("canvas_%i", k), TString::Format("canvas_%i", k), 800, 800);
      canvas1->cd(1);
      h1_bkg.at(0).at(k)->GetYaxis()->SetTitle("Entries");
      h1_bkg.at(0).at(k)->GetXaxis()->SetTitleOffset(1.2);
      h1_bkg.at(0).at(k)->GetXaxis()->SetRangeUser(0,150);;
      h1_bkg.at(0).at(k)->GetXaxis()->SetTitle(histo1d_titles[k]);

      for (int j = 0; j < h1_bkg.size(); j++)
	{

	  if (j == 0)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kBlue-2);
	      h1_bkg.at(j).at(k)->SetLineWidth(6);
	    }
	  if (j == 1)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kBlue);
	      h1_bkg.at(j).at(k)->SetLineWidth(4);
	    }
	  if (j == 2)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kBlue+2);
	      h1_bkg.at(j).at(k)->SetLineWidth(2);
	      h1_bkg.at(j).at(k)->SetLineStyle(1);
	    }
	  if (j == 3)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kBlue+2);
	      h1_bkg.at(j).at(k)->SetLineWidth(3);
	      h1_bkg.at(j).at(k)->SetLineStyle(2);
	    }

	  if (j == 4)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kGreen+2);
	      h1_bkg.at(j).at(k)->SetLineWidth(4);
	      h1_bkg.at(j).at(k)->SetLineStyle(1);
	    }
	  if (j == 5)
	    {
	      h1_bkg.at(j).at(k)->SetLineColor(kGreen-2);
	      h1_bkg.at(j).at(k)->SetLineWidth(3);
	      h1_bkg.at(j).at(k)->SetLineStyle(2);
	    }
		
		
	  if (j == h1_bkg.size()-1)
	    {
	      h1_bkg.at(j).at(k)->SetFillColor(kRed);
	      h1_bkg.at(j).at(k)->SetLineColor(kRed);
	      h1_bkg.at(j).at(k)->SetLineWidth(2);
	      h1_bkg.at(j).at(k)->SetFillStyle(3003);
	    }
	  
	  // h1_bkg.at(j).at(k)->Draw("histosame");

	  if (j<h1_bkg.size()-1) leg->AddEntry(h1_bkg.at(j).at(k), "#font[42]{" + title_samples2[j] + "}", "l");
	  else leg->AddEntry(h1_bkg.at(j).at(k), "#font[42]{" + title_samples2[j] + "}", "f");
	}

      h1_bkg.at(h1_bkg.size()-1).at(k)->Draw("histo");
      for (int j = 0; j < h1_bkg.size()-1; j++)    h1_bkg.at(j).at(k)->Draw("histosame");

      leg->SetFillStyle(0);
      leg->SetLineWidth(0);
      leg->SetLineColor(0);
      leg->SetBorderSize(0);

      LabelsReco();
      leg->Draw();
      canvas1->Print(TString::Format("plots/FCC_vs_ILC_cut%i_%s.eps",cuts,histonames[k-1].Data()));
    }




}


void plotsReco_2D( int cuts = 0, float lum = 900)
{

  folder = TString::Format("../results_FCC/selection_cuts%i", cuts);


  //-------------------------------
  // 2d histograms to be plotted
  std::vector<TString> histonames_2d = {"h_nch", "h_npfos", "h_costheta_energy", "h_mjjmin_mjjmax", "h_major_minor_thrust","h_ChargedMomDiff_ChargedMomSum","h_Emin_Emax","h_nvtx1_nvtx2"};

  //labels to be added to the 2d plots
  TString histo2d_titles_x[] = {
				"# tracks j_{1}",
      "# pfos j_{1}",
      "|cos #theta| most energetic #gamma_{cand}",
      "m_{jj}^{min} [GeV] (4-jet reco)",
      "T-major",
      "MomCharged / MomAll (jet1)",
      "Min. Jet Energy",
      "nvtx1"
      };

  TString histo2d_titles_y[] = {
				"# tracks j_{2}",
      "# pfos j_{2}",
      "E most energetic #gamma_{cand} [GeV]",
      "m_{jj}^{max} [GeV] (4-jet reco)",
      "T-minor",
      "MomCharged / MomAll (jet2)",
      "Max. Jet Energy",
      "nvtx2"
      };

   std::vector<TString> histonames = {"h_mjj", "h_mj1_mj2", "h_y23", "h_d23", "h_thrust", "h_major_thrust", "h_minor_thrust","h_sphericity"};
  //labels
  TString histo1d_titles[] = {
      "",
      "M_{j_{1}j_{2}} [GeV]",
      "M_{j_{1}}+M_{j_{2}} [GeV]",
      "y23",
      "d23 [GeV]",
      "T-principle",
      "T-major",
      "T-minor",
      "sphericity"};


  std::vector<std::vector<TH1F *>> h1_bkg;
  std::vector<std::vector<TH2F *>> h2_bkg;

  for (int isample = 0; isample < sizeof(samples2) / sizeof(TString); isample++)
  {
    std::vector<TH2F *> h2_bkg_temp = GetHisto2D(samples2[isample], lum, histonames_2d,false);
    h2_bkg.push_back(h2_bkg_temp);
    std::vector<TH1F *> h1_bkg_temp = GetHisto1D(samples2[isample], lum, histonames);
    h1_bkg.push_back(h1_bkg_temp);
  }

  SetQQbarStyle();
  TGaxis::SetMaxDigits(3);

 
  for (int k = 0; k < sizeof(histo2d_titles_y) / sizeof(TString); k++)
  {

    gStyle->SetPadRightMargin(0.2);
    TCanvas *canvas1 = new TCanvas(TString::Format("canvas2d_%i", k), TString::Format("canvas2d_%i", k), 1200, 600);
    canvas1->Divide(2, 1);
   
    for (int j = 0; j < h2_bkg.size(); j++)
    {
      canvas1->cd(j + 1);
      gPad->SetLogz();
      h2_bkg.at(j).at(k)->GetXaxis()->SetTitle(histo2d_titles_x[k]);
      h2_bkg.at(j).at(k)->GetYaxis()->SetTitle(histo2d_titles_y[k]);
      h2_bkg.at(j).at(k)->Draw("colz");
      QQBARLabel2(0.2, 0.85, TString::Format("#font[42]{%s, N_{total}=%i}", samples2[j].Data(), int(h1_bkg.at(j).at(1)->Integral())), kRed,0.04);

      Labels();
      canvas1->Print(TString::Format("plots/FCC_vs_ILC_cut%i_2D_%s.eps",cuts,histonames_2d[k].Data()));
    }
  }
}

void selection_plots_ISR_FCC_vs_ILC()
{

  float lum = 100;
  int pol = 0;
  cout << "Events for Polarization " << pol << " (0=left, 1=right, 2=80left,30right, 3=80right,30left) and Lum=" << lum << endl;
  for (int cuts = 1; cuts < 2; cuts++)
    {
      cout << cuts << " ";
      //plotsReco_1D(cuts,lum);
      plotsReco_2D(cuts,lum);
      // plotsProcLCWS2023(cuts);
    
    }
}
