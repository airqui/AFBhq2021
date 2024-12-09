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

void plotsReco( int cuts = 0, float lum = 900)
{

  folder = TString::Format("../results_cut6_npfo22_nchpfo15/selection_cuts%i", cuts);

  //-------------------------------
  // 1d histograms to be plotted.. they have to exist in the selection_XXX.root file
  std::vector<TString> histonames = {"h_mjj"};//, "h_mj1_mj2", "h_y23", "h_d23", "h_thrust", "h_major_thrust", "h_minor_thrust","h_sphericity"};
  //labels
  TString histo1d_titles[] = {
    "",
    "M_{j_{1}j_{2}} [GeV]"};
  
  std::vector<std::vector<TH1F *>> h1_bkg;

  for (int isample = 0; isample < sizeof(samples) / sizeof(TString)-1; isample++)
    {

      std::vector<TH1F *> h1_bkg_temp = GetHisto1D(samples[isample], lum, histonames);
      h1_bkg.push_back(h1_bkg_temp);
    }

  for (int i = 0; i < h1_bkg.size(); i++)
    {
      cout << samples[i] << ": " << h1_bkg.at(i).at(1)->Integral() << " ";
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
	  
	  h1_bkg.at(j).at(k)->Draw("histosame");

	  if (j<h1_bkg.size()-1) leg->AddEntry(h1_bkg.at(j).at(k), "#font[42]{" + title_samples[j] + "}", "l");
	  else leg->AddEntry(h1_bkg.at(j).at(k), "#font[42]{" + title_samples[j] + "}", "f");
	}
      leg->SetFillStyle(0);
      leg->SetLineWidth(0);
      leg->SetLineColor(0);
      leg->SetBorderSize(0);

      LabelsReco();
      leg->Draw();
      canvas1->Print(TString::Format("plots/%s.eps",histonames[k-1].Data()));
    }




}    
void selection_plots()
{

  float lum = 100;
  int pol = 0;
  cout << "Events for Polarization " << pol << " (0=left, 1=right, 2=80left,30right, 3=80right,30left) and Lum=" << lum << endl;
  for (int cuts = 0; cuts < 7; cuts++)
    {
      if(cuts==1) cuts=6;
      cout << cuts << " ";
      plotsReco(cuts,lum);
      // plotsProcLCWS2023(cuts);
    
    }
}
