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


int nbkgs=1;
//if we study ISR & VBF --> nbkgs=2



void plots_reco_C2Yield_PRL(float lum=2000, bool bkg = false, bool yield=false,TString folder_0="../results_new/", int cuts=7)
{

  std::vector<std::vector<TH1F *>> h1_signal;
  std::vector<std::vector<TH2F *>> h2_signal;
  std::vector<std::vector<TH1F *>> h1_bkg;
  std::vector<std::vector<TH2F *>> h2_bkg;

  std::vector<TString> histos_1d = {"ncharged"};
  std::vector<TString> histos_2d = {"S2_rapidity", "B2_rapidity"};



  //samples has to be the HV100,50,10,ISR
  for (int isample = 0; isample < sizeof(samples) / sizeof(TString) -1; isample++)
    {
    cout<<isample<<endl;
    //folder=TString::Format("%s/selection_cuts%i",folder_0.Data(),cuts);
    folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);

    //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
    std::vector<TH1F *> h1_temp = GetHisto1D(samples[isample], lum, histos_1d);
    h1_signal.push_back(h1_temp);
    //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
    std::vector<TH2F *> h2_temp = GetHisto2D(samples[isample], lum, histos_2d,true);
    for(int j=0; j<h2_temp.size();j++) {
      // h2_temp.at(j)->Rebin2D(2,4); //Rebin to help with short statistics
      if(yield==false) h2_temp.at(j)->Rebin2D(2,4); 
    }
    h2_signal.push_back(h2_temp);
  }


  for (int isample = sizeof(samples) / sizeof(TString); isample < sizeof(samples) / sizeof(TString)+1; isample++)
  {
    folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),20);
    std::vector<TH1F *> h1_temp = GetHisto1D(samples[isample-1], lum, histos_1d);
    h1_bkg.push_back(h1_temp);
    
    std::vector<TH2F *> h2_temp = GetHisto2D(samples[isample-1], lum, histos_2d,true);
    for(int j=0; j<h2_temp.size();j++) {
      //h2_temp.at(j)->Rebin2D(2,4);
      if(yield==false) h2_temp.at(j)->Rebin2D(2,4);
    }
    h2_bkg.push_back(h2_temp);
  }



  cout<<" EVENTS : ";
  for (int i = 0; i < h1_signal.size(); i++)
  {
    cout << samples[i] << ": " << h1_signal.at(i).at(1)->Integral() << " ";
  }
  cout << endl;

  SetQQbarStyle();
  gStyle->SetOptTitle(1);
  gStyle->SetTitleStyle(0);
  // gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  // gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleY(0.9);
  // gStyle->SetMarkerSize(1.5);
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(kRainBow);


  	
  for (int i = 0; i < h2_signal.size(); i++)    {

    cout<<" NUMERO DE EVENTOS: "<<samples[i]<<" ";
    float normaliz= h1_signal.at(h2_signal.size()-nbkgs).at(1)->Integral();
    cout<<normaliz<<" ";
    if(i!=(h2_signal.size()-nbkgs)) {
      normaliz+=h1_signal.at(i).at(1)->Integral();
      cout<<" + SM BKG: "<<normaliz;
    }
   
    cout<<endl;
   if(yield==false) h2_signal.at(i).at(0)->Divide(h2_bkg.at(0).at(1));    
  }

  if(yield==false) {  
    for (int j = 0; j < h2_signal.size(); j++)
      {


	TCanvas *canvas1 = new TCanvas(TString::Format("Correlations_cut%i_%s_reco", cuts,samples[j].Data()),TString::Format("Correlations_cut%i_%s_reco", cuts,samples[j].Data()), 800, 800);
	canvas1->cd();

	h2_signal.at(j).at(0)->Draw("surf1");
	h2_signal.at(j).at(0)->GetZaxis()->SetRangeUser(0,10);
	h2_signal.at(j).at(0)->GetXaxis()->SetTitleOffset(1.5);
	h2_signal.at(j).at(0)->GetYaxis()->SetTitleOffset(1.5);
	h2_signal.at(j).at(0)->GetZaxis()->SetTitleOffset(1.5);
	h2_signal.at(j).at(0)->GetZaxis()->SetTitle("C^{(2)}(#Deltay, #Delta#phi)");
	h2_signal.at(j).at(0)->GetYaxis()->SetTitle("#Delta#phi");
	h2_signal.at(j).at(0)->GetXaxis()->SetTitle("#Deltay");
	h2_signal.at(j).at(0)->GetXaxis()->SetRangeUser(-5,5);
	
	// h2_signal.at(j).at(0)->Draw("colz");
	
	QQBARLabel2(0.57, 0.05, TString::Format("#font[42]{%s}", title_samples[j].Data()), kRed + 2,0.055);
	
	//LabelsReco();
	if(bkg==true) canvas1->Print(TString::Format("Correlations_cut%i_%s_reco.eps", cuts,samples[j].Data()));
	else canvas1->Print(TString::Format("Correlations_cut%i_%s_reco_onlyHV.eps", cuts,samples[j].Data()));
      }
  } else
  {

    for(int step=0; step<2; step++) {

    TH1F *hyield[10];
    //TH1F *ratio_hyield[10];
    TGraph *g_ratio_yield[10];

    TCanvas *canvas1 = new TCanvas(TString::Format("yieldcanvas_range%i",step),TString::Format("yieldcanvas_range%i",step), 800, 800);
    canvas1->cd();
    
    TLegend *leg = new TLegend(0.4, 0.65, 0.75, 0.92); //(0.4,0.3,0.5,0.6);
    leg->SetTextSize(0.03);


    TF1 * minus[10];
    TF1 * plus[10];

    for (int j = 0; j < h2_signal.size(); j++)
      {

	minus[j] = new TF1(TString::Format("minus_j%i",j),"pol8",-TMath::Pi() / 2.,0);
	plus[j] = new TF1(TString::Format("plus_j%i",j),"pol8",0, 3.*TMath::Pi() / 2);


	int first=50;
	int bins=16;

	TH1F* htemp1;
	TH1F* htemp1_SM1;
	TH1F* htemp1_SM2;

        TH1F* htemp2;
        TH1F* htemp3;
	
	if(step==1) {
	  hyield[j] = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
	  htemp1 = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp1_%i",j), 1,first-step*bins);//, first-2*step*bins*2);
	  hyield[j]->Add(htemp1);

	} else {
	  hyield[j] = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first-bins,first+bins);//, first+2*step*bins*2);
	}

	if(j != (h2_signal.size()-1)) {
	  float na=h1_signal.at(j).at(1)->Integral();
	  float nb=h1_signal.at(h2_signal.size()-1).at(1)->Integral();
	  hyield[j]->Scale(na/(na+nb));

	  if(step==1) {
	    htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first+step*bins,first*2);//, first+2*step*bins*2);
	    htemp1_SM2 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp1_%i_SM2",j),1, first-step*bins);//, first-2*step*bins*2);
	    htemp1_SM1->Add(htemp1_SM2);
	  } else {
	    htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first-bins,first+bins);//, first+2*step*bins*2);
	  }
	  htemp1_SM1->Scale(nb/(na+nb));
	  hyield[j]->Add(htemp1_SM1);
	}
	
	htemp2 = (TH1F *)h2_bkg.at(0).at(1)->ProjectionY(TString::Format("htemp2_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
	htemp3 = (TH1F *)h2_bkg.at(0).at(1)->ProjectionY(TString::Format("htemp3_%i",j), 1,first-step*bins);//,first-2*step*bins*2);
	htemp2->Add(htemp3);

	
	for(int i=0; i<hyield[j]->GetNbinsX()+1; i++ ) {
	  if(htemp2->GetBinContent(i)!=0 && hyield[j]->GetBinContent(i)!=0 ) {
	    double a=hyield[j]->GetBinContent(i);
	    double b=htemp2->GetBinContent(i);
	    double ea=hyield[j]->GetBinError(i);
	    double eb=htemp2->GetBinError(i);
	    hyield[j]->SetBinError(i,(1./b)*sqrt(pow(ea,2)+pow(a*eb/b,2)));
	    hyield[j]->SetBinContent(i,a/b);//hyield[j]->GetBinContent(i)/htemp2->GetBinContent(i));
	    
	  } else {
	    hyield[j]->SetBinContent(i,0);
	    hyield[j]->SetBinError(i,0);
	  }
	}

	///	if(j==h2_signal.size()-1) {
	
	//TF1 * minus = new TF1("minus","pol8",-TMath::Pi() / 2.,0);
	hyield[j]->Fit(TString::Format("minus_j%i",j),"REMN");
	  
	  for(int i=0; i<hyield[j]->GetNbinsX()+1; i++ ) {
	    if(hyield[j]->GetBinCenter(i)<=0 ) {

	      float params[20]={0};
	      for(int ipol=0; ipol<9; ipol++) params[ipol]=minus[j]->GetParameter(ipol);
	      float value = params[0];
	      for(int ipol=1; ipol<9; ipol++) {
		value+=params[ipol]*pow(hyield[j]->GetBinCenter(i),ipol);
	      }
	      if(cuts==6 && j==6) hyield[j]->SetBinContent(i,value);
	    }
	  }

	  hyield[j]->Fit(TString::Format("plus_j%i",j),"REMN");

	  for(int i=0; i<hyield[j]->GetNbinsX()+1; i++ ) {
	    if(hyield[j]->GetBinCenter(i)>0 ) {

	      float params[20]={0};
	      for(int ipol=0; ipol<9; ipol++) params[ipol]=plus[j]->GetParameter(ipol);
	      float value = params[0];
	      for(int ipol=1; ipol<9; ipol++) {
		value+=params[ipol]*pow(hyield[j]->GetBinCenter(i),ipol);
	      }
	      if(cuts==6 && j==6) hyield[j]->SetBinContent(i,value);
	    }
	  }
//	}
	
	  if (j == 0)
        {
          hyield[j]->SetLineColor(kBlue-2);
          hyield[j]->SetLineWidth(6);
        }
        if (j == 1)
        {
          hyield[j]->SetLineColor(kBlue);
          hyield[j]->SetLineWidth(4);
        }
	if (j == 2)
        {
          hyield[j]->SetLineColor(kBlue+2);
          hyield[j]->SetLineWidth(2);
          hyield[j]->SetLineStyle(1);
        }
	if (j == 3)
        {
          hyield[j]->SetLineColor(kBlue+2);
          hyield[j]->SetLineWidth(3);
          hyield[j]->SetLineStyle(2);
        }
		
        if (j == 6)
        {
          hyield[j]->SetFillColor(kRed);
          hyield[j]->SetLineColor(kRed);
          hyield[j]->SetLineWidth(3);
          hyield[j]->SetFillStyle(0);

	}

	if (j == 4)
	  {
	    hyield[j]->SetLineColor(kGreen+2);
	    hyield[j]->SetLineWidth(4);
	    hyield[j]->SetLineStyle(1);
	  }
	      
	if (j == 5)
	  {
	    hyield[j]->SetLineColor(kGreen-2);
	    hyield[j]->SetLineWidth(3);
	    hyield[j]->SetLineStyle(2);
	  }
	
                
        



        hyield[j]->Draw("histoesame");
	if(step==0) {
	  hyield[j]->GetYaxis()->SetRangeUser(0, 5);
	  hyield[j]->GetYaxis()->SetTitle("Yield,   0<|#Delta y|<1.6");
	} else {
	  hyield[j]->GetYaxis()->SetRangeUser(0, 5);
          hyield[j]->GetYaxis()->SetTitle("Yield,   1.6<|#Delta y|<5.0");
	}
        if  (j == 6)
          leg->AddEntry(hyield[j], "#font[42]{" + title_samples[j] +"}", "f");
        else
          leg->AddEntry(hyield[j], "#font[42]{" + title_samples[j] + "}","l");// & "+title_samples[6]+ "}", "l");
      }
      leg->SetFillStyle(0);
      leg->SetLineWidth(0);
      leg->SetLineColor(0);
      leg->SetBorderSize(0);

      leg->Draw();

      LabelsReco();
      if(bkg==true) canvas1->Print(TString::Format("plots/Yield_reco_range%i_cut%i.eps",step,cuts));
      else canvas1->Print(TString::Format("plots/reco_onlyHV_range%i_cut%i.eps",step,cuts));

      /*TCanvas *canvas2 = new TCanvas(TString::Format("ratio_yieldcanvas_range%i",step),TString::Format("ratio_yieldcanvas_range%i",step), 800, 400);
      canvas2->cd();
      gPad->SetGridy();
      TLegend *leg2 = new TLegend(0.18, 0.55, 0.65, 0.9); //(0.4,0.3,0.5,0.6);
      leg2->SetTextSize(0.055);

       ///       for (int j = 0; j < h2_signal.size(); j++) {
      //	hyield[j]->Rebin(5);
      //	}

      
       for (int j = 0; j < h2_signal.size()-1; j++)
      {

	Int_t n = 0;
	Double_t x[100];
	Double_t y[100];

	for(int i=1; i<hyield[j]->GetNbinsX()+1; i++ ) {
	  x[i-1]=hyield[j]->GetBinCenter(i);
	  y[i-1]= (hyield[j]->GetBinContent(i) / hyield[4]->GetBinContent(i));///sqrt(pow(hyield[j]->GetBinError(i),2) + pow(hyield[h2_signal.size()-1]->GetBinError(i),2));
	  n++;
	}
	
	g_ratio_yield[j]=new TGraph(n,x,y);

	 if (j == 0)
        {
          g_ratio_yield[j]->SetLineColor(kBlue-2);
          g_ratio_yield[j]->SetLineWidth(6);
        }
        if (j == 1)
        {
          g_ratio_yield[j]->SetLineColor(kBlue);
          g_ratio_yield[j]->SetLineWidth(4);
        }
	if (j == 2)
        {
          g_ratio_yield[j]->SetLineColor(kBlue+2);
          g_ratio_yield[j]->SetLineWidth(2);
          g_ratio_yield[j]->SetLineStyle(1);
        }
	if (j == 3)
        {
          g_ratio_yield[j]->SetLineColor(kBlue+2);
          g_ratio_yield[j]->SetLineWidth(3);
          g_ratio_yield[j]->SetLineStyle(2);
        }
		
  	
        if(j==0) {
	  g_ratio_yield[j]->Draw("al");
	  g_ratio_yield[j]->GetYaxis()->SetTitle("Y_{HV+SM}/Y_{SM}");
	  g_ratio_yield[j]->GetXaxis()->SetTitle("#Delta#phi");
	  g_ratio_yield[j]->GetYaxis()->SetTitleOffset(0.8);
	  g_ratio_yield[j]->GetYaxis()->SetTitleSize(0.08);
	  g_ratio_yield[j]->GetYaxis()->SetRangeUser(0.5, 4);
	  g_ratio_yield[j]->GetYaxis()->SetLabelSize(0.06);
	  g_ratio_yield[j]->GetXaxis()->SetLabelSize(0.06);
	  g_ratio_yield[j]->GetXaxis()->SetTitleOffset(1.);
	  g_ratio_yield[j]->GetXaxis()->SetTitleSize(0.08);
	  
	} else g_ratio_yield[j]->Draw("l");

      }

       Int_t n = 0;
       Double_t x[100];
       Double_t y[100];
       Double_t ey[100];

       for(int i=1; i<hyield[4]->GetNbinsX()+1; i++ ) {
	 x[i-1]=hyield[4]->GetBinCenter(i);
	 y[i-1]= 1;
	 ey[i-1]=sqrt(hyield[4]->GetBinError(i)*hyield[4]->GetBinError(i)+hyield[0]->GetBinError(i)*hyield[0]->GetBinError(i));
	 n++;
       }
	
    TGraphErrors *g_ratio_err=new TGraphErrors(n,x,y,0,ey);

    g_ratio_err->SetFillColor(kRed);
    g_ratio_err->SetLineColor(kRed);
    g_ratio_err->SetLineWidth(0);
    g_ratio_err->SetFillStyle(3001);
    g_ratio_err->Draw("l3");

    for(int j=0; j<4; j++ ){
      leg2->AddEntry(g_ratio_yield[j], "#font[42]{" + title_samples[j] + " & "+title_samples[4]+ "}", "l");
    }
    leg2->AddEntry(g_ratio_err, "#font[42]{" + title_samples[4] + " (+ stat. unc. 2 ab^{-1})}", "f");
    
    leg2->SetFillStyle(1001);
    leg2->SetFillColor(0);
    leg2->SetLineWidth(0);
    leg2->SetLineColor(0);
    leg2->SetBorderSize(0);

    leg2->Draw();
    

    LabelsReco(0.066);
    if(step==0) QQBARLabel2(0.28, 0.43, TString::Format("#font[42]{Yield, 0<|#Delta y|<1.6}"), kGreen+2 ,0.065);
    if(step==1) QQBARLabel2(0.28, 0.43, TString::Format("#font[42]{Yield, 1.6<|#Delta y|<5}"), kGreen+2 ,0.065);

      if(bkg==true) canvas2->Print(TString::Format("ratio_Yield_reco_range%i_cut%i.eps",step,cuts));
      else canvas2->Print(TString::Format("ratio_Yield_reco_onlyHV_range%i_cut%i.eps",step,cuts));
      
      
*/
	}
      
  }
}



void C2_yieldplots(int cuts=6)
{

  TString folder_0 = "../results_cut6_npfo22_nchpfo15/";
    
  plots_reco_C2Yield_PRL( 2000, true, true,folder_0,cuts);
  //plots_reco_C2Yield_PRL( 2000, true, false,folder_0,cuts);

  //if second argument is false, only one proces is used
  //if it is itrue, S2= S2HV+S2SM using cuts 6
  // for both cases, B2=B2SM cuts 20

  //if the third argument is true --> we show Yield.
  //if false --> C2
  
}


