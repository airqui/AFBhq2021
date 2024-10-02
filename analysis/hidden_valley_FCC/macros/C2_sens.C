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
int nSM=4;
//if we study ISR & VBF --> nbkgs=2

std::vector<TH1F*> v_hyield (  std::vector<std::vector<TH1F *>> h1_signal, std::vector<std::vector<TH2F *>> h2_signal,  std::vector<std::vector<TH2F *>> h2_bkg, int step=0) {

  std::vector<TH1F*> result;
    
  for (int j = 0; j < h2_signal.size(); j++)
    {

      TH1F *hyield_temp;
      
      TF1 * minus;
      TF1 * plus;
    
      minus = new TF1(TString::Format("minus_j%i",j),"pol8",-TMath::Pi() / 2.,0);
      plus = new TF1(TString::Format("plus_j%i",j),"pol8",0, 3.*TMath::Pi() / 2);
	
	
      int first=50;
      int bins=16;
	
      TH1F* htemp1;
      TH1F* htemp1_SM1;
      TH1F* htemp1_SM2;

      TH1F* htemp2;
      TH1F* htemp3;
	
      if(step==1) {
	hyield_temp = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
	htemp1 = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp1_%i",j), 1,first-step*bins);//, first-2*step*bins*2);
	hyield_temp->Add(htemp1);

      } else {
	hyield_temp = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first-bins,first+bins);//, first+2*step*bins*2);
      }

      if(j != (h2_signal.size()-1)) {
	float na=h1_signal.at(j).at(1)->Integral();
	float nb=h1_signal.at(h2_signal.size()-1).at(1)->Integral();
	hyield_temp->Scale(na/(na+nb));

	if(step==1) {
	  htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first+step*bins,first*2);//, first+2*step*bins*2);
	  htemp1_SM2 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp1_%i_SM2",j),1, first-step*bins);//, first-2*step*bins*2);
	  htemp1_SM1->Add(htemp1_SM2);
	} else {
	  htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first-bins,first+bins);//, first+2*step*bins*2);
	}
	htemp1_SM1->Scale(nb/(na+nb));
	hyield_temp->Add(htemp1_SM1);
      }
	
      htemp2 = (TH1F *)h2_bkg.at(0).at(1)->ProjectionY(TString::Format("htemp2_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
      htemp3 = (TH1F *)h2_bkg.at(0).at(1)->ProjectionY(TString::Format("htemp3_%i",j), 1,first-step*bins);//,first-2*step*bins*2);
      htemp2->Add(htemp3);

	
      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(htemp2->GetBinContent(i)!=0 && hyield_temp->GetBinContent(i)!=0 ) {
	  double a=hyield_temp->GetBinContent(i);
	  double b=htemp2->GetBinContent(i);
	  double ea=hyield_temp->GetBinError(i);
	  double eb=htemp2->GetBinError(i);
	  hyield_temp->SetBinError(i,(1./b)*sqrt(pow(ea,2)+pow(a*eb/b,2)));
	  hyield_temp->SetBinContent(i,a/b);//hyield_temp->GetBinContent(i)/htemp2->GetBinContent(i));
	    
	} else {
	  hyield_temp->SetBinContent(i,0);
	  hyield_temp->SetBinError(i,0);
	}
      }

      ///	if(j==h2_signal.size()-1) {
	
      //TF1 * minus = new TF1("minus","pol8",-TMath::Pi() / 2.,0);
      hyield_temp->Fit(TString::Format("minus_j%i",j),"REMN");
	  
      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(hyield_temp->GetBinCenter(i)<=0 ) {

	  float params[20]={0};
	  for(int ipol=0; ipol<9; ipol++) params[ipol]=minus->GetParameter(ipol);
	  float value = params[0];
	  for(int ipol=1; ipol<9; ipol++) {
	    value+=params[ipol]*pow(hyield_temp->GetBinCenter(i),ipol);
	  }
	  //	  if(j==nSM) hyield_temp->SetBinContent(i,value);
	  // uncomment if you want to fit the SM contribution... 
	}
      }

      hyield_temp->Fit(TString::Format("plus_j%i",j),"REMN");

      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(hyield_temp->GetBinCenter(i)>0 ) {

	  float params[20]={0};
	  for(int ipol=0; ipol<9; ipol++) params[ipol]=plus->GetParameter(ipol);
	  float value = params[0];
	  for(int ipol=1; ipol<9; ipol++) {
	    value+=params[ipol]*pow(hyield_temp->GetBinCenter(i),ipol);
	  }
	  //if(j==nSM) hyield_temp->SetBinContent(i,value);
           // uncomment if you want to fit the SM contribution...    

	}
      }


      hyield_temp->SetName(TString::Format("%s_yield%i",samples[j].Data(),step));
 
      result.push_back(hyield_temp);

    }

  return result;
}


std::vector<TH1F*> v_hyield_2 (  std::vector<std::vector<TH1F *>> h1_signal, std::vector<std::vector<TH2F *>> h2_signal,  int step=0) {

  std::vector<TH1F*> result;
  
   
  for (int j = 0; j < h2_signal.size(); j++)
    {

      TH1F *hyield_temp;
      
      TF1 * minus;
      TF1 * plus;
    
      minus = new TF1(TString::Format("minus_j%i",j),"pol8",-TMath::Pi() / 2.,0);
      plus = new TF1(TString::Format("plus_j%i",j),"pol8",0, 3.*TMath::Pi() / 2);
	
	
      int first=50;
      int bins=16;
	
      TH1F* htemp1;
      TH1F* htemp1_SM1;
      TH1F* htemp1_SM2;

      TH1F* htemp2;
      TH1F* htemp3;
	
      if(step==1) {
	hyield_temp = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
	htemp1 = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp1_%i",j), 1,first-step*bins);//, first-2*step*bins*2);
	hyield_temp->Add(htemp1);

      } else {
	hyield_temp = (TH1F *)h2_signal.at(j).at(0)->ProjectionY(TString::Format("htemp_%i",j), first-bins,first+bins);//, first+2*step*bins*2);
      }

      if(j != (h2_signal.size()-1)) {
	float na=h1_signal.at(j).at(1)->Integral();
	float nb=h1_signal.at(h2_signal.size()-1).at(1)->Integral();
	hyield_temp->Scale(na/(na+nb));

	if(step==1) {
	  htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first+step*bins,first*2);//, first+2*step*bins*2);
	  htemp1_SM2 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp1_%i_SM2",j),1, first-step*bins);//, first-2*step*bins*2);
	  htemp1_SM1->Add(htemp1_SM2);
	} else {
	  htemp1_SM1 = (TH1F *)h2_signal.at(h2_signal.size()-1).at(0)->ProjectionY(TString::Format("htemp_%i_SM1",j), first-bins,first+bins);//, first+2*step*bins*2);
	}
	htemp1_SM1->Scale(nb/(na+nb));
	hyield_temp->Add(htemp1_SM1);
      }
	
      htemp2 = (TH1F *)h2_signal.at(0).at(1)->ProjectionY(TString::Format("htemp2_%i",j), first+step*bins,first*2);//, first+2*step*bins*2);
      htemp3 = (TH1F *)h2_signal.at(0).at(1)->ProjectionY(TString::Format("htemp3_%i",j), 1,first-step*bins);//,first-2*step*bins*2);
      htemp2->Add(htemp3);

	
      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(htemp2->GetBinContent(i)!=0 && hyield_temp->GetBinContent(i)!=0 ) {
	  double a=hyield_temp->GetBinContent(i);
	  double b=htemp2->GetBinContent(i);
	  double ea=hyield_temp->GetBinError(i);
	  double eb=htemp2->GetBinError(i);
	  hyield_temp->SetBinError(i,(1./b)*sqrt(pow(ea,2)+pow(a*eb/b,2)));
	  hyield_temp->SetBinContent(i,a/b);//hyield_temp->GetBinContent(i)/htemp2->GetBinContent(i));
	    
	} else {
	  hyield_temp->SetBinContent(i,0);
	  hyield_temp->SetBinError(i,0);
	}
      }

      ///	if(j==h2_signal.size()-1) {
	
      //TF1 * minus = new TF1("minus","pol8",-TMath::Pi() / 2.,0);
      hyield_temp->Fit(TString::Format("minus_j%i",j),"REMN");
	  
      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(hyield_temp->GetBinCenter(i)<=0 ) {

	  float params[20]={0};
	  for(int ipol=0; ipol<9; ipol++) params[ipol]=minus->GetParameter(ipol);
	  float value = params[0];
	  for(int ipol=1; ipol<9; ipol++) {
	    value+=params[ipol]*pow(hyield_temp->GetBinCenter(i),ipol);
	  }
	  if(j==nSM) hyield_temp->SetBinContent(i,value);
	}
      }

      hyield_temp->Fit(TString::Format("plus_j%i",j),"REMN");

      for(int i=0; i<hyield_temp->GetNbinsX()+1; i++ ) {
	if(hyield_temp->GetBinCenter(i)>0 ) {

	  float params[20]={0};
	  for(int ipol=0; ipol<9; ipol++) params[ipol]=plus->GetParameter(ipol);
	  float value = params[0];
	  for(int ipol=1; ipol<9; ipol++) {
	    value+=params[ipol]*pow(hyield_temp->GetBinCenter(i),ipol);
	  }
	  if(j==nSM) hyield_temp->SetBinContent(i,value);
	}
      }


      hyield_temp->SetName(TString::Format("yield%i",step));
 
      result.push_back(hyield_temp);

    }

  return result;
}

void histos_ref_C2Yield_PRL(float lum=2000, TString folder_0="../results_new/", int cuts=6, TString MC="_MC")
{

  std::vector<std::vector<TH1F *>> h1_signal;
  std::vector<std::vector<TH2F *>> h2_signal;

  std::vector<std::vector<TH1F *>> h1_bkg;
  std::vector<std::vector<TH2F *>> h2_bkg;

  std::vector<TString> histos_1d = {"ncharged"};
  std::vector<TString> histos_2d = {"S2"+MC+"_rapidity", "B2"+MC+"_rapidity"};

  //samples has to be the HV100,50,10,ISR
  for (int isample = 0; isample < sizeof(samples) / sizeof(TString)-1 ; isample++)
    {
      cout<<isample<<endl;
      //folder=TString::Format("%s/selection_cuts%i",folder_0.Data(),cuts);
      folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);

      //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
      std::vector<TH1F *> h1_temp = GetHisto1D(samples[isample], lum, histos_1d);
      h1_signal.push_back(h1_temp);
      //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
      std::vector<TH2F *> h2_temp= GetHisto2D(samples[isample], lum, histos_2d,true);

      h2_signal.push_back(h2_temp);
      //h2_signal_MC.push_back(h2_temp_MC);

    }

  folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),20);
  
  for (int isample = sizeof(samples) / sizeof(TString); isample < sizeof(samples) / sizeof(TString)+1; isample++)
    {
      std::vector<TH1F *> h1_temp = GetHisto1D(samples[isample-1], lum, histos_1d);//we read again the ISR sample, but this time with cuts 20                                                                 
      h1_bkg.push_back(h1_temp);
      std::vector<TH2F *> h2_temp = GetHisto2D(samples[isample-1], lum, histos_2d,true);
      
      h2_bkg.push_back(h2_temp);
    }
  

  TFile *MyFile = new TFile(TString::Format("rootfiles/yields%s_lum%i.root",MC.Data(),int(lum)), "RECREATE");
  for(int step=0; step<2; step++) {

    
    std::vector<TH1F*> hyield=v_hyield(h1_signal,h2_signal,h2_bkg,step);
      int linest=2;
    if(MC!="_MC") linest=1;
    for(int i=0; i<hyield.size(); i++)  hyield.at(i)->SetLineStyle(linest);

    MyFile->cd();
    for(int i=0; i<hyield.size(); i++)  hyield.at(i)->Write();
  }
      
}



void histos_ref_C2Yield_PRL_fragmentation(TString folder_0="../results_new/")
{

  //Function to draw the yield using S2 and B2 from the same type of events, QQBAR, with no cuts.
  //it is used to compare pythia and Herwig
  //Analysis at MC-particle level (not detector level)
  
  std::vector<std::vector<TH1F *>> h1_signal;
  std::vector<std::vector<TH2F *>> h2_signal;

  std::vector<std::vector<TH1F *>> h1_bkg;
  std::vector<std::vector<TH2F *>> h2_bkg;

  std::vector<TString> histos_1d = {"ncharged"};
  std::vector<TString> histos_2d = {"S2_MC_rapidity", "B2_MC_rapidity"};

  
  TString samples2[] = {"eetoDvDvbar_m62.5"};
  float lum=2000;
  //samples has to be the HV100,50,10,ISR
  for (int isample = 0; isample < sizeof(samples2) / sizeof(TString) ; isample++)
    {
      cout<<isample<<endl;
      //folder=TString::Format("%s/selection_cuts%i",folder_0.Data(),cuts);
      folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),nSM);

      //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
      cout<<folder<<endl;
      std::vector<TH1F *> h1_temp = GetHisto1D(samples2[isample], lum, histos_1d);
      h1_signal.push_back(h1_temp);
      //folder=TString::Format("%s/QCDcorrelations_cuts%i",folder_0.Data(),cuts);
      std::vector<TH2F *> h2_temp= GetHisto2D(samples2[isample], lum, histos_2d,true);

      h2_signal.push_back(h2_temp);
      //h2_signal_MC.push_back(h2_temp_MC);


      TFile *MyFile = new TFile(TString::Format("rootfiles/yields_%s.root",samples2[isample].Data()), "RECREATE");
      for(int step=0; step<2; step++) {
    
	std::vector<TH1F*> hyield=v_hyield_2(h1_signal,h2_signal,step);
	int linest=isample+1;
    
	for(int i=0; i<hyield.size(); i++)  hyield.at(i)->SetLineStyle(linest);

	MyFile->cd();
	for(int i=0; i<hyield.size(); i++)  hyield.at(i)->Write();
      }
    }
  
}


void plots_C2Yield_syst_fragmentation(float lum=2000, int yield=0, double factor_det=1, double factor_frag=1) {

  //sensitivity for
 
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

  gStyle->SetPadLeftMargin(0.18);
  
  TFile *MyFile =  new TFile(TString::Format("rootfiles/yields_lum%i.root",int(lum)));
  TFile *MyFileMC =  new TFile(TString::Format("rootfiles/yields_MC_lum%i.root",int(lum)));

  TFile *MyFile_pythia =  new TFile("rootfiles/yields_qqbar_Pythia.root");
  TFile *MyFile_herwig =  new TFile("rootfiles/yields_qqbar_Herwig.root");

  //*******************************************************
  // MCparticle level vs Reco level uncertainty
  // here we include the statistical uncertainty calculation of given luminosity
  TH1F *hyield[10], *hyieldMC[10];

  for(int j=0; j<nSM+1; j++) {
    //(TH1F *)f
    hyield[j] = (TH1F *)MyFile->Get(TString::Format("%s_yield%i",samples[j].Data(),yield));
    hyieldMC[j] = (TH1F*)MyFileMC->Get(TString::Format("%s_yield%i",samples[j].Data(),yield));
    cout<<TString::Format("%s_yield%i",samples[j].Data(),yield)<<endl;
    hyield[j]->Rebin(4);
    hyieldMC[j]->Rebin(4);
   
    if (j == 0)
      {
	hyield[j]->SetLineColor(kBlue-2);
	hyield[j]->SetLineWidth(4);
      }
    if (j == 1)
      {
	hyield[j]->SetLineColor(kBlue);
	hyield[j]->SetLineWidth(3);
      }
    if (j == 2)
      {
	hyield[j]->SetLineColor(kBlue+2);
	hyield[j]->SetLineWidth(2);
	//hyield[j]->SetLineStyle(1);
      }
    if (j == 3)
      {
	hyield[j]->SetLineColor(kBlue+2);
	hyield[j]->SetLineWidth(1);
	//hyield[j]->SetLineStyle(1);
      }

    /*      if (j == 4)
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
    */  

    
    if (j == nSM)
      {
	hyield[j]->SetFillColor(kRed);
	hyield[j]->SetLineColor(kRed);
	hyield[j]->SetLineWidth(3);
	hyield[j]->SetFillStyle(0);
      }




    if(yield==0) {
      hyield[j]->GetYaxis()->SetRangeUser(0, 5);
      hyield[j]->GetYaxis()->SetTitle("Yield,   0<|#Delta y|<1.6");
    } else {
      hyield[j]->GetYaxis()->SetRangeUser(0, 5);
      hyield[j]->GetYaxis()->SetTitle("Yield,   1.6<|#Delta y|<5.0");
    }

    for(int i1=0; i1<hyield[j]->GetNbinsX(); i1++) {
      float y1=hyield[j]->GetBinContent(i1);
      float y2=hyieldMC[j]->GetBinContent(i1);
      float err= y1-y2;
      if(fabs(err)<0.15 && i1>0 && i1<hyield[j]->GetNbinsX()-1) {
	cout<<"    eeeeeeeeeee "<<err<<" "<<endl;
	err=0.5* ( fabs(hyield[j]->GetBinContent(i1-1)-hyieldMC[j]->GetBinContent(i1-1)) + fabs(hyield[j]->GetBinContent(i1+1)-hyieldMC[j]->GetBinContent(i1+1)) );
	err=err*fabs(hyield[j]->GetBinContent(i1+1)-hyieldMC[j]->GetBinContent(i1+1))/(hyield[j]->GetBinContent(i1+1)-hyieldMC[j]->GetBinContent(i1+1));
      }
      if(hyield[j]->GetBinCenter(i1)>1.7 && hyield[j]->GetBinCenter(i1)<1.9 && yield==0) err*=2;
      if(hyield[j]->GetBinCenter(i1)>1.9 && hyield[j]->GetBinCenter(i1)<2.1 && yield==0) err*=1.3;

      cout<<hyield[j]->GetBinCenter(i1)<<" "<<i1<<" "<<y1<<" "<<y2<<" "<<err<<" "<<(y1-y2)/err<<endl;
      //      hyield[j]->SetBinError(i1,0);
      // err/=2.;
      if(factor_det<1)  hyield[j]->SetBinError(i1,sqrt(pow(hyield[j]->GetBinError(i1),2)+factor_det*2*pow(err,2)));
      else hyield[j]->SetBinError(i1,sqrt(pow(hyield[j]->GetBinError(i1),2)+pow(err,2)));

    }
  }


  //*******************************************************
  // Pythia vs Herwig uncertainty
  // ---------------------------
 TH1F *hyield_p, *hyield_h;

 //(TH1F *)f
 hyield_p = (TH1F *)MyFile_pythia->Get(TString::Format("yield%i",yield));
 hyield_h = (TH1F*)MyFile_herwig->Get(TString::Format("yield%i",yield));

 hyield_p->Rebin(4);
 hyield_h->Rebin(4);
   
 for(int i1=0; i1<hyield_p->GetNbinsX(); i1++) {
   float y1=hyield_p->GetBinContent(i1);
   float y2=hyield_h->GetBinContent(i1);
   float err= y1-y2;
   if(fabs(err)<0.05 && i1>0 && i1<hyield_p->GetNbinsX()-1) {
     //     cout<<"    eeeeeeeeeee "<<err<<" "<<endl;
     err=0.5* ( fabs(hyield_p->GetBinContent(i1-1)-hyield_h->GetBinContent(i1-1)) + fabs(hyield_p->GetBinContent(i1+1)-hyield_p->GetBinContent(i1+1)) );
     err=err*fabs(hyield_p->GetBinContent(i1+1)-hyield_h->GetBinContent(i1+1))/(hyield_p->GetBinContent(i1+1)-hyield_h->GetBinContent(i1+1));
   }
   //cout<<"pythia-herwig"<<i1<<" "<<hyield_p->GetBinCenter(i1)<<" "<<y1<<" "<<y2<<" "<<err<<" "<<(y1-y2)/err<<endl;
   hyield_p->SetBinError(i1,sqrt(pow(factor_frag*err/y1,2)));
 }
  

 //************************************************
 //PLOTTING
  TCanvas *canvas1 = new TCanvas(TString::Format("yieldcanvas_range%i",yield),TString::Format("yieldcanvas_range%i",yield), 800, 800);
  canvas1->cd();
    
  TLegend *leg = new TLegend(0.2, 0.72, 0.45, 0.92); //(0.4,0.3,0.5,0.6);
  leg->SetTextSize(0.035);
  
  for (int j = 0; j < nSM+1; j++) {
    hyield[j]->Draw("histoesame");
    //hyield_[j]->Draw("histoesame");
    
    if  (j == nSM) {
      leg->AddEntry(hyield[j], "#font[42]{" + title_samples[j] +"}", "f");
      //leg->AddEntry(hyield_[j], "#font[42]{" + title_samples[j] +"}", "f");
    }
    else {
      leg->AddEntry(hyield[j], "#font[42]{" + title_samples[j] + "}","l");// & "+title_samples[6]+ "}", "l");
      //leg->AddEntry(hyield_[j], "#font[42]{" + title_samples[j] + " & "+title_samples[4]+ "}", "l");
    }
    
  }
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  leg->SetBorderSize(0);
  
  leg->Draw();
  
  LabelsReco();
  //canvas1->Print(TString::Format("plots/Yield_reco_range%i_syst.eps",yield));
  // canvas1->Print(TString::Format("plots/Yield_reco_range%i_syst.pdf",yield));

  
  TCanvas *canvas2 = new TCanvas(TString::Format("ratio_yieldcanvas_range%i",yield),TString::Format("ratio_yieldcanvas_range%i",yield), 800, 800);
  canvas2->cd();
  //  gPad->SetGridy();
  TLegend *leg2 = new TLegend(0.25, 0.6, 0.55, 0.8); //(0.4,0.3,0.5,0.6);
  leg2->SetTextSize(0.03);

  TGraph * g_ratio_yield[10];

  double xg[2]={-2,5};
  double y3[2]={3,3};
  double y5[2]={5,5};

  TGraph * g_3sigma = new TGraph(2,xg,y3);
  TGraph * g_5sigma = new TGraph(2,xg,y5);

  for (int j = 0; j < nSM+1; j++)
    { 

      Int_t n = 0;
      Double_t x[100];
      Double_t y[100];

      for(int i=1; i<hyield[j]->GetNbinsX()+1; i++ ) {
	x[i-1]=hyield[j]->GetBinCenter(i);
	//uncertainty is the quadrature summ of SM stat uncertainty, MC-detector unc. and pythia-herwig unc.
	//the stat and MC-det uncertainty are both coded in hyield[4]->GetBinError(i)
	//the hyield_p histograms contains the relative unc pythia-herwig. We are interested in the total
	y[i-1]= fabs(hyield[j]->GetBinContent(i) -hyield[nSM]->GetBinContent(i))/ sqrt(pow(hyield[j]->GetBinError(i),2)+pow(hyield[nSM]->GetBinError(i),2)+pow(0.1*factor_frag*hyield_p->GetBinError(i)*hyield[nSM]->GetBinError(i),2)  ) ;
	
	//	if(factor_det==1 && factor_frag==0.1) y[i-1]= fabs(hyield[j]->GetBinContent(i) -hyield[4]->GetBinContent(i))/ sqrt(pow(hyield[j]->GetBinError(i),2)+pow(factor_det*0.5*hyield[4]->GetBinError(i),2)+pow(factor_frag*hyield_p->GetBinError(i)*hyield[j]->GetBinError(i),2)  ) ;

	n++;
      }
	
      g_ratio_yield[j]=new TGraph(n-2,x,y);

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

      /*
	      if (j == 4)
          {
            g_ratio_yield[j]->SetLineColor(kGreen+2);
            g_ratio_yield[j]->SetLineWidth(4);
            g_ratio_yield[j]->SetLineStyle(1);
          }
              
        if (j == 5)
          {
            g_ratio_yield[j]->SetLineColor(kGreen-2);
            g_ratio_yield[j]->SetLineWidth(3);
            g_ratio_yield[j]->SetLineStyle(2);
          }
      */
  	
      if(j==0) {
	g_ratio_yield[j]->Draw("ac");
	g_ratio_yield[j]->GetYaxis()->SetTitle("|Y_{HV+SM} - Y_{SM}| / #sigma_{Y}");
	g_ratio_yield[j]->GetXaxis()->SetTitle("#Delta#phi");
	g_ratio_yield[j]->GetYaxis()->SetTitleOffset(1.4);
	g_ratio_yield[j]->GetYaxis()->SetTitleSize(0.06);
	g_ratio_yield[j]->GetYaxis()->SetRangeUser(0, 22);
	g_ratio_yield[j]->GetYaxis()->SetLabelSize(0.04);
	g_ratio_yield[j]->GetXaxis()->SetLabelSize(0.04);
	g_ratio_yield[j]->GetXaxis()->SetTitleOffset(1.2);
	g_ratio_yield[j]->GetXaxis()->SetTitleSize(0.06);
	  
      } else g_ratio_yield[j]->Draw("c");

    }

  Int_t n = 0;
  Double_t x[100];
  Double_t y[100];
  Double_t ey[100];

  for(int i=1; i<hyield[nSM]->GetNbinsX()+1; i++ ) {
    x[i-1]=hyield[nSM]->GetBinCenter(i);
    y[i-1]= 1;
    ey[i-1]=sqrt(hyield[nSM]->GetBinError(i)*hyield[nSM]->GetBinError(i)+hyield[0]->GetBinError(i)*hyield[0]->GetBinError(i));
    n++;
  }
	
 
  for(int j=0; j<nSM; j++ ){
    leg2->AddEntry(g_ratio_yield[j], "#font[42]{" + title_samples[j] + "}","l");// & "+title_samples[6]+ "}", "l");
  }
  //leg2->AddEntry(g_ratio_err, "#font[42]{" + title_samples[4] + " (+ stat. unc. 2 ab^{-1})}", "f");
    
  leg2->SetFillStyle(1001);
  leg2->SetFillColor(0);
  leg2->SetLineWidth(0);
  leg2->SetLineColor(0);
  leg2->SetBorderSize(0);

  leg2->Draw();


  g_3sigma->SetLineColor(kOrange+1);
  g_3sigma->SetLineStyle(2);
  g_3sigma->SetLineWidth(2);
  g_3sigma->Draw("l");
    
  g_5sigma->SetLineColor(kRed+1);
  g_5sigma->SetLineStyle(2);
  g_5sigma->SetLineWidth(2);
  g_5sigma->Draw("l");
  QQBARLabel2(0.88, 0.35, "5#sigma", kRed+1 ,0.045);
  QQBARLabel2(0.88, 0.278, "3#sigma", kOrange+1 ,0.045);
  
    
  LabelsReco(0.04);

  //QQBARLabel2(0.23, 0.88, TString::Format("Stats#oplus Optimistic Syst (%i fb^{-1}#oplus #frac{det. modeling}{10} #oplus #frac{fragmentation}{10})",int(lum)), kBlack ,0.025);
  if(factor_det==1 && factor_frag==0.1) QQBARLabel2(0.23, 0.88, TString::Format("Stats#oplus Optimistic Syst (%i fb^{-1}#oplus det. modeling #oplus #frac{fragmentation}{10})",int(lum)), kBlack ,0.025);
  if(factor_det==0.1 && factor_frag==0.1) QQBARLabel2(0.23, 0.88, TString::Format("Stats#oplus Optimistic Syst (%i fb^{-1}#oplus #frac{det. modeling}{10} #oplus #frac{fragmentation}{10})",int(lum)), kBlack ,0.025);
  if(factor_det==0.1 && factor_frag==1) QQBARLabel2(0.23, 0.88, TString::Format("Stats#oplus Optimistic Syst (%i fb^{-1}#oplus #frac{det. modeling}{10} #oplus fragmentation)",int(lum)), kBlack ,0.025);
  if(factor_det==1 && factor_frag==1) QQBARLabel2(0.23, 0.88, TString::Format("Stats#oplus Optimistic Syst (%i fb^{-1}#oplus det. modeling #oplus fragmentation)",int(lum)), kBlack ,0.025);


  //QQBARLabel2(0.23, 0.88,"Syst (Frag)", kBlack ,0.025);                                                            

  if(yield==0) QQBARLabel2(0.23, 0.83, TString::Format("#font[42]{Yield, 0<|#Delta y|<1.6}"), kRed+1 ,0.035);
  if(yield==1) QQBARLabel2(0.23, 0.83, TString::Format("#font[42]{Yield, 1.6<|#Delta y|<5}"), kRed+1 ,0.035);

  TString extra="";
  if(factor_det==1 && factor_frag==0.1) extra="_FragDiv10";
  if(factor_det==0.1 && factor_frag==0.1) extra="_Detdiv10_Fragdiv10";
  if(factor_det==0.1 && factor_frag==1) extra="_Detdiv10";
  //if(factor_det==1 && factor_frag==1)


  canvas2->Print(TString::Format("plots/Sens_Yield_range%i_lum%i%s.eps",yield,int(lum),extra.Data()));
  canvas2->Print(TString::Format("plots/Sens_Yield_range%i_lum%i%s.pdf",yield,int(lum),extra.Data()));
  //canvas2->Print(TString::Format("diff_Yield_range%i.eps",yield));
  //canvas2->Print(TString::Format("diff_Yield_range%i.pdf",yield));

  
      
  
}




void C2_sens(int cuts=6)
{

  //fragmentation analysis
  TString folder_0 = "../results_FCC/";

  //study of Detector Level Systematics
  //comparing the yields calculated at MC-level vs the detector level with S and B as defined in the draft
  // the sensitivity is defined as |YSM-YBSM| / Error_YSM with Error_YSM=|YSM_MC - YSM_detector|
  //  TString folder_0 = "../results_cut6_npfo22_nchpfo15/";
  //int cuts=6;

  //THIS PART IS ONLY NEEDED ONCE --> it stores the histograms in a root file in the folder root-file
  for(float lum=100; lum<200; lum+=100) {
    histos_ref_C2Yield_PRL( lum, folder_0, cuts);
    histos_ref_C2Yield_PRL( lum, folder_0,cuts,"");

  }
     
  // THIS PART is the part that calculates the actual sensitivity
   for(float lum=100; lum<200; lum+=100) {
    for(int i=0; i<2; i++) {
      plots_C2Yield_syst_fragmentation(lum,i,1,1);
      plots_C2Yield_syst_fragmentation(lum,i,0.1,0.1);
      //plots_C2Yield_syst_fragmentation(lum,i,1,0.1);
      //plots_C2Yield_syst_fragmentation(lum,i,0.1,1);
    }
   }
  

  
}


