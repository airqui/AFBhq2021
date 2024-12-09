TString folder = "../results_FCC/";


TString samples[] = {
  "HV_240_mDv120_mqv100",
  "HV_240_mDv120_mqv10",
  "HV_240_mDv100_mqv50",
  "ISR_240",
  "qqSM_240"};


TString title_samples[] = {
  "#font[12]{m_{D_{v}}= 120 GeV, m_{q_{v}}=100 GeV}",
  "#font[12]{m_{D_{v}}= 120 GeV, m_{q_{v}}= 10 GeV}",
  "#font[12]{m_{D_{v}}= 100 GeV, m_{q_{v}}= 50 GeV}",
  "SM"};


// NOT USED --->
TString samples2[] = {
  "ISR_240",
  "ISR"};                                                                                                                                                                                                                                                                                                                                                                                  
TString title_samples2[] = { 
  "qqbar+ISR-240GeV, FCCdet",
  "qqbar+ISR-250GeV, ILCdet"};
//< --- NOT USED

void Labels()
{

  //QQBARLabel(ildx, 0.953, "");
  QQBARLabel2(0.04, 0.04,"Analysis at detector level, Pythia8+SGV", kOrange + 3);

}

void LabelsReco(Double_t textsize=0.03)
{

  //QQBARLabel(ildx, 0.953, "");
  QQBARLabel2(0.04, 0.04,"Pythia8+SGV (FCC240, ILD-for-FCC detector)", kOrange + 3,textsize);

}

std::vector<TH1F *> GetHisto1D(TString sample = "2f_hadronic_sample", float lum = 900, std::vector<TString> histonames={"h_mjj","h_mj1_mj2","h_y23","h_d23","h_thrust","h_major_thrust","h_minor_thrust"}, TString pol="unpol")
{


  std::vector<TH1F *> h1;
  std::vector<TH1F *> h0;

  TString filename = folder + "_" + sample +  ".root";
 
  TFile *f = new TFile(filename);
  TH1F *h_luminosity_cross_2f = (TH1F *)f->Get("h_costheta_nocuts");
  h1.push_back(h_luminosity_cross_2f);
  TH1F * htemp[20];
  for(int j=0; j<histonames.size(); j++) {
    htemp[j]= (TH1F *)f->Get(histonames.at(j));
    h1.push_back(htemp[j]);
  }

  
  float cross_sec_temp=0;
  for(int i=0; i<cross_section_pairs.size(); i++) {
    if(cross_section_pairs.at(i).first==sample+"_"+pol) cross_sec_temp=cross_section_pairs.at(i).second;
  }
  if(cross_sec_temp==0) {
    cout<<" ERRORR !! NO CROSS SECTION FOUND FOR SAMPLE "<<sample<<" ... in GetHisto1D"<<endl;
    return h0;
  }

  float luminosity_sample = h_luminosity_cross_2f->GetEntries() / cross_sec_temp;
  for (int j = 0; j < h1.size(); j++) 
    h1.at(j)->Scale(lum / luminosity_sample);


  return h1;

}



std::vector<TH2F *> GetHisto2D(TString sample = "2f_hadronic_sample", float lum = 900, std::vector<TString> histonames={"h_nch","h_npfos","h_costheta_energy","h_mw1_mw2","h_major_minor_thrust"},bool correl_norm=true, TString pol="unpol")
{

  std::vector<TH2F *> h1;
  std::vector<TH2F *> h0;


  TString filename = folder + "_" + sample +  ".root";
  cout<<filename<<endl;
  TFile *f = new TFile(filename);
  TH1F *h_luminosity_cross_2f = (TH1F *)f->Get("h_costheta_nocuts");
  TH2F * htemp[20];
  for(int j=0; j<histonames.size(); j++) {
    cout<<j<<" "<<histonames.at(j)<<endl;
    htemp[j]= (TH2F *)f->Get(histonames.at(j));
    h1.push_back(htemp[j]);
  }

  float cross_sec_temp=0;
  for(int i=0; i<cross_section_pairs.size(); i++) {
    if(cross_section_pairs.at(i).first==sample+"_"+pol) cross_sec_temp=cross_section_pairs.at(i).second;
  }
  if(cross_sec_temp==0) {
    cout<<" ERRORR !! NO CROSS SECTION FOUND FOR SAMPLE"<<sample<<" ... in GetHisto2D"<<endl;
    return h0;
  }

  cout<<h1.size()<<" "<<endl;
  float luminosity_sample = h_luminosity_cross_2f->GetEntries() / cross_sec_temp;
  if(correl_norm==true) {
    for (int j = 0; j < h1.size(); j++) {
      // h1.at(j)->Scale(1./h1.at(j)->GetEntries());
      for(int i1=0; i1<h1.at(j)->GetNbinsX(); i1++) 
        for(int i2=0; i2<h1.at(j)->GetNbinsY(); i2++) {
          h1.at(j)->SetBinError(i1+1,i2+1,sqrt(h1.at(j)->GetBinContent(i1+1,i2+1)*luminosity_sample/lum)/h1.at(j)->GetEntries());///h1.at(j)->GetEntries()* luminosity_sample/lum);
          h1.at(j)->SetBinContent(i1+1,i2+1,h1.at(j)->GetBinContent(i1+1,i2+1)/h1.at(j)->GetEntries());// luminosity_sample);
        }
    }
  } else {
       for (int j = 0; j < h1.size(); j++) {
        h1.at(j)->Scale(lum/luminosity_sample);
    }
  }

  cout<<sample<<" lum_sample:"<<luminosity_sample<<"  lum:"<<lum<<endl;

  return h1;


}


std::vector<TH2F *> GetHisto2D_2(TString sample = "2f_hadronic_sample", float lum = 900, std::vector<TString> histonames={"h_nch","h_npfos","h_costheta_energy","h_mw1_mw2","h_major_minor_thrust"},bool correl_norm=true, TString pol="unpol")
{

  std::vector<TH2F *> h1;
  std::vector<TH2F *> h0;


  TString filename = folder + "_" + sample +  ".root";
  TFile *f = new TFile(filename);
  TH1F *h_luminosity_cross_2f = (TH1F *)f->Get("h_costheta_nocuts");
  TH2F * htemp[20];
  for(int j=0; j<histonames.size(); j++) {
    htemp[j]= (TH2F *)f->Get(histonames.at(j));
    h1.push_back(htemp[j]);
  }

  float cross_sec_temp=0;
  for(int i=0; i<cross_section_pairs.size(); i++) {
    if(cross_section_pairs.at(i).first==sample+"_"+pol) cross_sec_temp=cross_section_pairs.at(i).second;
  }
  if(cross_sec_temp==0) {
    cout<<" ERRORR !! NO CROSS SECTION FOUND FOR SAMPLE"<<sample<<" ... in GetHisto2D"<<endl;
    return h0;
  }

  float luminosity_sample = h_luminosity_cross_2f->GetEntries() / cross_sec_temp;
  if(correl_norm==true) {
    for (int j = 0; j < h1.size(); j++) {
      // h1.at(j)->Scale(1./h1.at(j)->GetEntries());
      for(int i1=0; i1<h1.at(j)->GetNbinsX(); i1++) 
        for(int i2=0; i2<h1.at(j)->GetNbinsY(); i2++) {
          h1.at(j)->SetBinError(i1+1,i2+1,2.*sqrt(h1.at(j)->GetBinContent(i1+1,i2+1)*luminosity_sample/lum)/h1.at(j)->GetEntries());///h1.at(j)->GetEntries()* luminosity_sample/lum);
          h1.at(j)->SetBinContent(i1+1,i2+1,h1.at(j)->GetBinContent(i1+1,i2+1)/h1.at(j)->GetEntries());// luminosity_sample);
        }
    }
  } else {
       for (int j = 0; j < h1.size(); j++) {
        h1.at(j)->Scale(lum/luminosity_sample);
    }
  }

  cout<<sample<<" lum_sample:"<<luminosity_sample<<"  lum:"<<lum<<endl;

  return h1;


}
