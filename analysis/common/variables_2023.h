// Declaration of leaf types
   Float_t         mc_quark_E[2];
   Float_t         mc_quark_px[2];
   Float_t         mc_quark_py[2];
   Float_t         mc_quark_pz[2];
   Float_t         mc_quark_m[2];
   Float_t         mc_quark_pdg[2];
   Float_t         mc_quark_charge[2];
   Float_t         mc_ISR_E[2];
   Float_t         mc_ISR_px[2];
   Float_t         mc_ISR_py[2];
   Float_t         mc_ISR_pz[2];
   Float_t         mc_ISR_m[2];
   Float_t         mc_ISR_pdg[2];
   Float_t         mc_ISR_charge[2];
   Int_t           mc_quark_ps_n;
   Float_t         mc_quark_ps_E[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_px[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_py[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_pz[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_m[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_pdg[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_charge[1];   //[mc_quark_ps_n]
   Float_t         mc_quark_ps_y12;
   Float_t         mc_quark_ps_y23;
   Float_t         mc_quark_ps_d12;
   Float_t         mc_quark_ps_d23;
   Float_t         mc_quark_ps_jet_E[2];
   Float_t         mc_quark_ps_jet_px[2];
   Float_t         mc_quark_ps_jet_py[2];
   Float_t         mc_quark_ps_jet_pz[2];
   Int_t           mc_stable_n;
   Float_t         mc_stable_E[1];   //[mc_stable_n]
   Float_t         mc_stable_px[1];   //[mc_stable_n]
   Float_t         mc_stable_py[1];   //[mc_stable_n]
   Float_t         mc_stable_pz[1];   //[mc_stable_n]
   Float_t         mc_stable_m[1];   //[mc_stable_n]
   Int_t           mc_stable_pdg[1];   //[mc_stable_n]
   Float_t         mc_stable_charge[1];   //[mc_stable_n]
   Int_t           mc_stable_isoverlay[1];   //[mc_stable_n]
   Int_t           mc_stable_isisr[1];   //[mc_stable_n]
   Float_t         mc_stable_y12;
   Float_t         mc_stable_y23;
   Float_t         mc_stable_d12;
   Float_t         mc_stable_d23;
   Float_t         mc_stable_jet_E[2];
   Float_t         mc_stable_jet_px[2];
   Float_t         mc_stable_jet_py[2];
   Float_t         mc_stable_jet_pz[2];
   Float_t         mc_major_thrust_value;
   Float_t         mc_major_thrust_axis[3];
   Float_t         mc_minor_thrust_value;
   Float_t         mc_minor_thrust_axis[3];
   Float_t         mc_principle_thrust_value;
   Float_t         mc_principle_thrust_axis[3];
   Float_t         truejet_E[5];
   Float_t         truejet_px[5];
   Float_t         truejet_py[5];
   Float_t         truejet_pz[5];
   Int_t           truejet_type[5];
   Int_t           truejet_pdg[5];
   Float_t         jet_E[2];
   Float_t         jet_px[2];
   Float_t         jet_py[2];
   Float_t         jet_pz[2];
   Float_t         jet_btag[2];
   Float_t         jet_ctag[2];
   Float_t         jet_btag_2[2];
   Float_t         jet_ctag_2[2];
   Float_t         y23;
   Float_t         y12;
   Float_t         d23;
   Float_t         d12;
   Float_t         oblateness;
   Float_t         aplanarity;
   Float_t         major_thrust_value;
   Float_t         major_thrust_axis[3];
   Float_t         minor_thrust_value;
   Float_t         minor_thrust_axis[3];
   Float_t         principle_thrust_value;
   Float_t         principle_thrust_axis[3];
   Float_t         sphericity;
   Float_t         sphericity_tensor[3];
   Float_t         fourjet_E[4];
   Float_t         fourjet_px[4];
   Float_t         fourjet_py[4];
   Float_t         fourjet_pz[4];
   Int_t           pfo_n;
   Int_t           jet_nvtx;
   Int_t           pfo_n_j1;
   Int_t           jet_nvtx_j1;
   Int_t           pfo_n_j2;
   Int_t           jet_nvtx_j2;
   Int_t           pfo_match[2000];   //[pfo_n]
   Int_t           pfo_truejet_pdg[2000];   //[pfo_n]
   Int_t           pfo_truejet_type[2000];   //[pfo_n]
   Int_t           pfo_pdgcheat[2000];   //[pfo_n]
   Float_t         pfo_E[2000];   //[pfo_n]
   Float_t         pfo_px[2000];   //[pfo_n]
   Float_t         pfo_py[2000];   //[pfo_n]
   Float_t         pfo_pz[2000];   //[pfo_n]
   Float_t         pfo_m[2000];   //[pfo_n]
   Int_t           pfo_type[2000];   //[pfo_n]
   Int_t           pfo_isoverlay[2000];   //[pfo_n]
   Int_t           pfo_isisr[2000];   //[pfo_n]
   Int_t           pfo_vtx[2000];   //[pfo_n]
   Int_t           pfo_charge[2000];   //[pfo_n]
   Int_t           pfo_ntracks[2000];   //[pfo_n]
   Int_t           pfo_tpc_hits[2000];   //[pfo_n]
   Float_t         pfo_dedx[2000];   //[pfo_n]
   Float_t         pfo_dedxerror[2000];   //[pfo_n]
   Float_t         pfo_d0[2000];   //[pfo_n]
   Float_t         pfo_d0error[2000];   //[pfo_n]
   Float_t         pfo_z0[2000];   //[pfo_n]
   Float_t         pfo_z0error[2000];   //[pfo_n]
   Float_t         pfo_phi[2000];   //[pfo_n]
   Float_t         pfo_phierror[2000];   //[pfo_n]
   Float_t         pfo_omega[2000];   //[pfo_n]
   Float_t         pfo_omegaerror[2000];   //[pfo_n]
   Float_t         pfo_tanlambda[2000];   //[pfo_n]
   Float_t         pfo_tanlambdaerror[2000];   //[pfo_n]
   Float_t         pfo_chi2[2000];   //[pfo_n]
   Float_t         pfo_ndf[2000];   //[pfo_n]
   Float_t         pfo_piddedx_e_dedxdist[2000];   //[pfo_n]
   Float_t         pfo_piddedx_mu_dedxdist[2000];   //[pfo_n]
   Float_t         pfo_piddedx_pi_dedxdist[2000];   //[pfo_n]
   Float_t         pfo_piddedx_k_dedxdist[2000];   //[pfo_n]
   Float_t         pfo_piddedx_p_dedxdist[2000];   //[pfo_n]
   Float_t         pfo_piddedx_e_dedxdist_2[2000];   //[pfo_n]
   Float_t         pfo_piddedx_mu_dedxdist_2[2000];   //[pfo_n]
   Float_t         pfo_piddedx_pi_dedxdist_2[2000];   //[pfo_n]
   Float_t         pfo_piddedx_k_dedxdist_2[2000];   //[pfo_n]
   Float_t         pfo_piddedx_p_dedxdist_2[2000];   //[pfo_n]

   // List of branches
   TBranch        *b_mc_quark_E;   //!
   TBranch        *b_mc_quark_px;   //!
   TBranch        *b_mc_quark_py;   //!
   TBranch        *b_mc_quark_pz;   //!
   TBranch        *b_mc_quark_m;   //!
   TBranch        *b_mc_quark_pdg;   //!
   TBranch        *b_mc_quark_charge;   //!
   TBranch        *b_mc_ISR_E;   //!
   TBranch        *b_mc_ISR_px;   //!
   TBranch        *b_mc_ISR_py;   //!
   TBranch        *b_mc_ISR_pz;   //!
   TBranch        *b_mc_ISR_m;   //!
   TBranch        *b_mc_ISR_pdg;   //!
   TBranch        *b_mc_ISR_charge;   //!
   TBranch        *b_mc_quark_ps_n;   //!
   TBranch        *b_mc_quark_ps_E;   //!
   TBranch        *b_mc_quark_ps_px;   //!
   TBranch        *b_mc_quark_ps_py;   //!
   TBranch        *b_mc_quark_ps_pz;   //!
   TBranch        *b_mc_quark_ps_m;   //!
   TBranch        *b_mc_quark_ps_pdg;   //!
   TBranch        *b_mc_quark_ps_charge;   //!
   TBranch        *b_mc_quark_ps_y12;   //!
   TBranch        *b_mc_quark_ps_y23;   //!
   TBranch        *b_mc_quark_ps_d12;   //!
   TBranch        *b_mc_quark_ps_d23;   //!
   TBranch        *b_mc_quark_ps_jet_E;   //!
   TBranch        *b_mc_quark_ps_jet_px;   //!
   TBranch        *b_mc_quark_ps_jet_py;   //!
   TBranch        *b_mc_quark_ps_jet_pz;   //!
   TBranch        *b_mc_stable_n;   //!
   TBranch        *b_mc_stable_E;   //!
   TBranch        *b_mc_stable_px;   //!
   TBranch        *b_mc_stable_py;   //!
   TBranch        *b_mc_stable_pz;   //!
   TBranch        *b_mc_stable_m;   //!
   TBranch        *b_mc_stable_pdg;   //!
   TBranch        *b_mc_stable_charge;   //!
   TBranch        *b_mc_stable_isoverlay;   //!
   TBranch        *b_mc_stable_isisr;   //!
   TBranch        *b_mc_stable_y12;   //!
   TBranch        *b_mc_stable_y23;   //!
   TBranch        *b_mc_stable_d12;   //!
   TBranch        *b_mc_stable_d23;   //!
   TBranch        *b_mc_stable_jet_E;   //!
   TBranch        *b_mc_stable_jet_px;   //!
   TBranch        *b_mc_stable_jet_py;   //!
   TBranch        *b_mc_stable_jet_pz;   //!
   TBranch        *b_mc_major_thrust_value;   //!
   TBranch        *b_mc_major_thrust_axis;   //!
   TBranch        *b_mc_minor_thrust_value;   //!
   TBranch        *b_mc_minor_thrust_axis;   //!
   TBranch        *b_mc_principle_thrust_value;   //!
   TBranch        *b_mc_principle_thrust_axis;   //!
   TBranch        *b_truejet_E;   //!
   TBranch        *b_truejet_px;   //!
   TBranch        *b_truejet_py;   //!
   TBranch        *b_truejet_pz;   //!
   TBranch        *b_truejet_type;   //!
   TBranch        *b_truejet_pdg;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_btag;   //!
   TBranch        *b_jet_ctag;   //!
   TBranch        *b_jet_btag_2;   //!
   TBranch        *b_jet_ctag_2;   //!
   TBranch        *b_y23;   //!
   TBranch        *b_y12;   //!
   TBranch        *b_d23;   //!
   TBranch        *b_d12;   //!
   TBranch        *b_oblateness;   //!
   TBranch        *b_aplanarity;   //!
   TBranch        *b_major_thrust_value;   //!
   TBranch        *b_major_thrust_axis;   //!
   TBranch        *b_minor_thrust_value;   //!
   TBranch        *b_minor_thrust_axis;   //!
   TBranch        *b_principle_thrust_value;   //!
   TBranch        *b_principle_thrust_axis;   //!
   TBranch        *b_sphericity;   //!
   TBranch        *b_sphericity_tensor;   //!
   TBranch        *b_fourjet_E;   //!
   TBranch        *b_fourjet_px;   //!
   TBranch        *b_fourjet_py;   //!
   TBranch        *b_fourjet_pz;   //!
   TBranch        *b_pfo_n;   //!
   TBranch        *b_jet_nvtx;   //!
   TBranch        *b_pfo_n_j1;   //!
   TBranch        *b_jet_nvtx_j1;   //!
   TBranch        *b_pfo_n_j2;   //!
   TBranch        *b_jet_nvtx_j2;   //!
   TBranch        *b_pfo_match;   //!
   TBranch        *b_pfo_truejet_pdg;   //!
   TBranch        *b_pfo_truejet_type;   //!
   TBranch        *b_pfo_pdgcheat;   //!
   TBranch        *b_pfo_E;   //!
   TBranch        *b_pfo_px;   //!
   TBranch        *b_pfo_py;   //!
   TBranch        *b_pfo_pz;   //!
   TBranch        *b_pfo_m;   //!
   TBranch        *b_pfo_type;   //!
   TBranch        *b_pfo_isoverlay;   //!
   TBranch        *b_pfo_isisr;   //!
   TBranch        *b_pfo_vtx;   //!
   TBranch        *b_pfo_charge;   //!
   TBranch        *b_pfo_ntracks;   //!
   TBranch        *b_pfo_tpc_hits;   //!
   TBranch        *b_pfo_dedx;   //!
   TBranch        *b_pfo_dedxerror;   //!
   TBranch        *b_pfo_d0;   //!
   TBranch        *b_pfo_d0error;   //!
   TBranch        *b_pfo_z0;   //!
   TBranch        *b_pfo_z0error;   //!
   TBranch        *b_pfo_phi;   //!
   TBranch        *b_pfo_phierror;   //!
   TBranch        *b_pfo_omega;   //!
   TBranch        *b_pfo_omegaerror;   //!
   TBranch        *b_pfo_tanlambda;   //!
   TBranch        *b_pfo_tanlambdaerror;   //!
   TBranch        *b_pfo_chi2;   //!
   TBranch        *b_pfo_ndf;   //!
   TBranch        *b_pfo_piddedx_e_dedxdist;   //!
   TBranch        *b_pfo_piddedx_mu_dedxdist;   //!
   TBranch        *b_pfo_piddedx_pi_dedxdist;   //!
   TBranch        *b_pfo_piddedx_k_dedxdist;   //!
   TBranch        *b_pfo_piddedx_p_dedxdist;   //!
   TBranch        *b_pfo_piddedx_e_dedxdist_2;   //!
   TBranch        *b_pfo_piddedx_mu_dedxdist_2;   //!
   TBranch        *b_pfo_piddedx_pi_dedxdist_2;   //!
   TBranch        *b_pfo_piddedx_k_dedxdist_2;   //!
   TBranch        *b_pfo_piddedx_p_dedxdist_2;   //!
   
   


void SetBranchAddressFunction(TTree *fChain)
{

  fChain->SetBranchAddress("mc_quark_E", mc_quark_E, &b_mc_quark_E);
   fChain->SetBranchAddress("mc_quark_px", mc_quark_px, &b_mc_quark_px);
   fChain->SetBranchAddress("mc_quark_py", mc_quark_py, &b_mc_quark_py);
   fChain->SetBranchAddress("mc_quark_pz", mc_quark_pz, &b_mc_quark_pz);
   fChain->SetBranchAddress("mc_quark_m", mc_quark_m, &b_mc_quark_m);
   fChain->SetBranchAddress("mc_quark_pdg", mc_quark_pdg, &b_mc_quark_pdg);
   fChain->SetBranchAddress("mc_quark_charge", mc_quark_charge, &b_mc_quark_charge);
   fChain->SetBranchAddress("mc_ISR_E", mc_ISR_E, &b_mc_ISR_E);
   fChain->SetBranchAddress("mc_ISR_px", mc_ISR_px, &b_mc_ISR_px);
   fChain->SetBranchAddress("mc_ISR_py", mc_ISR_py, &b_mc_ISR_py);
   fChain->SetBranchAddress("mc_ISR_pz", mc_ISR_pz, &b_mc_ISR_pz);
   fChain->SetBranchAddress("mc_ISR_m", mc_ISR_m, &b_mc_ISR_m);
   fChain->SetBranchAddress("mc_ISR_pdg", mc_ISR_pdg, &b_mc_ISR_pdg);
   fChain->SetBranchAddress("mc_ISR_charge", mc_ISR_charge, &b_mc_ISR_charge);
   fChain->SetBranchAddress("mc_quark_ps_n", &mc_quark_ps_n, &b_mc_quark_ps_n);
   fChain->SetBranchAddress("mc_quark_ps_E", &mc_quark_ps_E, &b_mc_quark_ps_E);
   fChain->SetBranchAddress("mc_quark_ps_px", &mc_quark_ps_px, &b_mc_quark_ps_px);
   fChain->SetBranchAddress("mc_quark_ps_py", &mc_quark_ps_py, &b_mc_quark_ps_py);
   fChain->SetBranchAddress("mc_quark_ps_pz", &mc_quark_ps_pz, &b_mc_quark_ps_pz);
   fChain->SetBranchAddress("mc_quark_ps_m", &mc_quark_ps_m, &b_mc_quark_ps_m);
   fChain->SetBranchAddress("mc_quark_ps_pdg", &mc_quark_ps_pdg, &b_mc_quark_ps_pdg);
   fChain->SetBranchAddress("mc_quark_ps_charge", &mc_quark_ps_charge, &b_mc_quark_ps_charge);
   fChain->SetBranchAddress("mc_quark_ps_y12", &mc_quark_ps_y12, &b_mc_quark_ps_y12);
   fChain->SetBranchAddress("mc_quark_ps_y23", &mc_quark_ps_y23, &b_mc_quark_ps_y23);
   fChain->SetBranchAddress("mc_quark_ps_d12", &mc_quark_ps_d12, &b_mc_quark_ps_d12);
   fChain->SetBranchAddress("mc_quark_ps_d23", &mc_quark_ps_d23, &b_mc_quark_ps_d23);
   fChain->SetBranchAddress("mc_quark_ps_jet_E", mc_quark_ps_jet_E, &b_mc_quark_ps_jet_E);
   fChain->SetBranchAddress("mc_quark_ps_jet_px", mc_quark_ps_jet_px, &b_mc_quark_ps_jet_px);
   fChain->SetBranchAddress("mc_quark_ps_jet_py", mc_quark_ps_jet_py, &b_mc_quark_ps_jet_py);
   fChain->SetBranchAddress("mc_quark_ps_jet_pz", mc_quark_ps_jet_pz, &b_mc_quark_ps_jet_pz);
   fChain->SetBranchAddress("mc_stable_n", &mc_stable_n, &b_mc_stable_n);
   fChain->SetBranchAddress("mc_stable_E", &mc_stable_E, &b_mc_stable_E);
   fChain->SetBranchAddress("mc_stable_px", &mc_stable_px, &b_mc_stable_px);
   fChain->SetBranchAddress("mc_stable_py", &mc_stable_py, &b_mc_stable_py);
   fChain->SetBranchAddress("mc_stable_pz", &mc_stable_pz, &b_mc_stable_pz);
   fChain->SetBranchAddress("mc_stable_m", &mc_stable_m, &b_mc_stable_m);
   fChain->SetBranchAddress("mc_stable_pdg", &mc_stable_pdg, &b_mc_stable_pdg);
   fChain->SetBranchAddress("mc_stable_charge", &mc_stable_charge, &b_mc_stable_charge);
   fChain->SetBranchAddress("mc_stable_isoverlay", &mc_stable_isoverlay, &b_mc_stable_isoverlay);
   fChain->SetBranchAddress("mc_stable_isisr", &mc_stable_isisr, &b_mc_stable_isisr);
   fChain->SetBranchAddress("mc_stable_y12", &mc_stable_y12, &b_mc_stable_y12);
   fChain->SetBranchAddress("mc_stable_y23", &mc_stable_y23, &b_mc_stable_y23);
   fChain->SetBranchAddress("mc_stable_d12", &mc_stable_d12, &b_mc_stable_d12);
   fChain->SetBranchAddress("mc_stable_d23", &mc_stable_d23, &b_mc_stable_d23);
   fChain->SetBranchAddress("mc_stable_jet_E", mc_stable_jet_E, &b_mc_stable_jet_E);
   fChain->SetBranchAddress("mc_stable_jet_px", mc_stable_jet_px, &b_mc_stable_jet_px);
   fChain->SetBranchAddress("mc_stable_jet_py", mc_stable_jet_py, &b_mc_stable_jet_py);
   fChain->SetBranchAddress("mc_stable_jet_pz", mc_stable_jet_pz, &b_mc_stable_jet_pz);
   fChain->SetBranchAddress("mc_major_thrust_value", &mc_major_thrust_value, &b_mc_major_thrust_value);
   fChain->SetBranchAddress("mc_major_thrust_axis", mc_major_thrust_axis, &b_mc_major_thrust_axis);
   fChain->SetBranchAddress("mc_minor_thrust_value", &mc_minor_thrust_value, &b_mc_minor_thrust_value);
   fChain->SetBranchAddress("mc_minor_thrust_axis", mc_minor_thrust_axis, &b_mc_minor_thrust_axis);
   fChain->SetBranchAddress("mc_principle_thrust_value", &mc_principle_thrust_value, &b_mc_principle_thrust_value);
   fChain->SetBranchAddress("mc_principle_thrust_axis", mc_principle_thrust_axis, &b_mc_principle_thrust_axis);
   fChain->SetBranchAddress("truejet_E", truejet_E, &b_truejet_E);
   fChain->SetBranchAddress("truejet_px", truejet_px, &b_truejet_px);
   fChain->SetBranchAddress("truejet_py", truejet_py, &b_truejet_py);
   fChain->SetBranchAddress("truejet_pz", truejet_pz, &b_truejet_pz);
   fChain->SetBranchAddress("truejet_type", truejet_type, &b_truejet_type);
   fChain->SetBranchAddress("truejet_pdg", truejet_pdg, &b_truejet_pdg);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_px", jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet_py", jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet_pz", jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet_btag", jet_btag, &b_jet_btag);
   fChain->SetBranchAddress("jet_ctag", jet_ctag, &b_jet_ctag);
   fChain->SetBranchAddress("jet_btag_2", jet_btag_2, &b_jet_btag_2);
   fChain->SetBranchAddress("jet_ctag_2", jet_ctag_2, &b_jet_ctag_2);
   fChain->SetBranchAddress("y23", &y23, &b_y23);
   fChain->SetBranchAddress("y12", &y12, &b_y12);
   fChain->SetBranchAddress("d23", &d23, &b_d23);
   fChain->SetBranchAddress("d12", &d12, &b_d12);
   fChain->SetBranchAddress("oblateness", &oblateness, &b_oblateness);
   fChain->SetBranchAddress("aplanarity", &aplanarity, &b_aplanarity);
   fChain->SetBranchAddress("major_thrust_value", &major_thrust_value, &b_major_thrust_value);
   fChain->SetBranchAddress("major_thrust_axis", major_thrust_axis, &b_major_thrust_axis);
   fChain->SetBranchAddress("minor_thrust_value", &minor_thrust_value, &b_minor_thrust_value);
   fChain->SetBranchAddress("minor_thrust_axis", minor_thrust_axis, &b_minor_thrust_axis);
   fChain->SetBranchAddress("principle_thrust_value", &principle_thrust_value, &b_principle_thrust_value);
   fChain->SetBranchAddress("principle_thrust_axis", principle_thrust_axis, &b_principle_thrust_axis);
   fChain->SetBranchAddress("sphericity", &sphericity, &b_sphericity);
   fChain->SetBranchAddress("sphericity_tensor", sphericity_tensor, &b_sphericity_tensor);
   fChain->SetBranchAddress("fourjet_E", fourjet_E, &b_fourjet_E);
   fChain->SetBranchAddress("fourjet_px", fourjet_px, &b_fourjet_px);
   fChain->SetBranchAddress("fourjet_py", fourjet_py, &b_fourjet_py);
   fChain->SetBranchAddress("fourjet_pz", fourjet_pz, &b_fourjet_pz);
   fChain->SetBranchAddress("pfo_n", &pfo_n, &b_pfo_n);
   fChain->SetBranchAddress("jet_nvtx", &jet_nvtx, &b_jet_nvtx);
   fChain->SetBranchAddress("pfo_n_j1", &pfo_n_j1, &b_pfo_n_j1);
   fChain->SetBranchAddress("jet_nvtx_j1", &jet_nvtx_j1, &b_jet_nvtx_j1);
   fChain->SetBranchAddress("pfo_n_j2", &pfo_n_j2, &b_pfo_n_j2);
   fChain->SetBranchAddress("jet_nvtx_j2", &jet_nvtx_j2, &b_jet_nvtx_j2);
   fChain->SetBranchAddress("pfo_match", pfo_match, &b_pfo_match);
   fChain->SetBranchAddress("pfo_truejet_pdg", pfo_truejet_pdg, &b_pfo_truejet_pdg);
   fChain->SetBranchAddress("pfo_truejet_type", pfo_truejet_type, &b_pfo_truejet_type);
   fChain->SetBranchAddress("pfo_pdgcheat", pfo_pdgcheat, &b_pfo_pdgcheat);
   fChain->SetBranchAddress("pfo_E", pfo_E, &b_pfo_E);
   fChain->SetBranchAddress("pfo_px", pfo_px, &b_pfo_px);
   fChain->SetBranchAddress("pfo_py", pfo_py, &b_pfo_py);
   fChain->SetBranchAddress("pfo_pz", pfo_pz, &b_pfo_pz);
   fChain->SetBranchAddress("pfo_m", pfo_m, &b_pfo_m);
   fChain->SetBranchAddress("pfo_type", pfo_type, &b_pfo_type);
   fChain->SetBranchAddress("pfo_isoverlay", pfo_isoverlay, &b_pfo_isoverlay);
   fChain->SetBranchAddress("pfo_isisr", pfo_isisr, &b_pfo_isisr);
   fChain->SetBranchAddress("pfo_vtx", pfo_vtx, &b_pfo_vtx);
   fChain->SetBranchAddress("pfo_charge", pfo_charge, &b_pfo_charge);
   fChain->SetBranchAddress("pfo_ntracks", pfo_ntracks, &b_pfo_ntracks);
   fChain->SetBranchAddress("pfo_tpc_hits", pfo_tpc_hits, &b_pfo_tpc_hits);
   fChain->SetBranchAddress("pfo_dedx", pfo_dedx, &b_pfo_dedx);
   fChain->SetBranchAddress("pfo_dedxerror", pfo_dedxerror, &b_pfo_dedxerror);
   fChain->SetBranchAddress("pfo_d0", pfo_d0, &b_pfo_d0);
   fChain->SetBranchAddress("pfo_d0error", pfo_d0error, &b_pfo_d0error);
   fChain->SetBranchAddress("pfo_z0", pfo_z0, &b_pfo_z0);
   fChain->SetBranchAddress("pfo_z0error", pfo_z0error, &b_pfo_z0error);
   fChain->SetBranchAddress("pfo_phi", pfo_phi, &b_pfo_phi);
   fChain->SetBranchAddress("pfo_phierror", pfo_phierror, &b_pfo_phierror);
   fChain->SetBranchAddress("pfo_omega", pfo_omega, &b_pfo_omega);
   fChain->SetBranchAddress("pfo_omegaerror", pfo_omegaerror, &b_pfo_omegaerror);
   fChain->SetBranchAddress("pfo_tanlambda", pfo_tanlambda, &b_pfo_tanlambda);
   fChain->SetBranchAddress("pfo_tanlambdaerror", pfo_tanlambdaerror, &b_pfo_tanlambdaerror);
   fChain->SetBranchAddress("pfo_chi2", pfo_chi2, &b_pfo_chi2);
   fChain->SetBranchAddress("pfo_ndf", pfo_ndf, &b_pfo_ndf);
   fChain->SetBranchAddress("pfo_piddedx_e_dedxdist", pfo_piddedx_e_dedxdist, &b_pfo_piddedx_e_dedxdist);
   fChain->SetBranchAddress("pfo_piddedx_mu_dedxdist", pfo_piddedx_mu_dedxdist, &b_pfo_piddedx_mu_dedxdist);
   fChain->SetBranchAddress("pfo_piddedx_pi_dedxdist", pfo_piddedx_pi_dedxdist, &b_pfo_piddedx_pi_dedxdist);
   fChain->SetBranchAddress("pfo_piddedx_k_dedxdist", pfo_piddedx_k_dedxdist, &b_pfo_piddedx_k_dedxdist);
   fChain->SetBranchAddress("pfo_piddedx_p_dedxdist", pfo_piddedx_p_dedxdist, &b_pfo_piddedx_p_dedxdist);
   fChain->SetBranchAddress("pfo_piddedx_e_dedxdist_2", pfo_piddedx_e_dedxdist_2, &b_pfo_piddedx_e_dedxdist_2);
   fChain->SetBranchAddress("pfo_piddedx_mu_dedxdist_2", pfo_piddedx_mu_dedxdist_2, &b_pfo_piddedx_mu_dedxdist_2);
   fChain->SetBranchAddress("pfo_piddedx_pi_dedxdist_2", pfo_piddedx_pi_dedxdist_2, &b_pfo_piddedx_pi_dedxdist_2);
   fChain->SetBranchAddress("pfo_piddedx_k_dedxdist_2", pfo_piddedx_k_dedxdist_2, &b_pfo_piddedx_k_dedxdist_2);
   fChain->SetBranchAddress("pfo_piddedx_p_dedxdist_2", pfo_piddedx_p_dedxdist_2, &b_pfo_piddedx_p_dedxdist_2);
   

}
