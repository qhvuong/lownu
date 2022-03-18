void rate()
{
  TChain * tree = new TChain( "cafTree", "cafTree" );
  tree->Add("/nashome/q/qvuong/gen_data/ND_FHC_FV_00.root");

  TChain * meta = new TChain( "meta", "meta" );
  meta->Add( "/nashome/q/qvuong/gen_data/ND_FHC_FV_00.root" ); // make certain this is the exact same file(s)
  double total_pot_CC = 0.;
  double pot_CC;
  double yrPOT_CC = 1.1E21;
  meta->SetBranchAddress( "pot", &pot_CC );
  //const int Nfiles = meta->GetEntries();
  for( int ii = 0; ii < meta->GetEntries(); ++ii ) {
    meta->GetEntry(ii);
    total_pot_CC += pot_CC;
  }

  TH1D *hRate_CC_e0 = new TH1D("hRate_CC_e0","",50,0,40);
  TH1D *hRate_CC_e3 = new TH1D("hRate_CC_e3","",50,0,40);
  TH1D *hRate_CC_m0 = new TH1D("hRate_CC_m0","",50,0,40);
  TH1D *hRate_CC_m3 = new TH1D("hRate_CC_m3","",50,0,40);

  double Ev;
  double eP, eN, ePip, ePim, ePi0, eOther;
  int LepPDG;
  tree->SetBranchAddress( "Ev", &Ev );
  tree->SetBranchAddress( "eP", &eP );
  tree->SetBranchAddress( "eN", &eN );
  tree->SetBranchAddress( "ePip", &ePip );
  tree->SetBranchAddress( "ePim", &ePim );
  tree->SetBranchAddress( "ePi0", &ePi0 );
  tree->SetBranchAddress( "eOther", &eOther );
  tree->SetBranchAddress( "LepPDG", &LepPDG );

  const int N = tree->GetEntries();
  for( int i = 0; i < N; i++){
    tree->GetEntry(i);
    double nu = eP + eN + ePip + ePim + ePi0 + eOther;

    if(LepPDG == 11){
      if(nu<10.0) hRate_CC_e0->Fill(Ev);
      if(nu<0.3)  hRate_CC_e3->Fill(Ev);
    }
    if(LepPDG == 13){
      if(nu<10.0) hRate_CC_m0->Fill(Ev);
      if(nu<0.3)  hRate_CC_m3->Fill(Ev);
    }
  }

  hRate_CC_e0->SetStats(0);
  hRate_CC_e0->Scale(yrPOT_CC/total_pot_CC);
  hRate_CC_e0->SetTitle("nueCC Rate vs Ev");
  hRate_CC_e0->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hRate_CC_e0->GetYaxis()->SetTitle("Rate (N/1.1#times 10^{21}POT)");
  hRate_CC_e3->SetStats(0);
  hRate_CC_e3->Scale(yrPOT_CC/total_pot_CC);
  hRate_CC_e3->SetTitle("nueCC Rate vs Ev");
  hRate_CC_e3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hRate_CC_e3->GetYaxis()->SetTitle("Rate (N/1.1#times 10^{21}POT)");
  hRate_CC_m0->SetStats(0);
  hRate_CC_m0->Scale(yrPOT_CC/total_pot_CC);
  hRate_CC_m0->SetTitle("numuCC Rate vs Ev");
  hRate_CC_m0->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hRate_CC_m0->GetYaxis()->SetTitle("Rate (N/1.1#times 10^{21}POT)");
  hRate_CC_m3->SetStats(0);
  hRate_CC_m3->Scale(yrPOT_CC/total_pot_CC);
  hRate_CC_m3->SetTitle("numuCC Rate vs Ev");
  hRate_CC_m3->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hRate_CC_m3->GetYaxis()->SetTitle("Rate (N/1.1#times 10^{21}POT)");



  TChain * tree_nue = new TChain( "tree", "tree" );
  tree_nue->Add("/nashome/q/qvuong/gen_data/nueFHC_000.root");

  TChain * meta_nue = new TChain( "meta", "meta" );
  meta_nue->Add( "/nashome/q/qvuong/gen_data/nueFHC_000.root" ); // make certain this is the exact same file(s)
  double total_pot_nue = 0.;
  double pot_nue;
  double yrPOT_nue = 3.3E21;
  meta_nue->SetBranchAddress( "pot", &pot_nue );
  //const int Nfiles = meta->GetEntries();
  for( int ii = 0; ii < meta_nue->GetEntries(); ++ii ) {
    meta_nue->GetEntry(ii);
    total_pot_nue += pot_nue;
  }

  TH1D *hRate_nue = new TH1D("hRate_nue","",50,0,40);

  double Enu;
  tree_nue->SetBranchAddress( "Enu", &Enu );

  const int N_nue = tree_nue->GetEntries();
  for( int i = 0; i < N_nue; i++){
    tree_nue->GetEntry(i);
    hRate_nue->Fill(Enu);
  }

  hRate_nue->SetStats(0);
  hRate_nue->Scale(yrPOT_nue/total_pot_nue);
  hRate_nue->SetTitle("nu+e Rate vs Ev");
  hRate_nue->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hRate_nue->GetYaxis()->SetTitle("Rate (N/3.3#times 10^{21}POT)");

 
  TCanvas *c1 = new TCanvas("c1","",800,600);
  gPad->SetLogy();
  hRate_CC_m0->SetLineColor(kRed);
  hRate_CC_m0->Draw();
  hRate_CC_m3->SetLineColor(kBlue);
  hRate_CC_m3->Draw("same");
  c1->SaveAs("RateVsEv_CC_m.pdf");

  TCanvas *c2 = new TCanvas("c2","",800,600);
  gPad->SetLogy();
  hRate_CC_e0->SetLineColor(kRed);
  hRate_CC_e0->Draw();
  hRate_CC_e3->SetLineColor(kBlue);
  hRate_CC_e3->Draw("same");
  c2->SaveAs("RateVsEv_CC_e.pdf");

  TCanvas *c3 = new TCanvas("c3","",800,600);
  gPad->SetLogy();
  hRate_nue->SetLineColor(kRed);
  hRate_nue->Draw();
  c3->SaveAs("RateVsEv_nue.pdf");
}
