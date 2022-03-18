#include <string>

void lownu_prob_Elep()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  
  for(int i = 10; i<30; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  std::cout << "File number:" << i << "\n";
  } 
  // replace this with the location of your file
  // Note: this step is going to print a bunch of Error failure loading library. These are not actually errors, just ROOT being stupid


  // Get the number of "protons on target", which is essentially the number of neutrinos in this file
/*  TChain * meta = new TChain( "meta", "meta" );
  meta->Add( "/Users/qhuong/Documents/UofR/2021_Summer/data/ND_FHC_FV_00.root" ); // make certain this is the exact same file(s)
*/
  TH2D *pmueVsEv0 = new TH2D("pmueVsEv0","",100,0,40,100,0,0.0005);
  TH2D *pmueVsEv1 = new TH2D("pmueVsEv1","",100,0,40,100,0,0.015);
  TH2D *pmumuVsEv0 = new TH2D("pmumuVsEv0","",100,0,40,100,0.95,1);
  TH2D *pmumuVsEv1 = new TH2D("pmumuVsEv1","",100,0,40,100,0.75,1);
  TH2D *peeVsEv0 = new TH2D("peeVsEv0","",100,0,40,100,0.95,1);
  TH2D *peeVsEv1 = new TH2D("peeVsEv1","",100,0,40,100,0.75,1);

  TH1D *nmAcc0 = new TH1D("nmAcc0","",200,0,40); //numu acceptance
  TH1D *nmAcc3 = new TH1D("nmAcc3","",200,0,40); //numu acceptance
  TH1D *neAcc0 = new TH1D("neAcc0","",200,0,40); //nue acceptance
  TH1D *neAcc3 = new TH1D("neAcc3","",200,0,40); //nue acceptance
  TH1D *nmElep0 = new TH1D("nmEv0","",200,0,40); //numu acceptance
  TH1D *nmElep3 = new TH1D("nmEv3","",200,0,40); //numu acceptance
  TH1D *neElep0 = new TH1D("neEv0","",200,0,40); //nue acceptance
  TH1D *neElep3 = new TH1D("neEv3","",200,0,40); //nue acceptance

  TH1D *nmAcc0_1 = new TH1D("nmAcc0_1","",200,0,40); //numu acceptance
  TH1D *nmAcc3_1 = new TH1D("nmAcc3_1","",200,0,40); //numu acceptance
  TH1D *neAcc0_1 = new TH1D("neAcc0_1","",200,0,40); //nue acceptance
  TH1D *neAcc3_1 = new TH1D("neAcc3_1","",200,0,40); //nue acceptance

  //TH2D *nm_hElepVsEv0 = new TH2D("nm_hElepVsEv0","",100,0,10,100,0,10);
  //TH2D *nm_hElepVsEv3 = new TH2D("nm_hElepVsEv3","",100,0,10,100,0,10);
  // Most of them are weights related to systematic uncertainties
  
  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  int nuPDG; // PDG code of the neutrino, numu = 14, nue = 12, antineutrinos are negative
  double Ev; // the energy of the neutrino, in GeV
  int LepPDG; // PDG code of the final-state lepton, mu = 13, e = 11
  double LepE; // Total energy of the final-state lepton; note that true nu is not saved but nu = Ev - LepE
  double LepNuAngle; // angle between lepton and neutrino
  double eP, eN, ePip, ePim, ePi0, eOther; // energy in the final state due to different particles

  // information about the "reconstruction". We can talk more about what "reconstruction" means for this
  double Ev_reco, Elep_reco; // the measured neutrino energy and lepton energy. Reco nu = Ev_reco - Elep_reco
  int reco_numu, reco_nue; // = 1 if the reconstruction thinks it's a muon or electron
  int muon_contained, muon_tracker, muon_ecal; // different ways that the muon can be measured
  double Ehad_veto; // hadronic energy near the edge of the detector, which is a hint that we might not have measured all the energy

  
  tree->SetBranchAddress( "vtx_x", &vtx_x );
  tree->SetBranchAddress( "vtx_y", &vtx_y );
  tree->SetBranchAddress( "vtx_z", &vtx_z );
  tree->SetBranchAddress( "Ev", &Ev );
  tree->SetBranchAddress( "nuPDG", &nuPDG );
  tree->SetBranchAddress( "LepPDG", &LepPDG );
  tree->SetBranchAddress( "LepE", &LepE );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );
  tree->SetBranchAddress( "Ev_reco", &Ev_reco );
  tree->SetBranchAddress( "Elep_reco", &Elep_reco );
  tree->SetBranchAddress( "reco_numu", &reco_numu );
  tree->SetBranchAddress( "reco_nue", &reco_nue );
  tree->SetBranchAddress( "muon_contained", &muon_contained );
  tree->SetBranchAddress( "muon_tracker", &muon_tracker );
  tree->SetBranchAddress( "muon_ecal", &muon_ecal );
  tree->SetBranchAddress( "Ehad_veto", &Ehad_veto );
  tree->SetBranchAddress( "eP", &eP );
  tree->SetBranchAddress( "eN", &eN );
  tree->SetBranchAddress( "ePip", &ePip );
  tree->SetBranchAddress( "ePim", &ePim );
  tree->SetBranchAddress( "ePi0", &ePi0 );
  tree->SetBranchAddress( "eOther", &eOther );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );

  const int N = tree->GetEntries();
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    double nu = eP + eN + ePip + ePim + ePi0 + eOther;
    double reco_nu = Ev_reco - Elep_reco;

    if(LepPDG == 13){
      if(nu<10)  {nmElep0->Fill(LepE);}
      if(nu<0.3) {nmElep3->Fill(LepE);}}
    if(LepPDG == 11){
      if(nu<10)  {neElep0->Fill(LepE);}
      if(nu<0.3) {neElep3->Fill(LepE);}}

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
/*   
    if(LepPDG == 13){
      if(nu<10)  {nmElep0->Fill(LepE);}
      if(nu<0.3) {nmElep3->Fill(LepE);}}
    if(LepPDG == 11){
      if(nu<10)  {neElep0->Fill(LepE);}
      if(nu<0.3) {neElep3->Fill(LepE);}}
*/
    if(LepPDG == 13){
      //if(nu<10)  {nmElep0->Fill(LepE);}
      //if(nu<0.3) {nmElep3->Fill(LepE);}

      if(reco_numu && (muon_contained || muon_tracker || muon_ecal)){
        if(nu<10)  {nmAcc0->Fill(LepE);}
        if(nu<0.3) {nmAcc3->Fill(LepE);}
        if(reco_nu<10)  {nmAcc0_1->Fill(LepE);}
        if(reco_nu<0.3) {nmAcc3_1->Fill(LepE);}
      }
    }
    if(LepPDG == 11){
      //if(nu<10)  {neElep0->Fill(LepE);}
      //if(nu<0.3) {neElep3->Fill(LepE);}

      if(reco_nue){
        if(nu<10)  {neAcc0->Fill(LepE);}
        if(nu<0.3) {neAcc3->Fill(LepE);}
        if(reco_nu<10)  {neAcc0_1->Fill(LepE);}
        if(reco_nu<0.3) {neAcc3_1->Fill(LepE);}
      }
    }

  }


  nmAcc0->Divide(nmElep0);
  nmAcc0->SetTitle("numu acceptance (fiducial && recontruction && true_nu / true_nu)");
  nmAcc0->GetXaxis()->SetTitle("LepE (GeV)");
  //nmAcc0->GetYaxis()->SetTitle("reco_cut/true_cut");
  nmAcc3->Divide(nmElep3);
  neAcc0->Divide(neElep0);
  neAcc0->SetTitle("nue acceptance (fiducial && recontruction && true_nu / true_nu)");
  neAcc0->GetXaxis()->SetTitle("LepE (GeV)");
  //neAcc0->GetYaxis()->SetTitle("reco_cut/true_cut");
  neAcc3->Divide(neElep3);
  nmAcc0_1->Divide(nmElep0);
  nmAcc0_1->SetTitle("numu acceptance (fiducial && recontruction && reco_nu / true_nu)");
  nmAcc0_1->GetXaxis()->SetTitle("LepE (GeV)");
  //nmAcc0_1->GetYaxis()->SetTitle("reco_cut/true_cut");
  nmAcc3_1->Divide(nmElep3);
  neAcc0_1->Divide(neElep0);
  neAcc0_1->SetTitle("nue acceptance (fiducial && recontruction && reco_nu / true_nu)");
  neAcc0_1->GetXaxis()->SetTitle("LepE (GeV)");
  //neAcc0_1->GetYaxis()->SetTitle("reco_cut/true_cut");
  neAcc3_1->Divide(neElep3);

  TCanvas *c = new TCanvas("c","",1800,1200);
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLogy();
  nmAcc0->SetLineColor(kBlue);
  nmAcc0->Draw();
  nmAcc3->SetLineColor(kRed);
  nmAcc3->Draw("same");
  c->cd(2);
  gPad->SetLogy();
  neAcc0->SetLineColor(kBlue);
  neAcc0->Draw();
  neAcc3->SetLineColor(kRed);
  neAcc3->Draw("same");
  c->cd(3);
  gPad->SetLogy();
  nmAcc0_1->SetLineColor(kBlue);
  nmAcc0_1->Draw();
  nmAcc3_1->SetLineColor(kRed);
  nmAcc3_1->Draw("same");
  c->cd(4);
  gPad->SetLogy();
  neAcc0_1->SetLineColor(kBlue);
  neAcc0_1->Draw();
  neAcc3_1->SetLineColor(kRed);
  neAcc3_1->Draw("same");
  c->SaveAs("acc_LepE.png");
/*
  TF1 *PmueE = new TF1("PmueE","0.0004 * pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PmumuE = new TF1("PmumuE","1-2*sqrt(0.0004)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PeeE = new TF1("PeeE","1-2*sqrt(0.0004)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PmueE1 = new TF1("PmueE","0.01 * pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PmumuE1 = new TF1("PmumuE","1-2*sqrt(0.01)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PeeE1 = new TF1("PeeE","1-2*sqrt(0.01)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TCanvas *cpVsEv = new TCanvas("cpmueVsEv0","",1200,600);
  cpVsEv->Divide(3,2);
  cpVsEv->cd(1);
  gPad->SetLogx();
  gPad->SetLogz();
  pmueVsEv0->SetTitle("Pmue vs Ev (tmue = 0.01)");
  pmueVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  pmueVsEv0->GetYaxis()->SetTitle("P_{#mu e}");
  pmueVsEv0->Draw("colz");
  PmueE->Draw("same");
  cpVsEv->cd(2);
  gPad->SetLogx();
  gPad->SetLogz();
  pmumuVsEv0->SetTitle("Pmumu vs Ev (tmue = 0.01)");
  pmumuVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  pmumuVsEv0->GetYaxis()->SetTitle("P_{#mu#mu}");
  pmumuVsEv0->Draw("colz");
  PmumuE->Draw("same");
  cpVsEv->cd(3);
  gPad->SetLogx();
  gPad->SetLogz();
  peeVsEv0->SetTitle("Pee vs Ev (tmue = 0.01)");
  peeVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  peeVsEv0->GetYaxis()->SetTitle("P_{ee}");
  peeVsEv0->Draw("colz");
  PeeE->Draw("same");
  cpVsEv->cd(4);
  gPad->SetLogx();
  gPad->SetLogz();
  pmueVsEv1->SetTitle("Pmue vs Ev (tmue = 0.16)");
  pmueVsEv1->GetXaxis()->SetTitle("Ev (GeV)");
  pmueVsEv1->GetYaxis()->SetTitle("P_{#mu e}");
  pmueVsEv1->Draw("colz");
  PmueE1->Draw("same");
  cpVsEv->cd(5);
  gPad->SetLogx();
  gPad->SetLogz();
  pmumuVsEv1->SetTitle("Pmumu vs Ev (tmue = 0.16)");
  pmumuVsEv1->GetXaxis()->SetTitle("Ev (GeV)");
  pmumuVsEv1->GetYaxis()->SetTitle("P_{#mu#mu}");
  pmumuVsEv1->Draw("colz");
  PmumuE1->Draw("same");
  cpVsEv->cd(6);
  gPad->SetLogx();
  gPad->SetLogz();
  peeVsEv1->SetTitle("Pee vs Ev (tmue = 0.16)");
  peeVsEv1->GetXaxis()->SetTitle("Ev (GeV)");
  peeVsEv1->GetYaxis()->SetTitle("P_{ee}");
  peeVsEv1->Draw("colz");
  PeeE1->Draw("same");
  cpVsEv->SaveAs("pVsEv.pdf");
*/ 
}

