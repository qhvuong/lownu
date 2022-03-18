#include <string>

void lownu()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  
  for(int i = 100; i<150; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  std::cout << "File number:" << i << "\n";
  } 
  // replace this with the location of your file
  // Note: this step is going to print a bunch of Error failure loading library. These are not actually errors, just ROOT being stupid


  // Get the number of "protons on target", which is essentially the number of neutrinos in this file
/*  TChain * meta = new TChain( "meta", "meta" );
  meta->Add( "/Users/qhuong/Documents/UofR/2021_Summer/data/ND_FHC_FV_00.root" ); // make certain this is the exact same file(s)
*/
  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,14);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,14);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,14);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,14);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,14);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,14);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,14);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,14);
  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,14);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,14);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,14);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,14);

  const Int_t nbinsX = 150; const Int_t nbinsY = 140;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<80)                xEdges[i+1] = xEdges[i] + 0.05;
    else if(i>=80 && i<120) xEdges[i+1] = xEdges[i] + 0.1;
    else                    xEdges[i+1] = xEdges[i] + 0.2;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.1;
  }
  TH2D *m_hElepVsEv0 = new TH2D("m_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv3 = new TH2D("m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0 = new TH2D("e_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv3 = new TH2D("e_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
/*
  TH2D *pmueVsEv0 = new TH2D("pmueVsEv0","",100,0,30,100,0,0.0008);
  TH2D *pmueVsEv3 = new TH2D("pmueVsEv3","",100,0,30,100,0,0.0008);
  TH2D *pmumuVsEv0 = new TH2D("pmumuVsEv0","",100,0,30,100,0.95,1);
  TH2D *pmumuVsEv3 = new TH2D("pmumuVsEv3","",100,0,30,100,0.95,1);
  TH2D *peeVsEv0 = new TH2D("peeVsEv0","",100,0,30,100,0.95,1);
  TH2D *peeVsEv3 = new TH2D("peeVsEv3","",100,0,30,100,0.95,1);
*/

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

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

    double nu = eP + eN + ePip + ePim + ePi0 + eOther;
    double reco_nu = Ev_reco - Elep_reco;
      
    double Uee2 = 0.01;
    double Umm2 = 0.0016;
    double s2ee2 = 4*Uee2*(1 - Uee2); 
    double s2mm2 = 4*Umm2*(1 - Umm2); 
    double s2me2 = 4*Uee2*Umm2;
    double dm2 = 10; //eV2
    double L = 500; //m
    double del = 1.27/1E3*dm2*L/Ev;
    double pmue = s2me2 * pow(sin(del),2);
    double pmumu = 1 - s2mm2 * pow(sin(del),2);   
    double pee = 1 - s2ee2 * pow(sin(del),2);
    TRandom *r1 = new TRandom();
    double LepE_sm = r1->Gaus(LepE,0.05);

    if(LepPDG == 13){
    // Choose events that are reconstructed muons, where the muon is measured (i.e. doesn't exit detector), and where there isn't energy by the edges
    if( reco_numu )
    {
      //muon to muon ("survival")
      mm_hElep_w0->Fill(Elep_reco,pmumu);
      if(nu<0.3) mm_hElep_w3->Fill(Elep_reco,pmumu);

      //muon to electron
      me_hElep_w0->Fill(Elep_reco,pmue);
      if(nu<0.3) me_hElep_w3->Fill(LepE_sm,pmue);
      
      m_hElepVsEv0->Fill(Ev,Elep_reco);
      if(nu<0.3) m_hElepVsEv3->Fill(Ev,Elep_reco);
    }  
    }

    if(LepPDG == 11){
    
    if( reco_nue ){
      
      //electron to electron ("survival")
      ee_hElep_w0->Fill(LepE_sm,pee);
      if(nu<0.3) ee_hElep_w3->Fill(LepE_sm,pee);
      
      //electron to muon 
      em_hElep_w0->Fill(LepE_sm,pmue);
      if(nu<0.3) em_hElep_w3->Fill(LepE_sm,pmue);
      
      e_hElepVsEv0->Fill(Ev,LepE_sm);
      if(nu<0.3) e_hElepVsEv3->Fill(Ev,LepE_sm);
    }
    }
  }

  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  THStack *e_hsElep_w0 = new THStack("e_hsElep_w0","oscillated + unoscillated nue weighted Elep_reco");
  me_hElep_w0->SetMarkerStyle(21);
  me_hElep_w0->SetMarkerSize(0.5);
  me_hElep_w0->SetMarkerColor(kRed);
  me_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w0->Add(me_hElep_w0);
  ee_hElep_w0->SetMarkerStyle(21);
  ee_hElep_w0->SetMarkerSize(0.5);
  ee_hElep_w0->SetMarkerColor(kBlue);
  ee_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w0->Add(ee_hElep_w0);
  THStack *e_hsElep_w3 = new THStack("e_hsElep_w3","oscillated + unoscillated nue weighted Elep_reco (nu<0.3)");
  me_hElep_w3->SetMarkerStyle(21);
  me_hElep_w3->SetMarkerSize(0.5);
  me_hElep_w3->SetMarkerColor(kRed);
  me_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w3->Add(me_hElep_w3);
  ee_hElep_w3->SetMarkerStyle(21);
  ee_hElep_w3->SetMarkerSize(0.5);
  ee_hElep_w3->SetMarkerColor(kBlue);
  ee_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  e_hsElep_w3->Add(ee_hElep_w3);
  TCanvas *e_cElep_w = new TCanvas("e_cElep_w","",1600,600);
  e_cElep_w->Divide(2,1);
  e_cElep_w->cd(1);
  e_hsElep_w0->Draw();
  e_cElep_w->cd(2);
  e_hsElep_w3->Draw();
  e_cElep_w->SaveAs("e_Elep_reco_w_6.pdf");

  THStack *m_hsElep_w0 = new THStack("m_hsElep_w0","oscillated + unoscillated numu weighted Elep_reco");
  em_hElep_w0->SetMarkerStyle(21);
  em_hElep_w0->SetMarkerSize(0.5);
  em_hElep_w0->SetMarkerColor(kRed);
  em_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w0->Add(em_hElep_w0);
  mm_hElep_w0->SetMarkerStyle(21);
  mm_hElep_w0->SetMarkerSize(0.5);
  mm_hElep_w0->SetMarkerColor(kBlue);
  mm_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w0->Add(mm_hElep_w0);
  THStack *m_hsElep_w3 = new THStack("m_hsElep_w3","oscillated + unoscillated numu weighted Elep_reco (nu<0.3)");
  em_hElep_w3->SetMarkerStyle(21);
  em_hElep_w3->SetMarkerSize(0.5);
  em_hElep_w3->SetMarkerColor(kRed);
  em_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w3->Add(em_hElep_w3);
  mm_hElep_w3->SetMarkerStyle(21);
  mm_hElep_w3->SetMarkerSize(0.5);
  mm_hElep_w3->SetMarkerColor(kBlue);
  mm_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hsElep_w3->Add(mm_hElep_w3);
  TCanvas *m_cElep_w = new TCanvas("m_cElep_w","",1600,600);
  m_cElep_w->Divide(2,1);
  m_cElep_w->cd(1);
  m_hsElep_w0->Draw();
  m_cElep_w->cd(2);
  m_hsElep_w3->Draw();
  m_cElep_w->SaveAs("m_Elep_reco_w_6.pdf");
  
  e_hElep_w0->SetTitle("nue weighted Elep_reco");
  e_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hElep_w3->SetTitle("nue weighted Elep_reco (nu<0.3)");
  e_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w0->SetTitle("numu weighted Elep_reco");
  m_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w3->SetTitle("numu weighted Elep_reco (nu<0.3)");
  m_hElep_w3->GetXaxis()->SetTitle("Elep_reco");

  m_hElepVsEv0->SetTitle("numu Elep_reco Vs true Ev");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  m_hElepVsEv3->SetTitle("numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv0->SetTitle("nue Elep_reco Vs true Ev");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv3->SetTitle("nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",800,600);
  cElepVsEv->Divide(2,2);
  cElepVsEv->cd(1);
  m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  m_hElepVsEv3->Draw("colz");
  cElepVsEv->cd(3);
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(4);
  e_hElepVsEv3->Draw("colz");
  cElepVsEv->SaveAs("ElepVsEv_6.pdf");
/*
  TF1 *PmueE = new TF1("PmueE","0.0003 * pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PmumuE = new TF1("PmumuE","1-2*sqrt(0.0003)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TF1 *PeeE = new TF1("PeeE","1-2*sqrt(0.0003)*pow(sin(1.27/1E3*4*500/x),2)",0,30);
  TCanvas *cpVsEv = new TCanvas("cpmueVsEv0","",900,900);
  cpVsEv->SetLogx();
  cpVsEv->SetLogz();
  cpVsEv->Divide(2,3);
  cpVsEv->cd(1);
  pmueVsEv0->SetTitle("Pmue vs Ev");
  pmueVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  pmueVsEv0->GetYaxis()->SetTitle("P_{#mu e}");
  pmueVsEv0->Draw("colz");
  PmueE->Draw("same");
  cpVsEv->cd(2);
  pmueVsEv3->SetTitle("Pmue vs Ev (nu<0.3)");
  pmueVsEv3->GetXaxis()->SetTitle("Ev (GeV)");
  pmueVsEv3->GetYaxis()->SetTitle("P_{#mu e}");
  pmueVsEv3->Draw("colz");
  PmueE->Draw("same");
  cpVsEv->cd(3);
  pmumuVsEv0->SetTitle("Pmumu vs Ev");
  pmumuVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  pmumuVsEv0->GetYaxis()->SetTitle("P_{#mu#mu}");
  pmumuVsEv0->Draw("colz");
  PmumuE->Draw("same");
  cpVsEv->cd(4);
  pmumuVsEv3->SetTitle("Pmumu vs Ev (nu<0.3)");
  pmumuVsEv3->GetXaxis()->SetTitle("Ev (GeV)");
  pmumuVsEv3->GetYaxis()->SetTitle("P_{#mu#mu}");
  pmumuVsEv3->Draw("colz");
  PmumuE->Draw("same");
  cpVsEv->cd(5);
  peeVsEv0->SetTitle("Pee vs Ev");
  peeVsEv0->GetXaxis()->SetTitle("Ev (GeV)");
  peeVsEv0->GetYaxis()->SetTitle("P_{ee}");
  peeVsEv0->Draw("colz");
  PeeE->Draw("same");
  cpVsEv->cd(6);
  peeVsEv3->SetTitle("Pee vs Ev (nu<0.3)");
  peeVsEv3->GetXaxis()->SetTitle("Ev (GeV)");
  peeVsEv3->GetYaxis()->SetTitle("P_{ee}");
  peeVsEv3->Draw("colz");
  PeeE->Draw("same");
  cpVsEv->SaveAs("pVsEv.pdf");
*/ 

  TFile *out = new TFile("/nashome/q/qvuong/gen_data/output_6.root","RECREATE");
  e_hElep_w0->Write();
  e_hElep_w3->Write();
  me_hElep_w0->Write();
  me_hElep_w3->Write();
  ee_hElep_w0->Write();
  ee_hElep_w3->Write();
  m_hElep_w0->Write();
  m_hElep_w3->Write();
  em_hElep_w0->Write();
  em_hElep_w3->Write();
  mm_hElep_w0->Write();
  mm_hElep_w3->Write();
  m_hElepVsEv0->Write();
  m_hElepVsEv3->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv3->Write();
  out->Close();

}

