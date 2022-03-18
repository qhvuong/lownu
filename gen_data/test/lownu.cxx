#include <iostream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TChain.h>

main()
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
  //tree->Add("../data/ND_FHC_FV_00.root");

  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,40);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,16);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,40);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,16);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,40);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,16);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,40);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,40);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,40);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0_ft = new TH1D("e_hElep_w0_ft","",100,0,40);
  TH1D *e_hElep_w3_ft = new TH1D("e_hElep_w3_ft","",100,0,16);
  TH1D *m_hElep_w0_ft = new TH1D("m_hElep_w0_ft","",100,0,40);
  TH1D *m_hElep_w3_ft = new TH1D("m_hElep_w3_ft","",100,0,16);

  const Int_t nbinsX = 280; const Int_t nbinsY = 100;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else                     xEdges[i+1] = xEdges[i] + 0.2;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.16;
  }
  
  const Int_t nbinsX0 = 400; const Int_t nbinsY0 = 100;
  Double_t xEdges0[nbinsX0+1], yEdges0[nbinsY0+1];
  xEdges0[0]=yEdges0[0]=0;
  for(int i=0; i<nbinsX0+1; i++)
  {
    if(i<200)                xEdges0[i+1] = xEdges0[i] + 0.02;
    else if(i>=200 && i<240) xEdges0[i+1] = xEdges0[i] + 0.1;
    //else if(i>=240 && i<280) xEdges0[i+1] = xEdges0[i] + 0.2;
    else                     xEdges0[i+1] = xEdges0[i] + 0.2;
  }
  for(int i=0; i<nbinsY0+1; i++)
  {
    yEdges0[i+1]=yEdges0[i]+0.4;
  }
   
  TH2D *m_hElepVsEv0 = new TH2D("m_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *m_hElepVsEv3 = new TH2D("m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepVsEv0 = new TH2D("nc_m_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *nc_m_hElepVsEv3 = new TH2D("nc_m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0 = new TH2D("e_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *e_hElepVsEv3 = new TH2D("e_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);

  TH2D *e_ThetaVsElep0 = new TH2D("e_ThetaVsElep0","",400,0.5,60,20,0,1200);
  TH2D *e_ThetaVsElep3 = new TH2D("e_ThetaVsElep3","",400,0.5,60,20,0,1200);
  TH2D *m_ThetaVsElep0 = new TH2D("m_ThetaVsElep0","",400,0.5,60,20,0,1200);
  TH2D *m_ThetaVsElep3 = new TH2D("m_ThetaVsElep3","",400,0.5,60,20,0,1200);

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

  double nu, reco_nu;
  double del, pmue, pmumu, pee;
  double LepE_sm;

  double Uee2 = 0.04;
  double Umm2 = 0.01;
  double s2ee2 = 4*Uee2*(1 - Uee2); 
  double s2mm2 = 4*Umm2*(1 - Umm2); 
  double s2me2 = 4*Uee2*Umm2;
  double dm2 = 6.0; //eV2
  double L = 500; //m

  TRandom *r1 = new TRandom(8888);
  TRandom3 *rando = new TRandom3(8888); 
 
  const int N = tree->GetEntries();
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
    //if( Ev < 0.2 ) continue;

    nu = eP + eN + ePip + ePim + ePi0 + eOther;
    reco_nu = Ev_reco - Elep_reco;
      
    del = 1.27/1E3*dm2*L/Ev;
    pmue = s2me2 * pow(sin(del),2);
    pmumu = 1 - s2mm2 * pow(sin(del),2);   
    pee = 1 - s2ee2 * pow(sin(del),2);
    LepE_sm = r1->Gaus(LepE,0.05);

    if(LepPDG == 13){

      if(nu<10.0) m_ThetaVsElep0->Fill(LepE_sm,LepNuAngle*1E3);
      if(nu<0.3)  m_ThetaVsElep3->Fill(LepE_sm,LepNuAngle*1E3);

      //target: muon to electron 
      if(nu<10.0) me_hElep_w0->Fill(LepE_sm,pmue);
      if(nu<0.3)  me_hElep_w3->Fill(LepE_sm,pmue);
        
      //templates: muon without cuts
      if(nu<10.0) nc_m_hElepVsEv0->Fill(Ev,LepE_sm);    
      if(nu<0.3)  nc_m_hElepVsEv3->Fill(Ev,LepE_sm);
    
      if( reco_numu && (muon_contained || muon_tracker || muon_ecal))
      {
        //target: muon to muon ("survival")
        if(nu<10.0) mm_hElep_w0->Fill(Elep_reco,pmumu);    
        if(nu<0.3)  mm_hElep_w3->Fill(Elep_reco,pmumu);
      
        //templates: muon with cuts
        if(nu<10.0) m_hElepVsEv0->Fill(Ev,Elep_reco);    
        if(nu<0.3)  m_hElepVsEv3->Fill(Ev,Elep_reco);
      }  
    }

    if(LepPDG == 11){
      
      if(nu<10.0) e_ThetaVsElep0->Fill(LepE_sm,LepNuAngle*1E3);
      if(nu<0.3)  e_ThetaVsElep3->Fill(LepE_sm,LepNuAngle*1E3);
    
      if(LepE_sm*LepNuAngle*LepNuAngle*1E3>3.){
      if( reco_nue ){
      
      //target: electron to electron ("survival")
      if(nu<10.0) ee_hElep_w0->Fill(LepE_sm,pee);
      if(nu<0.3)  ee_hElep_w3->Fill(LepE_sm,pee);
      
      //target: electron to muon 
      if(nu<10.0) em_hElep_w0->Fill(LepE_sm,pmue);
      if(nu<0.3)  em_hElep_w3->Fill(LepE_sm,pmue);
      
      //templates: eletron
      if(nu<10.0) e_hElepVsEv0->Fill(Ev,LepE_sm);
      if(nu<0.3)  e_hElepVsEv3->Fill(Ev,LepE_sm);
      }
      }
    }
  }

  
  //target: electron = muon-to-electron + electron-to-electron
  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);

  //target: muon = electron-to-muon + muon-to-muon
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);
 
 
  for( int bx = 1; bx <= m_hElep_w0->GetNbinsX(); bx++ ){
    double mean = m_hElep_w0->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    m_hElep_w0_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= m_hElep_w3->GetNbinsX(); bx++ ){
    double mean = m_hElep_w3->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    m_hElep_w3_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= e_hElep_w0->GetNbinsX(); bx++ ){
    double mean = e_hElep_w0->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    e_hElep_w0_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= e_hElep_w3->GetNbinsX(); bx++ ){
    double mean = e_hElep_w3->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    e_hElep_w3_ft->AddBinContent(bx, fluctuated_bin_content);
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  e_hElep_w0_ft->SetStats(0);
  e_hElep_w0_ft->SetTitle("fluctuated nueCC target (nu<10.0)");
  e_hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep_w3_ft->SetStats(0);
  e_hElep_w3_ft->SetTitle("fluctuated nueCC target (nu<0.3)");
  e_hElep_w3_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep_w0_ft->SetStats(0);
  m_hElep_w0_ft->SetTitle("fluctuated numuCC target (nu<10.0)");
  m_hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep_w3_ft->SetStats(0);
  m_hElep_w3_ft->SetTitle("fluctuated numuCC target (nu<0.3)");
  m_hElep_w3_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cTarget = new TCanvas("cTarget","",1600,1200);
  cTarget->Divide(2,2);
  cTarget->cd(1);
  e_hElep_w0_ft->Draw();
  cTarget->cd(2);
  e_hElep_w3_ft->Draw();
  cTarget->cd(3);
  m_hElep_w0_ft->Draw();
  cTarget->cd(4);
  m_hElep_w3_ft->Draw();
  cTarget->SaveAs("CC_target_ft_2.pdf");

  THStack *e_hsElep_w0 = new THStack("e_hsElep_w0","oscillated + unoscillated nue weighted Elep_reco (nu<10.0)");
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
  e_cElep_w->SaveAs("CC_e_Elep_reco_w_2.pdf");

  THStack *m_hsElep_w0 = new THStack("m_hsElep_w0","oscillated + unoscillated numu weighted Elep_reco (nu<10.0)");
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
  m_cElep_w->SaveAs("CC_m_Elep_reco_w_2.pdf");
  
  e_hElep_w0->SetTitle("nue weighted Elep_reco (nu<10.0)");
  e_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hElep_w3->SetTitle("nue weighted Elep_reco (nu<0.3)");
  e_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w0->SetTitle("numu weighted Elep_reco (nu<10.0)");
  m_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w3->SetTitle("numu weighted Elep_reco (nu<0.3)");
  m_hElep_w3->GetXaxis()->SetTitle("Elep_reco");

  m_hElepVsEv0->SetTitle("numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  m_hElepVsEv3->SetTitle("numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  nc_m_hElepVsEv0->SetTitle("numu Elep_reco Vs true Ev (nu<10.0) without reconstruction cuts");
  nc_m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  nc_m_hElepVsEv3->SetTitle("numu Elep_reco Vs true Ev (nu<0.3) without reconstruction cuts");
  nc_m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv0->SetTitle("nue Elep_reco Vs true Ev (nu<10.0)");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv3->SetTitle("nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",1200,600);
  cElepVsEv->Divide(3,2);
  cElepVsEv->cd(1);
  m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  nc_m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(3);
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(4);
  m_hElepVsEv3->Draw("colz");
  cElepVsEv->cd(5);
  nc_m_hElepVsEv3->Draw("colz");
  cElepVsEv->cd(6);
  e_hElepVsEv3->Draw("colz");
  cElepVsEv->SaveAs("CC_ElepVsEv_2.pdf");
/*
  TF1 *cut = new TF1("cut","sqrt(3*1E3/x)",0.5,60);
  m_ThetaVsElep0->SetTitle("numuCC Theta Vs Elep (nu<10.0)");
  m_ThetaVsElep0->GetXaxis()->SetTitle("Elep (GeV)");
  m_ThetaVsElep0->GetYaxis()->SetTitle("LepNuAngle (mrad)");
  m_ThetaVsElep3->SetTitle("numuCC Theta Vs Elep (nu<0.3)");
  m_ThetaVsElep3->GetXaxis()->SetTitle("Elep (GeV)");
  m_ThetaVsElep3->GetYaxis()->SetTitle("LepNuAngle (mrad)");
  e_ThetaVsElep0->SetTitle("nueCC Theta Vs Elep (nu<10.0)");
  e_ThetaVsElep0->GetXaxis()->SetTitle("Elep (GeV)");
  e_ThetaVsElep0->GetYaxis()->SetTitle("LepNuAngle (mrad)");
  e_ThetaVsElep3->SetTitle("nueCC Theta Vs Elep (nu<0.3)");
  e_ThetaVsElep3->GetXaxis()->SetTitle("Elep (GeV)");
  e_ThetaVsElep3->GetYaxis()->SetTitle("LepNuAngle (mrad)");
  TCanvas *cThetaVsElep = new TCanvas("cThetaVsElep","",800,600);
  cThetaVsElep->Divide(2,2);
  cThetaVsElep->cd(1);
  gPad->SetLogx();
  //gPad->SetLogy();
  m_ThetaVsElep0->Draw("colz");
  cut->Draw("same");
  cThetaVsElep->cd(2);
  gPad->SetLogx();
  e_ThetaVsElep0->Draw("colz");
  cut->Draw("same");
  cThetaVsElep->cd(3);
  gPad->SetLogx();
  //gPad->SetLogy();
  m_ThetaVsElep3->Draw("colz");
  cut->Draw("same");
  cThetaVsElep->cd(4);
  gPad->SetLogx();
  //gPad->SetLogy();
  e_ThetaVsElep3->Draw("colz");
  cut->Draw("same");
  cThetaVsElep->SaveAs("CC_ThetaVsElep.pdf");
*/
  TFile *out = new TFile("/nashome/q/qvuong/gen_data/output_2.root","RECREATE");
  e_hElep_w0->Write();
  e_hElep_w3->Write();
  e_hElep_w0_ft->Write();
  e_hElep_w3_ft->Write();
  me_hElep_w0->Write();
  me_hElep_w3->Write();
  ee_hElep_w0->Write();
  ee_hElep_w3->Write();
  m_hElep_w0->Write();
  m_hElep_w3->Write();
  m_hElep_w0_ft->Write();
  m_hElep_w3_ft->Write();
  em_hElep_w0->Write();
  em_hElep_w3->Write();
  mm_hElep_w0->Write();
  mm_hElep_w3->Write();
  m_hElepVsEv0->Write();
  m_hElepVsEv3->Write();
  nc_m_hElepVsEv0->Write();
  nc_m_hElepVsEv3->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv3->Write();
  out->Close();

}

