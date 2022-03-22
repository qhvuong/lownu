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
#include <TLegend.h>

int main()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  TChain * meta = new TChain( "meta", "meta" );
  
  for(int i = 200; i<250; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  meta->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) ); // make certain this is the exact same file(s)
  std::cout << "File number:" << i << "\n";
  } 

  double total_pot = 0.;
  double pot;
  meta->SetBranchAddress( "pot", &pot );
  const int Nfiles = meta->GetEntries();
  for( int ii = 0; ii < Nfiles; ++ii ) {
    meta->GetEntry(ii);
    total_pot += pot;
  }
  double yrPOT = 1.1E21;
  double scalePOT = yrPOT/total_pot;

  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,16);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,16);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,16);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,16);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,16);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,16);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,16);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0  = new TH1D("e_hElep_w0","",100,0,16);
  TH1D *e_hElep_w3  = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0  = new TH1D("m_hElep_w0","",100,0,16);
  TH1D *m_hElep_w3  = new TH1D("m_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0_ft = new TH1D("e_hElep_w0_ft","",100,0,16);
  TH1D *e_hElep_w3_ft = new TH1D("e_hElep_w3_ft","",100,0,16);
  TH1D *m_hElep_w0_ft = new TH1D("m_hElep_w0_ft","",100,0,16);
  TH1D *m_hElep_w3_ft = new TH1D("m_hElep_w3_ft","",100,0,16);

  TH1D *mm_hEvReco_w0 = new TH1D("mm_hEvReco_w0","",100,0,20);
  TH1D *mm_hEvReco_w3 = new TH1D("mm_hEvReco_w3","",100,0,20);
  TH1D *ee_hEvReco_w0 = new TH1D("ee_hEvReco_w0","",100,0,20);
  TH1D *ee_hEvReco_w3 = new TH1D("ee_hEvReco_w3","",100,0,20);
  TH1D *me_hEvReco_w0 = new TH1D("me_hEvReco_w0","",100,0,20);
  TH1D *me_hEvReco_w3 = new TH1D("me_hEvReco_w3","",100,0,20);
  TH1D *em_hEvReco_w0 = new TH1D("em_hEvReco_w0","",100,0,20);
  TH1D *em_hEvReco_w3 = new TH1D("em_hEvReco_w3","",100,0,20);
  TH1D *e_hEvReco_w0  = new TH1D("e_hEvReco_w0","",100,0,20);
  TH1D *e_hEvReco_w3  = new TH1D("e_hEvReco_w3","",100,0,20);
  TH1D *m_hEvReco_w0  = new TH1D("m_hEvReco_w0","",100,0,20);
  TH1D *m_hEvReco_w3  = new TH1D("m_hEvReco_w3","",100,0,20);
  TH1D *e_hEvReco_w0_ft = new TH1D("e_hEvReco_w0_ft","",100,0,20);
  TH1D *e_hEvReco_w3_ft = new TH1D("e_hEvReco_w3_ft","",100,0,20);
  TH1D *m_hEvReco_w0_ft = new TH1D("m_hEvReco_w0_ft","",100,0,20);
  TH1D *m_hEvReco_w3_ft = new TH1D("m_hEvReco_w3_ft","",100,0,20);

  const Int_t nbinsX = 480; const Int_t nbinsY = 100;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else                     xEdges[i+1] = xEdges[i] + 1.0;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.16;
  }

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

  TH2D *m_hElepVsEv0    = new TH2D("m_hElepVsEv0","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv3    = new TH2D("m_hElepVsEv3","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepVsEv0 = new TH2D("nc_m_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepVsEv3 = new TH2D("nc_m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0    = new TH2D("e_hElepVsEv0","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv3    = new TH2D("e_hElepVsEv3","",   nbinsX,xEdges,nbinsY,yEdges);

  TH2D *m_hElepVsEv0_cov = new TH2D("m_hElepVsEv0_cov","",19,mubins,100,0,16);
  TH2D *m_hElepVsEv3_cov = new TH2D("m_hElepVsEv3_cov","",19,mubins,100,0,16);
  TH2D *e_hElepVsEv0_cov = new TH2D("e_hElepVsEv0_cov","",7, ebins,100,0,16);
  TH2D *e_hElepVsEv3_cov = new TH2D("e_hElepVsEv3_cov","",7, ebins,100,0,16);

  TH2D *m_hEvRecoVsEv0    = new TH2D("m_hEvRecoVsEv0","",   nbinsX,xEdges,100,0,20);
  TH2D *m_hEvRecoVsEv3    = new TH2D("m_hEvRecoVsEv3","",   nbinsX,xEdges,100,0,20);
  TH2D *nc_m_hEvRecoVsEv0 = new TH2D("nc_m_hEvRecoVsEv0","",nbinsX,xEdges,100,0,20);
  TH2D *nc_m_hEvRecoVsEv3 = new TH2D("nc_m_hEvRecoVsEv3","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv0    = new TH2D("e_hEvRecoVsEv0","",   nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv1    = new TH2D("e_hEvRecoVsEv1","",   nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv2    = new TH2D("e_hEvRecoVsEv2","",   nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv3    = new TH2D("e_hEvRecoVsEv3","",   nbinsX,xEdges,100,0,20);

  TH2D *m_hEvRecoVsEv0_cov = new TH2D("m_hEvRecoVsEv0_cov","",19,mubins,100,0,20);
  TH2D *m_hEvRecoVsEv3_cov = new TH2D("m_hEvRecoVsEv3_cov","",19,mubins,100,0,20);
  TH2D *e_hEvRecoVsEv0_cov = new TH2D("e_hEvRecoVsEv0_cov","",7, ebins,100,0,20);
  TH2D *e_hEvRecoVsEv1_cov = new TH2D("e_hEvRecoVsEv1_cov","",7, ebins,100,0,20);
  TH2D *e_hEvRecoVsEv2_cov = new TH2D("e_hEvRecoVsEv2_cov","",7, ebins,100,0,20);
  TH2D *e_hEvRecoVsEv3_cov = new TH2D("e_hEvRecoVsEv3_cov","",7, ebins,100,0,20);

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
  //double Ev_reco, Elep_reco; // the measured neutrino energy and lepton energy. Reco nu = Ev_reco - Elep_reco
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
  //tree->SetBranchAddress( "Ev_reco", &Ev_reco );
  //tree->SetBranchAddress( "Elep_reco", &Elep_reco );
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
  double LepE_sm, Ev_sm;

  double Uee2 = 0.01;
  double Umm2 = 0.0016;
  double s2ee2 = 4*Uee2*(1 - Uee2); 
  double s2mm2 = 4*Umm2*(1 - Umm2); 
  double s2me2 = 4*Uee2*Umm2;
  double dm2 = 1.3; //eV2
  double L = 500; //m

  TRandom *r1 = new TRandom(8888);
  TRandom3 *rando = new TRandom3(8888); 
 
  const int N = tree->GetEntries();
  //const int N = 1000;
  //double scalePOT = 1.;
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
    //if( Ev < 0.2 ) continue;

    //nu = eP + eN + ePip + ePim + ePi0 + eOther;
      
    del = 1.27/1E3*dm2*L/Ev;
    pmue = s2me2 * pow(sin(del),2);
    pmumu = 1 - s2mm2 * pow(sin(del),2);   
    pee = 1 - s2ee2 * pow(sin(del),2);
    LepE_sm = r1->Gaus(LepE,0.05);
    Ev_sm   = r1->Gaus(Ev,0.05);
    nu = Ev_sm - LepE_sm; //reco_nu

    if(LepPDG == 13){

      if(nu<10.0) m_ThetaVsElep0->Fill(LepE_sm,LepNuAngle*1E3);
      if(nu<0.3)  m_ThetaVsElep3->Fill(LepE_sm,LepNuAngle*1E3);

      //target: muon to electron 
      if(nu<10.0) { me_hElep_w0->Fill(LepE_sm,pmue); me_hEvReco_w0->Fill(Ev_sm,pmue); }
      if(nu<0.3)  { me_hElep_w3->Fill(LepE_sm,pmue); me_hEvReco_w3->Fill(Ev_sm,pmue); }
        
      //templates: muon without cuts
      if(nu<10.0) { nc_m_hElepVsEv0->Fill(Ev,LepE_sm); nc_m_hEvRecoVsEv0->Fill(Ev,Ev_sm); }   
      if(nu<0.3)  { nc_m_hElepVsEv3->Fill(Ev,LepE_sm); nc_m_hEvRecoVsEv3->Fill(Ev,Ev_sm); }
    
      if( reco_numu && (muon_contained || muon_tracker || muon_ecal))
      {
        //target: muon to muon ("survival")
        if(nu<10.0) { mm_hElep_w0->Fill(LepE_sm,pmumu); mm_hEvReco_w0->Fill(Ev_sm,pmumu); }   
        if(nu<0.3)  { mm_hElep_w3->Fill(LepE_sm,pmumu); mm_hEvReco_w3->Fill(Ev_sm,pmumu); }
      
        //templates: muon with cuts
        if(nu<10.0) { m_hElepVsEv0->Fill(Ev,LepE_sm); m_hElepVsEv0_cov->Fill(Ev,LepE_sm);
		      m_hEvRecoVsEv0->Fill(Ev,Ev_sm); m_hEvRecoVsEv0_cov->Fill(Ev,Ev_sm); }    
        if(nu<0.3)  { m_hElepVsEv3->Fill(Ev,LepE_sm); m_hElepVsEv3_cov->Fill(Ev,LepE_sm); 
		      m_hEvRecoVsEv3->Fill(Ev,Ev_sm); m_hEvRecoVsEv3_cov->Fill(Ev,Ev_sm); }
      }  
    }

    if(LepPDG == 11){
      
      if(nu<10.0) e_ThetaVsElep0->Fill(LepE_sm,LepNuAngle*1E3);
      if(nu<0.3)  e_ThetaVsElep3->Fill(LepE_sm,LepNuAngle*1E3);
    
      if(LepE_sm*LepNuAngle*LepNuAngle*1E3>3.){
      if( reco_nue ){
      
      //target: electron to electron ("survival")
      if(nu<10.0) { ee_hElep_w0->Fill(LepE_sm,pee); ee_hEvReco_w0->Fill(Ev_sm,pee); }
      if(nu<0.3)  { ee_hElep_w3->Fill(LepE_sm,pee); ee_hEvReco_w3->Fill(Ev_sm,pee); }
      
      //target: electron to muon 
      if(nu<10.0) { em_hElep_w0->Fill(LepE_sm,pmue); em_hEvReco_w0->Fill(Ev_sm,pmue); }
      if(nu<0.3)  { em_hElep_w3->Fill(LepE_sm,pmue); em_hEvReco_w3->Fill(Ev_sm,pmue); }
      
      //templates: eletron
      if(nu<10.0) { e_hElepVsEv0->Fill(Ev,LepE_sm);   e_hElepVsEv0_cov->Fill(Ev,LepE_sm); 
		    e_hEvRecoVsEv0->Fill(Ev,Ev_sm); e_hEvRecoVsEv0_cov->Fill(Ev,Ev_sm); }
      if(nu<0.3)  { e_hElepVsEv3->Fill(Ev,LepE_sm);   e_hElepVsEv3_cov->Fill(Ev,LepE_sm); 
		    e_hEvRecoVsEv3->Fill(Ev,Ev_sm); e_hEvRecoVsEv3_cov->Fill(Ev,Ev_sm); }
      }
      }
    }
  }
  
  ee_hElep_w0->Scale(scalePOT); em_hElep_w0->Scale(scalePOT);
  ee_hElep_w3->Scale(scalePOT); em_hElep_w3->Scale(scalePOT);
  mm_hElep_w0->Scale(scalePOT); me_hElep_w0->Scale(scalePOT);
  mm_hElep_w3->Scale(scalePOT); me_hElep_w3->Scale(scalePOT);
  ee_hEvReco_w0->Scale(scalePOT); em_hEvReco_w0->Scale(scalePOT);
  ee_hEvReco_w3->Scale(scalePOT); em_hEvReco_w3->Scale(scalePOT);
  mm_hEvReco_w0->Scale(scalePOT); me_hEvReco_w0->Scale(scalePOT);
  mm_hEvReco_w3->Scale(scalePOT); me_hEvReco_w3->Scale(scalePOT);

  e_hElepVsEv0->Scale(scalePOT);    e_hElepVsEv3->Scale(scalePOT);
  m_hElepVsEv0->Scale(scalePOT);    m_hElepVsEv3->Scale(scalePOT);
  nc_m_hElepVsEv0->Scale(scalePOT); nc_m_hElepVsEv3->Scale(scalePOT);
  e_hEvRecoVsEv0->Scale(scalePOT);    e_hEvRecoVsEv3->Scale(scalePOT);
  m_hEvRecoVsEv0->Scale(scalePOT);    m_hEvRecoVsEv3->Scale(scalePOT);
  nc_m_hEvRecoVsEv0->Scale(scalePOT); nc_m_hEvRecoVsEv3->Scale(scalePOT);

  e_hElepVsEv0_cov->Scale(scalePOT); e_hElepVsEv3_cov->Scale(scalePOT);
  m_hElepVsEv0_cov->Scale(scalePOT); m_hElepVsEv3_cov->Scale(scalePOT);
  e_hEvRecoVsEv0_cov->Scale(scalePOT); e_hEvRecoVsEv3_cov->Scale(scalePOT);
  m_hEvRecoVsEv0_cov->Scale(scalePOT); m_hEvRecoVsEv3_cov->Scale(scalePOT);
  
  //target: electron = muon-to-electron + electron-to-electron
  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);
  e_hEvReco_w0->Add(me_hEvReco_w0); e_hEvReco_w0->Add(ee_hEvReco_w0);
  e_hEvReco_w3->Add(me_hEvReco_w3); e_hEvReco_w3->Add(ee_hEvReco_w3);

  //target: muon = electron-to-muon + muon-to-muon
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);
  m_hEvReco_w0->Add(em_hEvReco_w0); m_hEvReco_w0->Add(mm_hEvReco_w0);
  m_hEvReco_w3->Add(em_hEvReco_w3); m_hEvReco_w3->Add(mm_hEvReco_w3);
 
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

  for( int bx = 1; bx <= m_hEvReco_w0->GetNbinsX(); bx++ ){
    double mean = m_hEvReco_w0->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    m_hEvReco_w0_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= m_hEvReco_w3->GetNbinsX(); bx++ ){
    double mean = m_hEvReco_w3->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    m_hEvReco_w3_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= e_hEvReco_w0->GetNbinsX(); bx++ ){
    double mean = e_hEvReco_w0->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    e_hEvReco_w0_ft->AddBinContent(bx, fluctuated_bin_content);
  }
  for( int bx = 1; bx <= e_hEvReco_w3->GetNbinsX(); bx++ ){
    double mean = e_hEvReco_w3->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    e_hEvReco_w3_ft->AddBinContent(bx, fluctuated_bin_content);
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  e_hElep_w0_ft->SetStats(0);
  e_hElep_w0_ft->SetTitle("fluctuated nueCC Elep target (nu<10.0)");
  e_hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep_w0_ft->GetYaxis()->SetTitle("counts/yr");
  e_hElep_w0_ft->SetLineColor(kBlack);
  e_hElep_w3_ft->SetStats(0);
  e_hElep_w3_ft->SetTitle("fluctuated nueCC Elep target (nu<0.3)");
  e_hElep_w3_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep_w3_ft->GetYaxis()->SetTitle("counts/yr");
  e_hElep_w3_ft->SetLineColor(kBlack);
  m_hElep_w0_ft->SetStats(0);
  m_hElep_w0_ft->SetTitle("fluctuated numuCC Elep target (nu<10.0)");
  m_hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep_w0_ft->GetYaxis()->SetTitle("counts/yr");
  m_hElep_w0_ft->SetLineColor(kBlack);
  m_hElep_w3_ft->SetStats(0);
  m_hElep_w3_ft->SetTitle("fluctuated numuCC Elep target (nu<0.3)");
  m_hElep_w3_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep_w3_ft->GetYaxis()->SetTitle("counts/yr");
  m_hElep_w3_ft->SetLineColor(kBlack);

  me_hElep_w0->SetMarkerStyle(21);
  me_hElep_w0->SetMarkerSize(0.5);
  me_hElep_w0->SetMarkerColor(kBlue);
  ee_hElep_w0->SetMarkerStyle(21);
  ee_hElep_w0->SetMarkerSize(0.5);
  ee_hElep_w0->SetMarkerColor(kRed);
  em_hElep_w0->SetMarkerStyle(21);
  em_hElep_w0->SetMarkerSize(0.5);
  em_hElep_w0->SetMarkerColor(kBlue);
  mm_hElep_w0->SetMarkerStyle(21);
  mm_hElep_w0->SetMarkerSize(0.5);
  mm_hElep_w0->SetMarkerColor(kRed);

  me_hElep_w3->SetMarkerStyle(21);
  me_hElep_w3->SetMarkerSize(0.5);
  me_hElep_w3->SetMarkerColor(kBlue);
  ee_hElep_w3->SetMarkerStyle(21);
  ee_hElep_w3->SetMarkerSize(0.5);
  ee_hElep_w3->SetMarkerColor(kRed);
  em_hElep_w3->SetMarkerStyle(21);
  em_hElep_w3->SetMarkerSize(0.5);
  em_hElep_w3->SetMarkerColor(kRed);
  mm_hElep_w3->SetMarkerStyle(21);
  mm_hElep_w3->SetMarkerSize(0.5);
  mm_hElep_w3->SetMarkerColor(kRed);

  THStack *e_hsElep_w0 = new THStack("e_hsElep_w0","oscillated + unoscillated nue weighted Elep_reco (nu<10.0); Elep_reco (GeV); counts/yr");
  e_hsElep_w0->Add(me_hElep_w0);
  e_hsElep_w0->Add(ee_hElep_w0);
  THStack *e_hsElep_w3 = new THStack("e_hsElep_w3","oscillated + unoscillated nue weighted Elep_reco (nu<0.3); Elep_reco (GeV); counts/yr");
  e_hsElep_w3->Add(me_hElep_w3);
  e_hsElep_w3->Add(ee_hElep_w3);
  THStack *m_hsElep_w0 = new THStack("m_hsElep_w0","oscillated + unoscillated numu weighted Elep_reco (nu<10.0); Elep_reco (GeV); counts/yr");
  m_hsElep_w0->Add(em_hElep_w0);
  m_hsElep_w0->Add(mm_hElep_w0);
  THStack *m_hsElep_w3 = new THStack("m_hsElep_w3","oscillated + unoscillated numu weighted Elep_reco (nu<0.3); Elep_reco (GeV); counts/yr");
  m_hsElep_w3->Add(em_hElep_w3);
  m_hsElep_w3->Add(mm_hElep_w3);

  TCanvas *cElep_Target = new TCanvas("cElep_Target","",1800,1400);
  cElep_Target->Divide(2,2);
  cElep_Target->cd(1);
  e_hElep_w0_ft->Draw();
  e_hsElep_w0->Draw("same");
  TLegend *le0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  le0->AddEntry(e_hElep_w0_ft,"nueCC target");
  le0->AddEntry(me_hElep_w0,"oscillated mu->e");
  le0->AddEntry(ee_hElep_w0,"oscilated + unoscillated");
  le0->Draw();
  cElep_Target->cd(2);
  e_hElep_w3_ft->Draw();
  e_hsElep_w3->Draw("same");
  TLegend *le3 = new TLegend(0.65, 0.70, 0.9, 0.9);
  le3->AddEntry(e_hElep_w3_ft,"nueCC target");
  le3->AddEntry(me_hElep_w3,"oscillated mu->e");
  le3->AddEntry(ee_hElep_w3,"oscilated + unoscillated");
  le3->Draw();
  cElep_Target->cd(3);
  m_hElep_w0_ft->Draw();
  m_hsElep_w0->Draw("same");
  TLegend *lm0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lm0->AddEntry(m_hElep_w0_ft,"numuCC target");
  lm0->AddEntry(em_hElep_w0,"oscillated e->mu");
  lm0->AddEntry(mm_hElep_w0,"oscilated + unoscillated");
  lm0->Draw();
  cElep_Target->cd(4);
  m_hElep_w3_ft->Draw();
  m_hsElep_w3->Draw("same");
  TLegend *lm3 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lm3->AddEntry(m_hElep_w3_ft,"numuCC target");
  lm3->AddEntry(em_hElep_w3,"oscillated e->mu");
  lm3->AddEntry(mm_hElep_w3,"oscilated + unoscillated");
  lm3->Draw();
  cElep_Target->SaveAs("true_CC_Elep_target_1.png");

  e_hEvReco_w0_ft->SetStats(0);
  e_hEvReco_w0_ft->SetTitle("fluctuated nueCC Ev_reco target (nu<10.0)");
  e_hEvReco_w0_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco_w0_ft->GetYaxis()->SetTitle("counts/yr");
  e_hEvReco_w0_ft->SetLineColor(kBlack);
  e_hEvReco_w3_ft->SetStats(0);
  e_hEvReco_w3_ft->SetTitle("fluctuated nueCC Ev_reco target (nu<0.3)");
  e_hEvReco_w3_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco_w3_ft->GetYaxis()->SetTitle("counts/yr");
  e_hEvReco_w3_ft->SetLineColor(kBlack);
  m_hEvReco_w0_ft->SetStats(0);
  m_hEvReco_w0_ft->SetTitle("fluctuated numuCC Ev_reco target (nu<10.0)");
  m_hEvReco_w0_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco_w0_ft->GetYaxis()->SetTitle("counts/yr");
  m_hEvReco_w0_ft->SetLineColor(kBlack);
  m_hEvReco_w3_ft->SetStats(0);
  m_hEvReco_w3_ft->SetTitle("fluctuated numuCC Ev_reco target (nu<0.3)");
  m_hEvReco_w3_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco_w3_ft->GetYaxis()->SetTitle("counts/yr");
  m_hEvReco_w3_ft->SetLineColor(kBlack);

  me_hEvReco_w0->SetMarkerStyle(21);
  me_hEvReco_w0->SetMarkerSize(0.5);
  me_hEvReco_w0->SetMarkerColor(kBlue);
  ee_hEvReco_w0->SetMarkerStyle(21);
  ee_hEvReco_w0->SetMarkerSize(0.5);
  ee_hEvReco_w0->SetMarkerColor(kRed);
  em_hEvReco_w0->SetMarkerStyle(21);
  em_hEvReco_w0->SetMarkerSize(0.5);
  em_hEvReco_w0->SetMarkerColor(kBlue);
  mm_hEvReco_w0->SetMarkerStyle(21);
  mm_hEvReco_w0->SetMarkerSize(0.5);
  mm_hEvReco_w0->SetMarkerColor(kRed);

  me_hEvReco_w3->SetMarkerStyle(21);
  me_hEvReco_w3->SetMarkerSize(0.5);
  me_hEvReco_w3->SetMarkerColor(kBlue);
  ee_hEvReco_w3->SetMarkerStyle(21);
  ee_hEvReco_w3->SetMarkerSize(0.5);
  ee_hEvReco_w3->SetMarkerColor(kRed);
  em_hEvReco_w3->SetMarkerStyle(21);
  em_hEvReco_w3->SetMarkerSize(0.5);
  em_hEvReco_w3->SetMarkerColor(kBlue);
  mm_hEvReco_w3->SetMarkerStyle(21);
  mm_hEvReco_w3->SetMarkerSize(0.5);
  mm_hEvReco_w3->SetMarkerColor(kRed);

  THStack *e_hsEv_w0 = new THStack("e_hsEv_w0","oscillated + unoscillated nue weighted Ev_reco (nu<10.0); Ev_reco (GeV); counts/yr");
  e_hsEv_w0->Add(me_hEvReco_w0);
  e_hsEv_w0->Add(ee_hEvReco_w0);
  THStack *e_hsEv_w3 = new THStack("e_hsEv_w3","oscillated + unoscillated nue weighted Ev_reco (nu<0.3); Ev_reco (GeV); counts/yr");
  e_hsEv_w3->Add(me_hEvReco_w3);
  e_hsEv_w3->Add(ee_hEvReco_w3);
  THStack *m_hsEv_w0 = new THStack("m_hsEv_w0","oscillated + unoscillated numu weighted Ev_reco (nu<10.0); Ev_reco (GeV); counts/yr");
  m_hsEv_w0->Add(em_hEvReco_w0);
  m_hsEv_w0->Add(mm_hEvReco_w0);
  THStack *m_hsEv_w3 = new THStack("m_hsEv_w3","oscillated + unoscillated numu weighted Ev_reco (nu<0.3); Elep_reco (GeV); counts/yr");
  m_hsEv_w3->Add(em_hEvReco_w3);
  m_hsEv_w3->Add(mm_hEvReco_w3);

  TCanvas *cEv_Target = new TCanvas("cEv_Target","",1800,1400);
  cEv_Target->Divide(2,2);
  cEv_Target->cd(1);
  e_hEvReco_w0_ft->Draw();
  e_hsEv_w0->Draw("same");
  TLegend *leEv0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv0->AddEntry(e_hEvReco_w0_ft,"nueCC target");
  leEv0->AddEntry(me_hEvReco_w0,"oscillated mu->e");
  leEv0->AddEntry(ee_hEvReco_w0,"oscilated + unoscillated");
  leEv0->Draw();
  cEv_Target->cd(2);
  e_hEvReco_w3_ft->Draw();
  e_hsEv_w3->Draw("same");
  TLegend *leEv3 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv3->AddEntry(e_hEvReco_w3_ft,"nueCC target");
  leEv3->AddEntry(me_hEvReco_w3,"oscillated mu->e");
  leEv3->AddEntry(ee_hEvReco_w3,"oscilated + unoscillated");
  leEv3->Draw();
  cEv_Target->cd(3);
  m_hEvReco_w0_ft->Draw();
  m_hsEv_w0->Draw("same");
  TLegend *lmEv0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv0->AddEntry(m_hEvReco_w0_ft,"numuCC target");
  lmEv0->AddEntry(em_hEvReco_w0,"oscillated e->mu");
  lmEv0->AddEntry(mm_hEvReco_w0,"oscilated + unoscillated");
  lmEv0->Draw();
  cEv_Target->cd(4);
  m_hEvReco_w3_ft->Draw();
  m_hsEv_w3->Draw("same");
  TLegend *lmEv3 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv3->AddEntry(m_hEvReco_w3_ft,"numuCC target");
  lmEv3->AddEntry(em_hEvReco_w3,"oscillated e->mu");
  lmEv3->AddEntry(mm_hEvReco_w3,"oscilated + unoscillated");
  lmEv3->Draw();
  cEv_Target->SaveAs("true_CC_Ev_target_1.png");

  TCanvas *cEv_Target1 = new TCanvas("cEv_Target1","",1800,1400);
  cEv_Target1->Divide(2,2);
  cEv_Target1->cd(1);
  gPad->SetLogy();
  e_hEvReco_w0_ft->Draw();
  e_hsEv_w0->Draw("same");
  TLegend *leEv01 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv01->AddEntry(e_hEvReco_w0_ft,"nueCC target");
  leEv01->AddEntry(me_hEvReco_w0,"oscillated mu->e");
  leEv01->AddEntry(ee_hEvReco_w0,"oscilated + unoscillated");
  leEv01->Draw();
  cEv_Target1->cd(2);
  gPad->SetLogy();
  e_hEvReco_w3_ft->Draw();
  e_hsEv_w3->Draw("same");
  TLegend *leEv31 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv31->AddEntry(e_hEvReco_w3_ft,"nueCC target");
  leEv31->AddEntry(me_hEvReco_w3,"oscillated mu->e");
  leEv31->AddEntry(ee_hEvReco_w3,"oscilated + unoscillated");
  leEv31->Draw();
  cEv_Target1->cd(3);
  gPad->SetLogy();
  m_hEvReco_w0_ft->Draw();
  m_hsEv_w0->Draw("same");
  TLegend *lmEv01 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv01->AddEntry(m_hEvReco_w0_ft,"numuCC target");
  lmEv01->AddEntry(em_hEvReco_w0,"oscillated e->mu");
  lmEv01->AddEntry(mm_hEvReco_w0,"oscilated + unoscillated");
  lmEv01->Draw();
  cEv_Target1->cd(4);
  gPad->SetLogy();
  m_hEvReco_w3_ft->Draw();
  m_hsEv_w3->Draw("same");
  TLegend *lmEv31 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv31->AddEntry(m_hEvReco_w3_ft,"numuCC target");
  lmEv31->AddEntry(em_hEvReco_w3,"oscillated e->mu");
  lmEv31->AddEntry(mm_hEvReco_w3,"oscilated + unoscillated");
  lmEv31->Draw();
  cEv_Target1->SaveAs("true_CC_Ev_target_1_Log.png");

  e_hElep_w0->SetTitle("nue weighted Elep_reco (nu<10.0)");
  e_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  e_hElep_w3->SetTitle("nue weighted Elep_reco (nu<0.3)");
  e_hElep_w3->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w0->SetTitle("numu weighted Elep_reco (nu<10.0)");
  m_hElep_w0->GetXaxis()->SetTitle("Elep_reco");
  m_hElep_w3->SetTitle("numu weighted Elep_reco (nu<0.3)");
  m_hElep_w3->GetXaxis()->SetTitle("Elep_reco");

  e_hEvReco_w0->SetTitle("nue weighted Ev_reco (nu<10.0)");
  e_hEvReco_w0->GetXaxis()->SetTitle("Ev_reco");
  e_hEvReco_w3->SetTitle("nue weighted Ev_reco (nu<0.3)");
  e_hEvReco_w3->GetXaxis()->SetTitle("Ev_reco");
  m_hEvReco_w0->SetTitle("numu weighted Ev_reco (nu<10.0)");
  m_hEvReco_w0->GetXaxis()->SetTitle("Ev_reco");
  m_hEvReco_w3->SetTitle("numu weighted Ev_reco (nu<0.3)");
  m_hEvReco_w3->GetXaxis()->SetTitle("Ev_reco");

  m_hElepVsEv0->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv3->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  nc_m_hElepVsEv0->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0) without reconstruction cuts");
  nc_m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  nc_m_hElepVsEv3->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3) without reconstruction cuts");
  nc_m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0)");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv3->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",2700,1400);
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
  cElepVsEv->SaveAs("true_CC_Elep_templates_1.png");

  m_hEvRecoVsEv0->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0)");
  m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv3->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3)");
  m_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  nc_m_hEvRecoVsEv0->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0) without reconstruction cuts");
  nc_m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  nc_m_hEvRecoVsEv3->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3) without reconstruction cuts");
  nc_m_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0->SetTitle("CC nue Ev_reco Vs true Ev (nu<10.0)");
  e_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv3->SetTitle("CC nue Ev_reco Vs true Ev (nu<0.3)");
  e_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv = new TCanvas("cEvVsEv","",2700,1400);
  cEvVsEv->Divide(3,2);
  cEvVsEv->cd(1);
  m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(2);
  nc_m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(3);
  e_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(4);
  m_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->cd(5);
  nc_m_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->cd(6);
  e_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->SaveAs("true_CC_Ev_templates_1.png");

  m_hElepVsEv0_cov->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv3_cov->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0_cov->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0)");
  e_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv3_cov->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv3_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov = new TCanvas("cElepVsEv_cov","",1800,1400);
  cElepVsEv_cov->Divide(2,2);
  cElepVsEv_cov->cd(1);
  m_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(2);
  e_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(3);
  m_hElepVsEv3_cov->Draw("colz");
  cElepVsEv_cov->cd(4);
  e_hElepVsEv3_cov->Draw("colz");
  cElepVsEv_cov->SaveAs("cov_CC_Elep_1.png");

  m_hEvRecoVsEv0_cov->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0)");
  m_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv3_cov->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3)");
  m_hEvRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv3_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0_cov->SetTitle("CC nue Ev_reco Vs true Ev (nu<10.0)");
  e_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv3_cov->SetTitle("CC nue Ev_reco Vs true Ev (nu<0.3)");
  e_hEvRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv3_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv_cov = new TCanvas("cEvVsEv_cov","",1800,1400);
  cEvVsEv_cov->Divide(2,2);
  cEvVsEv_cov->cd(1);
  m_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(2);
  e_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(3);
  m_hEvRecoVsEv3_cov->Draw("colz");
  cEvVsEv_cov->cd(4);
  e_hEvRecoVsEv3_cov->Draw("colz");
  cEvVsEv_cov->SaveAs("cov_CC_Ev_1.png");

  TFile *out = new TFile("/dune/app/users/qvuong/lownu/gen_data/CC/output_1.root","RECREATE");
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
  m_hElepVsEv0_cov->Write();
  m_hElepVsEv3_cov->Write();
  e_hElepVsEv0_cov->Write();
  e_hElepVsEv3_cov->Write();

  e_hEvReco_w0->Write();
  e_hEvReco_w3->Write();
  e_hEvReco_w0_ft->Write();
  e_hEvReco_w3_ft->Write();
  me_hEvReco_w0->Write();
  me_hEvReco_w3->Write();
  ee_hEvReco_w0->Write();
  ee_hEvReco_w3->Write();
  m_hEvReco_w0->Write();
  m_hEvReco_w3->Write();
  m_hEvReco_w0_ft->Write();
  m_hEvReco_w3_ft->Write();
  em_hEvReco_w0->Write();
  em_hEvReco_w3->Write();
  mm_hEvReco_w0->Write();
  mm_hEvReco_w3->Write();
  m_hEvRecoVsEv0->Write();
  m_hEvRecoVsEv3->Write();
  nc_m_hEvRecoVsEv0->Write();
  nc_m_hEvRecoVsEv3->Write();
  e_hEvRecoVsEv0->Write();
  e_hEvRecoVsEv3->Write();
  m_hEvRecoVsEv0_cov->Write();
  m_hEvRecoVsEv3_cov->Write();
  e_hEvRecoVsEv0_cov->Write();
  e_hEvRecoVsEv3_cov->Write();
  out->Close();
  
  return(0);

}

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
