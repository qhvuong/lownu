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

int main()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  
  for(int i = 150; i<200; i++){
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

  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,16);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,16);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,16);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,16);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,16);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,16);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,16);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,16);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,16);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0_ft = new TH1D("e_hElep_w0_ft","",100,0,16);
  TH1D *e_hElep_w3_ft = new TH1D("e_hElep_w3_ft","",100,0,16);
  TH1D *m_hElep_w0_ft = new TH1D("m_hElep_w0_ft","",100,0,16);
  TH1D *m_hElep_w3_ft = new TH1D("m_hElep_w3_ft","",100,0,16);

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

  TH2D *m_hElepVsEv0 = new TH2D("m_hElepVsEv0","",19,mubins,nbinsY,yEdges);
  TH2D *m_hElepVsEv3 = new TH2D("m_hElepVsEv3","",19,mubins,nbinsY,yEdges);
  TH2D *nc_m_hElepVsEv0 = new TH2D("nc_m_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepVsEv3 = new TH2D("nc_m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0 = new TH2D("e_hElepVsEv0","",7,ebins,nbinsY,yEdges);
  TH2D *e_hElepVsEv3 = new TH2D("e_hElepVsEv3","",7,ebins,nbinsY,yEdges);

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
  //const int N = 3E7;
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
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  m_hElepVsEv0->SetTitle("numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  m_hElepVsEv3->SetTitle("numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  e_hElepVsEv0->SetTitle("nue Elep_reco Vs true Ev (nu<10.0)");
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
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(3);
  m_hElepVsEv3->Draw("colz");
  cElepVsEv->cd(4);
  e_hElepVsEv3->Draw("colz");
  cElepVsEv->SaveAs("covmtr_CC_ElepVsEv.pdf");
  
  TFile *out = new TFile("/nashome/q/qvuong/gen_data/CC/CC_covmtr.root","RECREATE");
  m_hElepVsEv0->Write();
  m_hElepVsEv3->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv3->Write();
  out->Close();

  return(0);
}

