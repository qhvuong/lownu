#include <iostream>
#include <string>
#include <TF1.h>
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
  TChain * tree = new TChain( "tree", "tree" );
  TChain * meta = new TChain( "meta", "meta" );

  for(int i = 10; i<15; i++){
  tree->Add( Form("/pnfs/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_0%d.root",i) );
  meta->Add( Form("/pnfs/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_0%d.root",i) );
  std::cout << "Nue File number:" << i << "\n";
  }
  
  double total_pot = 0.;
  double pot;
  meta->SetBranchAddress( "pot", &pot );
  const int Nfiles = meta->GetEntries();
  //const int Nfiles = 10000;
  for( int ii = 0; ii < Nfiles; ++ii ) {
    meta->GetEntry(ii);
    total_pot += pot;
  }
  double yrPOT = 1.1E21;
  double scalePOT = yrPOT/total_pot;
  std::cout << "scalePOT = " << scalePOT << "\n";

  //std::cout << total_pot << "\t" << yrPOT*total_pot/1.1E21 << "\t" << scalePOT << "\n";

  TH1D *mm_hElep0_w = new TH1D("mm_hElep0_w","",100,0,16);
  TH1D *mm_hElep1_w = new TH1D("mm_hElep1_w","",100,0,16);
  TH1D *mm_hElep2_w = new TH1D("mm_hElep2_w","",100,0,16);
  TH1D *ee_hElep0_w = new TH1D("ee_hElep0_w","",100,0,16);
  TH1D *ee_hElep1_w = new TH1D("ee_hElep1_w","",100,0,16);
  TH1D *ee_hElep2_w = new TH1D("ee_hElep2_w","",100,0,16);
  TH1D *me_hElep0_w = new TH1D("me_hElep0_w","",100,0,16);
  TH1D *me_hElep1_w = new TH1D("me_hElep1_w","",100,0,16);
  TH1D *me_hElep2_w = new TH1D("me_hElep2_w","",100,0,16);
  TH1D *em_hElep0_w = new TH1D("em_hElep0_w","",100,0,16);
  TH1D *em_hElep1_w = new TH1D("em_hElep1_w","",100,0,16);
  TH1D *em_hElep2_w = new TH1D("em_hElep2_w","",100,0,16);
  TH1D *e_hElep0_w = new TH1D("e_hElep0_w","",100,0,16);
  TH1D *e_hElep1_w = new TH1D("e_hElep1_w","",100,0,16);
  TH1D *e_hElep2_w = new TH1D("e_hElep2_w","",100,0,16);
  TH1D *m_hElep0_w = new TH1D("m_hElep0_w","",100,0,16);
  TH1D *m_hElep1_w = new TH1D("m_hElep1_w","",100,0,16);
  TH1D *m_hElep2_w = new TH1D("m_hElep2_w","",100,0,16);
  TH1D *hElep0_w   = new TH1D("hElep0_w","",100,0,16);
  TH1D *hElep1_w   = new TH1D("hElep1_w","",100,0,16);
  TH1D *hElep2_w   = new TH1D("hElep2_w","",100,0,16);
  TH1D *os_hElep0_w   = new TH1D("os_hElep0_w","",100,0,16);
  TH1D *os_hElep1_w   = new TH1D("os_hElep1_w","",100,0,16);
  TH1D *os_hElep2_w   = new TH1D("os_hElep2_w","",100,0,16);
  TH1D *unos_hElep0_w = new TH1D("unos_hElep0_w","",100,0,16);
  TH1D *unos_hElep1_w = new TH1D("unos_hElep1_w","",100,0,16);
  TH1D *unos_hElep2_w = new TH1D("unos_hElep2_w","",100,0,16);

  TH1D *hElep0_w_ft = new TH1D("hElep0_w_ft","",100,0,16);
  TH1D *hElep1_w_ft = new TH1D("hElep1_w_ft","",100,0,16);
  TH1D *hElep2_w_ft = new TH1D("hElep2_w_ft","",100,0,16);

  TH2D *BthetaVsEe = new TH2D("BthetaVsEe","",400,0.5,60,20,0,100);  

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

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
   
  TH2D *m_hElepVsEv0   = new TH2D("m_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv1   = new TH2D("m_hElepVsEv1","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv2   = new TH2D("m_hElepVsEv2","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv0_w = new TH2D("m_hElepVsEv0_w","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv1_w = new TH2D("m_hElepVsEv1_w","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepVsEv2_w = new TH2D("m_hElepVsEv2_w","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0   = new TH2D("e_hElepVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv1   = new TH2D("e_hElepVsEv1","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv2   = new TH2D("e_hElepVsEv2","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv0_w = new TH2D("e_hElepVsEv0_w","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv1_w = new TH2D("e_hElepVsEv1_w","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepVsEv2_w = new TH2D("e_hElepVsEv2_w","",nbinsX,xEdges,nbinsY,yEdges);

  TH2D *m_hElepVsEv0_cov = new TH2D("m_hElepVsEv0_cov","",19,mubins,100,0,16);
  TH2D *m_hElepVsEv1_cov = new TH2D("m_hElepVsEv1_cov","",19,mubins,100,0,16);
  TH2D *m_hElepVsEv2_cov = new TH2D("m_hElepVsEv2_cov","",19,mubins,100,0,16);
  TH2D *e_hElepVsEv0_cov = new TH2D("e_hElepVsEv0_cov","",7,ebins,100,0,16);
  TH2D *e_hElepVsEv1_cov = new TH2D("e_hElepVsEv1_cov","",7,ebins,100,0,16);
  TH2D *e_hElepVsEv2_cov = new TH2D("e_hElepVsEv2_cov","",7,ebins,100,0,16);

  TH2D *m_hEvRecoVsEv0 = new TH2D("m_hEvRecoVsEv0","",nbinsX,xEdges,100,0,20);  //Etheta2 < 3MeV
  TH2D *m_hEvRecoVsEv1 = new TH2D("m_hEvRecoVsEv1","",nbinsX,xEdges,100,0,20);  //Etheta2 < 1MeV
  TH2D *m_hEvRecoVsEv2 = new TH2D("m_hEvRecoVsEv2","",nbinsX,xEdges,100,0,20);  //Etheta2 < 0.5MeV
  TH2D *e_hEvRecoVsEv0 = new TH2D("e_hEvRecoVsEv0","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv1 = new TH2D("e_hEvRecoVsEv1","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv2 = new TH2D("e_hEvRecoVsEv2","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv0   = new TH2D("hEvRecoVsEv0","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv1   = new TH2D("hEvRecoVsEv1","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv2   = new TH2D("hEvRecoVsEv2","",nbinsX,xEdges,100,0,20);
  TH2D *m_hEvRecoVsEv0_w = new TH2D("m_hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,20);
  TH2D *m_hEvRecoVsEv1_w = new TH2D("m_hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,20);
  TH2D *m_hEvRecoVsEv2_w = new TH2D("m_hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv0_w = new TH2D("e_hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv1_w = new TH2D("e_hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv2_w = new TH2D("e_hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv0_w   = new TH2D("hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv1_w   = new TH2D("hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,20);
  TH2D *hEvRecoVsEv2_w   = new TH2D("hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,20);

  TH1D *hEvRes0 = new TH1D("hEvRes0","",100,-2,2);
  TH1D *hEvRes1 = new TH1D("hEvRes1","",100,-2,2);
  TH1D *hEvRes2 = new TH1D("hEvRes2","",100,-2,2);

  TH1D *m_hEvReco0 = new TH1D("m_hEvReco0","",100,0,20);
  TH1D *m_hEvReco1 = new TH1D("m_hEvReco1","",100,0,20);
  TH1D *m_hEvReco2 = new TH1D("m_hEvReco2","",100,0,20);
  TH1D *mm_hEvReco0 = new TH1D("mm_hEvReco0","",100,0,20);
  TH1D *mm_hEvReco1 = new TH1D("mm_hEvReco1","",100,0,20);
  TH1D *mm_hEvReco2 = new TH1D("mm_hEvReco2","",100,0,20);
  TH1D *me_hEvReco0 = new TH1D("me_hEvReco0","",100,0,20);
  TH1D *me_hEvReco1 = new TH1D("me_hEvReco1","",100,0,20);
  TH1D *me_hEvReco2 = new TH1D("me_hEvReco2","",100,0,20);
  TH1D *e_hEvReco0 = new TH1D("e_hEvReco0","",100,0,20);
  TH1D *e_hEvReco1 = new TH1D("e_hEvReco1","",100,0,20);
  TH1D *e_hEvReco2 = new TH1D("e_hEvReco2","",100,0,20);
  TH1D *ee_hEvReco0 = new TH1D("ee_hEvReco0","",100,0,20);
  TH1D *ee_hEvReco1 = new TH1D("ee_hEvReco1","",100,0,20);
  TH1D *ee_hEvReco2 = new TH1D("ee_hEvReco2","",100,0,20);
  TH1D *em_hEvReco0 = new TH1D("em_hEvReco0","",100,0,20);
  TH1D *em_hEvReco1 = new TH1D("em_hEvReco1","",100,0,20);
  TH1D *em_hEvReco2 = new TH1D("em_hEvReco2","",100,0,20);
  TH1D *hEvReco0 = new TH1D("hEvReco0","",100,0,20);
  TH1D *hEvReco1 = new TH1D("hEvReco1","",100,0,20);
  TH1D *hEvReco2 = new TH1D("hEvReco2","",100,0,20);
  TH1D *os_hEvReco0 = new TH1D("os_hEvReco0","",100,0,20);
  TH1D *os_hEvReco1 = new TH1D("os_hEvReco1","",100,0,20);
  TH1D *os_hEvReco2 = new TH1D("os_hEvReco2","",100,0,20);
  TH1D *unos_hEvReco0 = new TH1D("unos_hEvReco0","",100,0,20);
  TH1D *unos_hEvReco1 = new TH1D("unos_hEvReco1","",100,0,20);
  TH1D *unos_hEvReco2 = new TH1D("unos_hEvReco2","",100,0,20);
  TH1D *hEvReco0_ft = new TH1D("hEvReco0_ft","",100,0,20);
  TH1D *hEvReco1_ft = new TH1D("hEvReco1_ft","",100,0,20);
  TH1D *hEvReco2_ft = new TH1D("hEvReco2_ft","",100,0,20);

  TH2D *m_hEvRecoVsEv0_cov = new TH2D("m_hEvRecoVsEv0_cov","",19,mubins,100,0,20);
  TH2D *m_hEvRecoVsEv1_cov = new TH2D("m_hEvRecoVsEv1_cov","",19,mubins,100,0,20);
  TH2D *m_hEvRecoVsEv2_cov = new TH2D("m_hEvRecoVsEv2_cov","",19,mubins,100,0,20);
  TH2D *e_hEvRecoVsEv0_cov = new TH2D("e_hEvRecoVsEv0_cov","",7,ebins,100,0,20);
  TH2D *e_hEvRecoVsEv1_cov = new TH2D("e_hEvRecoVsEv1_cov","",7,ebins,100,0,20);
  TH2D *e_hEvRecoVsEv2_cov = new TH2D("e_hEvRecoVsEv2_cov","",7,ebins,100,0,20);

  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  double Enu, E[2], p[2]; // the energy of the neutrino, in GeV
  double px[2], py[2], pz[2];
  double best_px[2], best_py[2], best_pz[2];
  int pdg[2];
  double BthetaX, BthetaY, Btheta;
  double PthetaX, PthetaY, Ptheta;
  double BthetaX_sm, BthetaY_sm, Btheta_sm;

  tree->SetBranchAddress( "Enu", &Enu );
  tree->SetBranchAddress( "E", &E );
  tree->SetBranchAddress( "pdg", &pdg );
  tree->SetBranchAddress( "px", &px );
  tree->SetBranchAddress( "py", &py );
  tree->SetBranchAddress( "pz", &pz );
  tree->SetBranchAddress( "best_px", &best_px );
  tree->SetBranchAddress( "best_py", &best_py );
  tree->SetBranchAddress( "best_pz", &best_pz );

  TF1 *tsmear1 = new TF1( "tsmear1", "3.29 + 3.485*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmear2 = new TF1( "tsmear2", "10.287 + 4.889*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmearRatio = new TF1( "tsmearRatio", "0.039 + 0.551*pow(x,-1.) - 0.268*pow(x,-0.5)", 0., 999.9 );
  TF1 *doubleGaus = new TF1( "dg", "[0]*TMath::Exp(-0.5*pow(x/[1],2)) + [2]*TMath::Exp(-0.5*pow(x/[3],2))", -1000., 1000. );

  TRandom3 *rando = new TRandom3(8888);

  double s2tW = 0.23;
  double m_C_LL = -1./2. + s2tW; 
  double m_C_LR = s2tW;
  double e_C_LL = 1./2. + s2tW;
  double e_C_LR = s2tW;
  double y, sigma_m, sigma_e;

  double Uee2, Umm2, dm2;

  int para = 1;
  if(para == 1) {
  Uee2 = 0.01;
  Umm2 = 0.0016;
  dm2 = 1.3; //eV2
  }
  if(para == 2) {
  Uee2 = 0.04;
  Umm2 = 0.01;
  dm2 = 6.0; //eV2
  }
  double s2ee2 = 4*Uee2*(1 - Uee2); 
  double s2mm2 = 4*Umm2*(1 - Umm2); 
  double s2me2 = 4*Uee2*Umm2;
  double L = 500; //m
  double del, pmumu, pmue, pee;

  double me = 510;  //keV
  double Ev_reco;

  const int N = tree->GetEntries();
  //const int N = 10000;
  //scalePOT = 10;

  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 1000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);


    //for(int i=0; i<2; i++){
    if( pdg[1]==14 || pdg[1]==12 ){
    //if(E[i]<.5) continue;

    BthetaX = atan2(best_px[0],best_pz[0])*1E3;
    BthetaY = atan2(best_py[0],best_pz[0])*1E3;
    Btheta = atan2(sqrt(best_py[0]*best_py[0]+best_px[0]*best_px[0]),best_pz[0])*1E3;

    double p1 = tsmear1->Eval(E[0]);
    double p3 = tsmear2->Eval(E[0]);
    double ratio = tsmearRatio->Eval(E[0]);

    doubleGaus->SetParameter( 1, p1 );
    doubleGaus->SetParameter( 3, p3 );
    doubleGaus->SetParameter( 0, 1. );
    doubleGaus->SetParameter( 2, ratio );

    BthetaX_sm = BthetaX + doubleGaus->GetRandom();
    BthetaY_sm = BthetaY + doubleGaus->GetRandom();
    Btheta_sm = atan(sqrt(tan(BthetaX_sm/1E3)*tan(BthetaX_sm/1E3) + tan(BthetaY_sm/1E3)*tan(BthetaY_sm/1E3)))*1E3;    
      
    del = 1.27/1E3*dm2*L/Enu;
    pmue = s2me2 * pow(sin(del),2);
    pmumu = 1 - s2mm2 * pow(sin(del),2);   
    pee = 1 - s2ee2 * pow(sin(del),2);

    Ev_reco = E[0]/(1-E[0]*Btheta_sm*Btheta_sm/(2*me));

    for(int i=0; i<2; i++){
      p[i]=sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
    }
    y = E[0]/Enu;
    sigma_m = (m_C_LL*m_C_LL + m_C_LR*m_C_LR*(1.-y)*(1.-y));
    sigma_e = (e_C_LL*e_C_LL + e_C_LR*e_C_LR*(1.-y)*(1.-y));

    BthetaVsEe->Fill(E[0],Btheta_sm);  
  
    if(E[0]*Btheta_sm*Btheta_sm/1E3<3.){
      hEvRes0->Fill((Ev_reco-Enu)/Enu);
      if(pdg[1] == 14){
        //muon to muon ("survival")
        mm_hElep0_w->Fill(E[0],pmumu);
        mm_hEvReco0->Fill(Ev_reco,pmumu);
        
        //muon to electron
        me_hElep0_w->Fill(E[0],pmue*sigma_e/sigma_m);
        me_hEvReco0->Fill(Ev_reco,pmue*sigma_e/sigma_m);

        //muon template weighted     
        m_hElepVsEv0_w->Fill(Enu,E[0],sigma_e/sigma_m);

        m_hEvRecoVsEv0_w->Fill(Enu,Ev_reco,sigma_e/sigma_m);

        //muon template     
        m_hElepVsEv0->Fill(Enu,E[0]);
        m_hElepVsEv0_cov->Fill(Enu,E[0]);

        m_hEvRecoVsEv0->Fill(Enu,Ev_reco);
        m_hEvRecoVsEv0_cov->Fill(Enu,Ev_reco);
      }

      if(pdg[1] == 12){
        //electron to electron ("survival")
        ee_hElep0_w->Fill(E[0],pee);
        ee_hEvReco0->Fill(Ev_reco,pee);
      
        //electron to muon 
        em_hElep0_w->Fill(E[0],pmue*sigma_m/sigma_e);
        em_hEvReco0->Fill(Ev_reco,pmue*sigma_m/sigma_e);
      
        //electron template weighted
        e_hElepVsEv0_w->Fill(Enu,E[0],sigma_m/sigma_e);

        e_hEvRecoVsEv0_w->Fill(Enu,Ev_reco,sigma_m/sigma_e);

        //electron template
        e_hElepVsEv0->Fill(Enu,E[0]);
        e_hElepVsEv0_cov->Fill(Enu,E[0]);

        e_hEvRecoVsEv0->Fill(Enu,Ev_reco);
        e_hEvRecoVsEv0_cov->Fill(Enu,Ev_reco);
      }
    }

    if(E[0]*Btheta_sm*Btheta_sm/1E3<0.8){
      hEvRes1->Fill((Ev_reco-Enu)/Enu);
      if(pdg[1] == 14){
        //muon to muon ("survival")
        mm_hElep1_w->Fill(E[0],pmumu);
        mm_hEvReco1->Fill(Ev_reco,pmumu);
        
        //muon to electron
        me_hElep1_w->Fill(E[0],pmue*sigma_e/sigma_m);
        me_hEvReco1->Fill(Ev_reco,pmue*sigma_e/sigma_m);

        //muon template weighted     
        m_hElepVsEv1_w->Fill(Enu,E[0],sigma_e/sigma_m);

        m_hEvRecoVsEv1_w->Fill(Enu,Ev_reco,sigma_e/sigma_m);

        //muon template     
        m_hElepVsEv1->Fill(Enu,E[0]);
        m_hElepVsEv1_cov->Fill(Enu,E[0]);

        m_hEvRecoVsEv1->Fill(Enu,Ev_reco);
        m_hEvRecoVsEv1_cov->Fill(Enu,Ev_reco);
      }

      if(pdg[1] == 12){
        //electron to electron ("survival")
        ee_hElep1_w->Fill(E[0],pee);
        ee_hEvReco1->Fill(Ev_reco,pee);
      
        //electron to muon 
        em_hElep1_w->Fill(E[0],pmue*sigma_m/sigma_e);
        em_hEvReco1->Fill(Ev_reco,pmue*sigma_m/sigma_e);
      
        //electron template weighted
        e_hElepVsEv1_w->Fill(Enu,E[0],sigma_m/sigma_e);

        e_hEvRecoVsEv1_w->Fill(Enu,Ev_reco,sigma_m/sigma_e);

        //electron template
        e_hElepVsEv1->Fill(Enu,E[0]);
        e_hElepVsEv1_cov->Fill(Enu,E[0]);

        e_hEvRecoVsEv1->Fill(Enu,Ev_reco);
        e_hEvRecoVsEv1_cov->Fill(Enu,Ev_reco);
      }
    }

    if(E[0]*Btheta_sm*Btheta_sm/1E3<0.5){
      hEvRes2->Fill((Ev_reco-Enu)/Enu);
      if(pdg[1] == 14){
        //muon to muon ("survival")
        mm_hElep2_w->Fill(E[0],pmumu);
        mm_hEvReco2->Fill(Ev_reco,pmumu);
        
        //muon to electron
        me_hElep2_w->Fill(E[0],pmue*sigma_e/sigma_m);
        me_hEvReco2->Fill(Ev_reco,pmue*sigma_e/sigma_m);

        //muon template weighted     
        m_hElepVsEv2_w->Fill(Enu,E[0],sigma_e/sigma_m);

        m_hEvRecoVsEv2_w->Fill(Enu,Ev_reco,sigma_e/sigma_m);

        //muon template     
        m_hElepVsEv2->Fill(Enu,E[0]);
        m_hElepVsEv2_cov->Fill(Enu,E[0]);

        m_hEvRecoVsEv2->Fill(Enu,Ev_reco);
        m_hEvRecoVsEv2_cov->Fill(Enu,Ev_reco);
      }

      if(pdg[1] == 12){
        //electron to electron ("survival")
        ee_hElep2_w->Fill(E[0],pee);
        ee_hEvReco2->Fill(Ev_reco,pee);
      
        //electron to muon 
        em_hElep2_w->Fill(E[0],pmue*sigma_m/sigma_e);
        em_hEvReco2->Fill(Ev_reco,pmue*sigma_m/sigma_e);
      
        //electron template weighted
        e_hElepVsEv2_w->Fill(Enu,E[0],sigma_m/sigma_e);

        e_hEvRecoVsEv2_w->Fill(Enu,Ev_reco,sigma_m/sigma_e);

        //electron template
        e_hElepVsEv2->Fill(Enu,E[0]);
        e_hElepVsEv2_cov->Fill(Enu,E[0]);

        e_hEvRecoVsEv2->Fill(Enu,Ev_reco);
        e_hEvRecoVsEv2_cov->Fill(Enu,Ev_reco);
      }
    }
  }
  }

  hEvRes0->Scale(scalePOT);
  hEvRes1->Scale(scalePOT);
  hEvRes2->Scale(scalePOT);

  me_hElep0_w->Scale(scalePOT); mm_hElep0_w->Scale(scalePOT);
  ee_hElep0_w->Scale(scalePOT); em_hElep0_w->Scale(scalePOT);
  me_hElep1_w->Scale(scalePOT); mm_hElep1_w->Scale(scalePOT);
  ee_hElep1_w->Scale(scalePOT); em_hElep1_w->Scale(scalePOT);
  me_hElep2_w->Scale(scalePOT); mm_hElep2_w->Scale(scalePOT);
  ee_hElep2_w->Scale(scalePOT); em_hElep2_w->Scale(scalePOT);

  me_hEvReco0->Scale(scalePOT); mm_hEvReco0->Scale(scalePOT);
  ee_hEvReco0->Scale(scalePOT); em_hEvReco0->Scale(scalePOT);
  me_hEvReco1->Scale(scalePOT); mm_hEvReco1->Scale(scalePOT);
  ee_hEvReco1->Scale(scalePOT); em_hEvReco1->Scale(scalePOT);
  me_hEvReco2->Scale(scalePOT); mm_hEvReco2->Scale(scalePOT);
  ee_hEvReco2->Scale(scalePOT); em_hEvReco2->Scale(scalePOT);

  m_hElepVsEv0_w->Scale(scalePOT); e_hElepVsEv0_w->Scale(scalePOT);
  m_hElepVsEv0->Scale(scalePOT);   e_hElepVsEv0->Scale(scalePOT);
  m_hElepVsEv1_w->Scale(scalePOT); e_hElepVsEv1_w->Scale(scalePOT);
  m_hElepVsEv1->Scale(scalePOT);   e_hElepVsEv1->Scale(scalePOT);
  m_hElepVsEv2_w->Scale(scalePOT); e_hElepVsEv2_w->Scale(scalePOT);
  m_hElepVsEv2->Scale(scalePOT);   e_hElepVsEv2->Scale(scalePOT);

  m_hEvRecoVsEv0_w->Scale(scalePOT); e_hEvRecoVsEv0_w->Scale(scalePOT);
  m_hEvRecoVsEv0->Scale(scalePOT);   e_hEvRecoVsEv0->Scale(scalePOT);
  m_hEvRecoVsEv1_w->Scale(scalePOT); e_hEvRecoVsEv1_w->Scale(scalePOT);
  m_hEvRecoVsEv1->Scale(scalePOT);   e_hEvRecoVsEv1->Scale(scalePOT);
  m_hEvRecoVsEv2_w->Scale(scalePOT); e_hEvRecoVsEv2_w->Scale(scalePOT);
  m_hEvRecoVsEv2->Scale(scalePOT);   e_hEvRecoVsEv2->Scale(scalePOT);

  m_hElepVsEv0_cov->Scale(scalePOT);   e_hElepVsEv0_cov->Scale(scalePOT);
  m_hElepVsEv1_cov->Scale(scalePOT);   e_hElepVsEv1_cov->Scale(scalePOT);
  m_hElepVsEv2_cov->Scale(scalePOT);   e_hElepVsEv2_cov->Scale(scalePOT);
  m_hEvRecoVsEv0_cov->Scale(scalePOT); e_hEvRecoVsEv0_cov->Scale(scalePOT);
  m_hEvRecoVsEv1_cov->Scale(scalePOT); e_hEvRecoVsEv1_cov->Scale(scalePOT);
  m_hEvRecoVsEv2_cov->Scale(scalePOT); e_hEvRecoVsEv2_cov->Scale(scalePOT);

  e_hElep0_w->Add(me_hElep0_w);    e_hElep0_w->Add(ee_hElep0_w);
  m_hElep0_w->Add(em_hElep0_w);    m_hElep0_w->Add(mm_hElep0_w);
  hElep0_w->Add(e_hElep0_w);       hElep0_w->Add(m_hElep0_w);  
  os_hElep0_w->Add(me_hElep0_w);   os_hElep0_w->Add(em_hElep0_w);
  unos_hElep0_w->Add(ee_hElep0_w); unos_hElep0_w->Add(mm_hElep0_w);

  e_hElep1_w->Add(me_hElep1_w);    e_hElep1_w->Add(ee_hElep1_w);
  m_hElep1_w->Add(em_hElep1_w);    m_hElep1_w->Add(mm_hElep1_w);
  hElep1_w->Add(e_hElep1_w);       hElep1_w->Add(m_hElep1_w);  
  os_hElep1_w->Add(me_hElep1_w);   os_hElep1_w->Add(em_hElep1_w);
  unos_hElep1_w->Add(ee_hElep1_w); unos_hElep1_w->Add(mm_hElep1_w);

  e_hElep2_w->Add(me_hElep2_w);    e_hElep2_w->Add(ee_hElep2_w);
  m_hElep2_w->Add(em_hElep2_w);    m_hElep2_w->Add(mm_hElep2_w);
  hElep2_w->Add(e_hElep2_w);       hElep2_w->Add(m_hElep2_w);  
  os_hElep2_w->Add(me_hElep2_w);   os_hElep2_w->Add(em_hElep2_w);
  unos_hElep2_w->Add(ee_hElep2_w); unos_hElep2_w->Add(mm_hElep2_w);

  e_hEvReco0->Add(me_hEvReco0);    e_hEvReco0->Add(ee_hEvReco0);
  m_hEvReco0->Add(em_hEvReco0);    m_hEvReco0->Add(mm_hEvReco0);
  hEvReco0->Add(e_hEvReco0);       hEvReco0->Add(m_hEvReco0);
  os_hEvReco0->Add(me_hEvReco0);   os_hEvReco0->Add(em_hEvReco0);
  unos_hEvReco0->Add(ee_hEvReco0); unos_hEvReco0->Add(mm_hEvReco0);

  e_hEvReco1->Add(me_hEvReco1);    e_hEvReco1->Add(ee_hEvReco1);
  m_hEvReco1->Add(em_hEvReco1);    m_hEvReco1->Add(mm_hEvReco1);
  hEvReco1->Add(e_hEvReco1);       hEvReco1->Add(m_hEvReco1);
  os_hEvReco1->Add(me_hEvReco1);   os_hEvReco1->Add(em_hEvReco1);
  unos_hEvReco1->Add(ee_hEvReco1); unos_hEvReco1->Add(mm_hEvReco1);

  e_hEvReco2->Add(me_hEvReco2);    e_hEvReco2->Add(ee_hEvReco2);
  m_hEvReco2->Add(em_hEvReco2);    m_hEvReco2->Add(mm_hEvReco2);
  hEvReco2->Add(e_hEvReco2);       hEvReco2->Add(m_hEvReco2);
  os_hEvReco2->Add(me_hEvReco2);   os_hEvReco2->Add(em_hEvReco2);
  unos_hEvReco2->Add(ee_hEvReco2); unos_hEvReco2->Add(mm_hEvReco2);

  for( int bx = 1; bx <= hElep0_w->GetNbinsX(); bx++ ){
    double mean_Elep0 = hElep0_w->GetBinContent(bx);
    double mean_Elep1 = hElep1_w->GetBinContent(bx);
    double mean_Elep2 = hElep2_w->GetBinContent(bx);
    double ft_Elep0 = rando->Poisson(mean_Elep0);
    double ft_Elep1 = rando->Poisson(mean_Elep1);
    double ft_Elep2 = rando->Poisson(mean_Elep2);
    hElep0_w_ft->AddBinContent(bx, ft_Elep0);
    hElep1_w_ft->AddBinContent(bx, ft_Elep1);
    hElep2_w_ft->AddBinContent(bx, ft_Elep2);
  }
  for( int bx = 1; bx <= hEvReco0->GetNbinsX(); bx++ ){
    double mean_Ev0 = hEvReco0->GetBinContent(bx);
    double mean_Ev1 = hEvReco1->GetBinContent(bx);
    double mean_Ev2 = hEvReco2->GetBinContent(bx);
    double ft_Ev0 = rando->Poisson(mean_Ev0);
    double ft_Ev1 = rando->Poisson(mean_Ev1);
    double ft_Ev2 = rando->Poisson(mean_Ev2);
    hEvReco0_ft->AddBinContent(bx, ft_Ev0);
    hEvReco1_ft->AddBinContent(bx, ft_Ev1);
    hEvReco2_ft->AddBinContent(bx, ft_Ev2);
  }

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  e_hElep0_w->SetStats(0);
  e_hElep0_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 3.0MeV)");
  e_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_w->SetMarkerStyle(8);
  e_hElep0_w->SetMarkerSize(1.0);
  e_hElep0_w->SetMarkerColor(kBlack);
  m_hElep0_w->SetStats(0);
  m_hElep0_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 3.0MeV)");
  m_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_w->SetMarkerStyle(8);
  m_hElep0_w->SetMarkerSize(1.0);
  m_hElep0_w->SetMarkerColor(kBlack);
  e_hElep1_w->SetStats(0);
  e_hElep1_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 0.8MeV)");
  e_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep1_w->SetMarkerStyle(8);
  e_hElep1_w->SetMarkerSize(1.0);
  e_hElep1_w->SetMarkerColor(kBlack);
  m_hElep1_w->SetStats(0);
  m_hElep1_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 0.8MeV)");
  m_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep1_w->SetMarkerStyle(8);
  m_hElep1_w->SetMarkerSize(1.0);
  m_hElep1_w->SetMarkerColor(kBlack);
  e_hElep2_w->SetStats(0);
  e_hElep2_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 0.5MeV)");
  e_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep2_w->SetMarkerStyle(8);
  e_hElep2_w->SetMarkerSize(1.0);
  e_hElep2_w->SetMarkerColor(kBlack);
  m_hElep2_w->SetStats(0);
  m_hElep2_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 0.5MeV)");
  m_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep2_w->SetMarkerStyle(8);
  m_hElep2_w->SetMarkerSize(1.0);
  m_hElep2_w->SetMarkerColor(kBlack);

  hElep0_w_ft->SetStats(0);
  hElep0_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 3.0MeV)");
  hElep0_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hElep0_w_ft->SetLineColor(kBlack);
  hElep1_w_ft->SetStats(0);
  hElep1_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 0.8MeV)");
  hElep1_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hElep1_w_ft->SetLineColor(kBlack);
  hElep2_w_ft->SetStats(0);
  hElep2_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 0.5MeV)");
  hElep2_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hElep2_w_ft->SetLineColor(kBlack);

  me_hElep0_w->SetStats(0);
  me_hElep0_w->SetMarkerStyle(21);
  me_hElep0_w->SetMarkerSize(0.5);
  me_hElep0_w->SetMarkerColor(kBlue);
  ee_hElep0_w->SetStats(0);
  ee_hElep0_w->SetMarkerStyle(21);
  ee_hElep0_w->SetMarkerSize(0.5);
  ee_hElep0_w->SetMarkerColor(kRed);
  em_hElep0_w->SetStats(0);
  em_hElep0_w->SetMarkerStyle(21);
  em_hElep0_w->SetMarkerSize(0.5);
  em_hElep0_w->SetMarkerColor(kBlue);
  mm_hElep0_w->SetStats(0);
  mm_hElep0_w->SetMarkerStyle(21);
  mm_hElep0_w->SetMarkerSize(0.5);
  mm_hElep0_w->SetMarkerColor(kRed);

  me_hElep1_w->SetStats(0);
  me_hElep1_w->SetMarkerStyle(21);
  me_hElep1_w->SetMarkerSize(0.5);
  me_hElep1_w->SetMarkerColor(kBlue);
  ee_hElep1_w->SetStats(0);
  ee_hElep1_w->SetMarkerStyle(21);
  ee_hElep1_w->SetMarkerSize(0.5);
  ee_hElep1_w->SetMarkerColor(kRed);
  em_hElep1_w->SetStats(0);
  em_hElep1_w->SetMarkerStyle(21);
  em_hElep1_w->SetMarkerSize(0.5);
  em_hElep1_w->SetMarkerColor(kBlue);
  mm_hElep1_w->SetStats(0);
  mm_hElep1_w->SetMarkerStyle(21);
  mm_hElep1_w->SetMarkerSize(0.5);
  mm_hElep1_w->SetMarkerColor(kRed);

  me_hElep2_w->SetStats(0);
  me_hElep2_w->SetMarkerStyle(21);
  me_hElep2_w->SetMarkerSize(0.5);
  me_hElep2_w->SetMarkerColor(kBlue);
  ee_hElep2_w->SetStats(0);
  ee_hElep2_w->SetMarkerStyle(21);
  ee_hElep2_w->SetMarkerSize(0.5);
  ee_hElep2_w->SetMarkerColor(kRed);
  em_hElep2_w->SetStats(0);
  em_hElep2_w->SetMarkerStyle(21);
  em_hElep2_w->SetMarkerSize(0.5);
  em_hElep2_w->SetMarkerColor(kBlue);
  mm_hElep2_w->SetStats(0);
  mm_hElep2_w->SetMarkerStyle(21);
  mm_hElep2_w->SetMarkerSize(0.5);
  mm_hElep2_w->SetMarkerColor(kRed);

  THStack *e_hsElep0_w = new THStack("e_hsElep0_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 3.0MeV); Elep_reco (GeV); count");
  e_hsElep0_w->Add(me_hElep0_w);
  e_hsElep0_w->Add(ee_hElep0_w);
  THStack *m_hsElep0_w = new THStack("m_hsElep0_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 3.0MeV); Elep_reco (GeV); count");
  m_hsElep0_w->Add(em_hElep0_w);
  m_hsElep0_w->Add(mm_hElep0_w);
  THStack *e_hsElep1_w = new THStack("e_hsElep1_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 0.8MeV); Elep_reco (GeV); count");
  e_hsElep1_w->Add(me_hElep1_w);
  e_hsElep1_w->Add(ee_hElep1_w);
  THStack *m_hsElep1_w = new THStack("m_hsElep1_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 0.8MeV); Elep_reco (GeV); count");
  m_hsElep1_w->Add(em_hElep1_w);
  m_hsElep1_w->Add(mm_hElep1_w);
  THStack *e_hsElep2_w = new THStack("e_hsElep2_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 0.5MeV); Elep_reco (GeV); count");
  e_hsElep2_w->Add(me_hElep2_w);
  e_hsElep2_w->Add(ee_hElep2_w);
  THStack *m_hsElep2_w = new THStack("m_hsElep2_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 0.5MeV); Elep_reco (GeV); count");
  m_hsElep2_w->Add(em_hElep2_w);
  m_hsElep2_w->Add(mm_hElep2_w);

  TCanvas *cTarget_Elep = new TCanvas("cTarget_Elep","",2700,1400);
  cTarget_Elep->Divide(3,2);
  cTarget_Elep->cd(1);
  e_hElep0_w->Draw();
  e_hsElep0_w->Draw("same");
  TLegend *le0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  le0->AddEntry(e_hElep0_w,"true target");
  le0->AddEntry(me_hElep0_w,"true oscillated mu->e");
  le0->AddEntry(ee_hElep0_w,"true oscillated + unoscillated");
  le0->Draw();
  cTarget_Elep->cd(2);
  e_hElep1_w->Draw();
  e_hsElep1_w->Draw("same");
  TLegend *le1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  le1->AddEntry(e_hElep1_w,"true target");
  le1->AddEntry(me_hElep1_w,"true oscillated mu->e");
  le1->AddEntry(ee_hElep1_w,"true oscillated + unoscillated");
  le1->Draw();
  cTarget_Elep->cd(3);
  e_hElep2_w->Draw();
  e_hsElep2_w->Draw("same");
  TLegend *le2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  le2->AddEntry(e_hElep2_w,"true target");
  le2->AddEntry(me_hElep2_w,"true oscillated mu->e");
  le2->AddEntry(ee_hElep2_w,"true oscillated + unoscillated");
  le2->Draw();
  cTarget_Elep->cd(4);
  m_hElep0_w->Draw();
  m_hsElep0_w->Draw("same");
  TLegend *lm0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lm0->AddEntry(m_hElep0_w,"true target");
  lm0->AddEntry(em_hElep0_w,"true oscillated e->mu");
  lm0->AddEntry(mm_hElep0_w,"true oscillated + unoscillated");
  lm0->Draw();
  cTarget_Elep->cd(5);
  m_hElep1_w->Draw();
  m_hsElep1_w->Draw("same");
  TLegend *lm1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lm1->AddEntry(m_hElep1_w,"true target");
  lm1->AddEntry(em_hElep1_w,"true oscillated e->mu");
  lm1->AddEntry(mm_hElep1_w,"true oscillated + unoscillated");
  lm1->Draw();
  cTarget_Elep->cd(6);
  m_hElep2_w->Draw();
  m_hsElep2_w->Draw("same");
  TLegend *lm2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lm2->AddEntry(m_hElep2_w,"true target");
  lm2->AddEntry(em_hElep2_w,"true oscillated e->mu");
  lm2->AddEntry(mm_hElep2_w,"true oscillated + unoscillated");
  lm2->Draw();
  cTarget_Elep->SaveAs(Form("true_nue_Elep_target_%d.png",para));

  os_hElep0_w->SetStats(0);
  os_hElep0_w->SetMarkerStyle(21);
  os_hElep0_w->SetMarkerSize(0.5);
  os_hElep0_w->SetMarkerColor(kBlue);
  unos_hElep0_w->SetStats(0);
  unos_hElep0_w->SetMarkerStyle(21);
  unos_hElep0_w->SetMarkerSize(0.5);
  unos_hElep0_w->SetMarkerColor(kRed);
  hElep0_w->SetStats(0);
  hElep0_w->SetMarkerStyle(21);
  hElep0_w->SetMarkerSize(0.5);
  hElep0_w->SetMarkerColor(kRed);

  os_hElep1_w->SetStats(0);
  os_hElep1_w->SetMarkerStyle(21);
  os_hElep1_w->SetMarkerSize(0.5);
  os_hElep1_w->SetMarkerColor(kBlue);
  unos_hElep1_w->SetStats(0);
  unos_hElep1_w->SetMarkerStyle(21);
  unos_hElep1_w->SetMarkerSize(0.5);
  unos_hElep1_w->SetMarkerColor(kRed);
  hElep1_w->SetStats(0);
  hElep1_w->SetMarkerStyle(21);
  hElep1_w->SetMarkerSize(0.5);
  hElep1_w->SetMarkerColor(kRed);

  os_hElep2_w->SetStats(0);
  os_hElep2_w->SetMarkerStyle(21);
  os_hElep2_w->SetMarkerSize(0.5);
  os_hElep2_w->SetMarkerColor(kBlue);
  unos_hElep2_w->SetStats(0);
  unos_hElep2_w->SetMarkerStyle(21);
  unos_hElep2_w->SetMarkerSize(0.5);
  unos_hElep2_w->SetMarkerColor(kRed);
  hElep2_w->SetStats(0);
  hElep2_w->SetMarkerStyle(21);
  hElep2_w->SetMarkerSize(0.5);
  hElep2_w->SetMarkerColor(kRed);

  TCanvas *cElep_w_sum = new TCanvas("cElep_w_sum","",2700,700);
  cElep_w_sum->Divide(3,1);
  cElep_w_sum->cd(1);
  hElep0_w_ft->Draw();
  os_hElep0_w->Draw("same");
  hElep0_w->Draw("same");
  TLegend *l0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  l0->AddEntry(hElep0_w_ft,"fluctuated target");
  l0->AddEntry(os_hElep0_w,"true oscillated mu->e & e->mu");
  l0->AddEntry(hElep0_w,"true target");
  l0->Draw();
  cElep_w_sum->cd(2);
  hElep1_w_ft->Draw();
  os_hElep1_w->Draw("same");
  hElep1_w->Draw("same");
  TLegend *l1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  l1->AddEntry(hElep1_w_ft,"fluctuated target");
  l1->AddEntry(os_hElep1_w,"true oscillated mu->e & e->mu");
  l1->AddEntry(hElep1_w,"true target");
  l1->Draw();
  cElep_w_sum->cd(3);
  hElep2_w_ft->Draw();
  os_hElep2_w->Draw("same");
  hElep2_w->Draw("same");
  TLegend *l2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  l2->AddEntry(hElep2_w_ft,"fluctuated target");
  l2->AddEntry(os_hElep2_w,"true oscillated mu->e & e->mu");
  l2->AddEntry(hElep2_w,"true target");
  l2->Draw();
  cElep_w_sum->SaveAs(Form("true_nue_Elep_target_%d_sum.png",para));


  m_hElepVsEv0->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv1->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.8MeV)");
  m_hElepVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv1->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv1->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.8MeV)");
  e_hElepVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv1->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv2->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hElepVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv2->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hElepVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",2700,1400);
  cElepVsEv->Divide(3,2);
  cElepVsEv->cd(1);
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  e_hElepVsEv1->Draw("colz");
  cElepVsEv->cd(3);
  e_hElepVsEv2->Draw("colz");
  cElepVsEv->cd(4);
  m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(5);
  m_hElepVsEv1->Draw("colz");
  cElepVsEv->cd(6);
  m_hElepVsEv2->Draw("colz");
  cElepVsEv->SaveAs(Form("true_nue_Elep_templates_%d.png",para));


  m_hElepVsEv0_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv1_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.8MeV)");
  m_hElepVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv1_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv1_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.8MeV)");
  e_hElepVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv1_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv2_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hElepVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv2_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv2_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hElepVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv2_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov = new TCanvas("cElepVsEv_cov","",2700,1400);
  cElepVsEv_cov->Divide(3,2);
  cElepVsEv_cov->cd(1);
  e_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(2);
  e_hElepVsEv1_cov->Draw("colz");
  cElepVsEv_cov->cd(3);
  e_hElepVsEv2_cov->Draw("colz");
  cElepVsEv_cov->cd(4);
  m_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(5);
  m_hElepVsEv1_cov->Draw("colz");
  cElepVsEv_cov->cd(6);
  m_hElepVsEv2_cov->Draw("colz");
  cElepVsEv_cov->SaveAs(Form("cov_nue_Elep_%d.png",para));


  e_hEvReco0->SetStats(0);
  e_hEvReco0->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 3.0MeV)");
  e_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco0->SetMarkerStyle(8);
  e_hEvReco0->SetMarkerSize(1.0);
  e_hEvReco0->SetMarkerColor(kBlack);
  m_hEvReco0->SetStats(0);
  m_hEvReco0->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 3.0MeV)");
  m_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco0->SetMarkerStyle(8);
  m_hEvReco0->SetMarkerSize(1.0);
  m_hEvReco0->SetMarkerColor(kBlack);
  e_hEvReco1->SetStats(0);
  e_hEvReco1->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 0.8MeV)");
  e_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco1->SetMarkerStyle(8);
  e_hEvReco1->SetMarkerSize(1.0);
  e_hEvReco1->SetMarkerColor(kBlack);
  m_hEvReco1->SetStats(0);
  m_hEvReco1->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 0.8MeV)");
  m_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco1->SetMarkerStyle(8);
  m_hEvReco1->SetMarkerSize(1.0);
  m_hEvReco1->SetMarkerColor(kBlack);
  e_hEvReco2->SetStats(0);
  e_hEvReco2->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 0.5MeV)");
  e_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco2->SetMarkerStyle(8);
  e_hEvReco2->SetMarkerSize(1.0);
  e_hEvReco2->SetMarkerColor(kBlack);
  m_hEvReco2->SetStats(0);
  m_hEvReco2->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 0.5MeV)");
  m_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco2->SetMarkerStyle(8);
  m_hEvReco2->SetMarkerSize(1.0);
  m_hEvReco2->SetMarkerColor(kBlack);


  hEvReco0_ft->SetStats(0);
  hEvReco0_ft->SetTitle("fluctuated nu+e Ev_reco target (Etheta2 < 3.0MeV)");
  hEvReco0_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hEvReco0_ft->SetLineColor(kBlack);
  hEvReco1_ft->SetStats(0);
  hEvReco1_ft->SetTitle("fluctuated nu+e Ev_reco target (Etheta2 < 0.8MeV)");
  hEvReco1_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hEvReco1_ft->SetLineColor(kBlack);
  hEvReco2_ft->SetStats(0);
  hEvReco2_ft->SetTitle("fluctuated nu+e Ev_reco target (Etheta2 < 0.5MeV)");
  hEvReco2_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hEvReco2->SetMarkerStyle(8);
  hEvReco2->SetMarkerSize(1.0);
  hEvReco2->SetMarkerColor(kRed);
  hEvReco2_ft->SetLineColor(kBlack);

  me_hEvReco0->SetStats(0);
  me_hEvReco0->SetMarkerStyle(21);
  me_hEvReco0->SetMarkerSize(0.5);
  me_hEvReco0->SetMarkerColor(kBlue);
  ee_hEvReco0->SetStats(0);
  ee_hEvReco0->SetMarkerStyle(21);
  ee_hEvReco0->SetMarkerSize(0.5);
  ee_hEvReco0->SetMarkerColor(kRed);
  em_hEvReco0->SetStats(0);
  em_hEvReco0->SetMarkerStyle(21);
  em_hEvReco0->SetMarkerSize(0.5);
  em_hEvReco0->SetMarkerColor(kBlue);
  mm_hEvReco0->SetStats(0);
  mm_hEvReco0->SetMarkerStyle(21);
  mm_hEvReco0->SetMarkerSize(0.5);
  mm_hEvReco0->SetMarkerColor(kRed);

  me_hEvReco1->SetStats(0);
  me_hEvReco1->SetMarkerStyle(21);
  me_hEvReco1->SetMarkerSize(0.5);
  me_hEvReco1->SetMarkerColor(kBlue);
  ee_hEvReco1->SetStats(0);
  ee_hEvReco1->SetMarkerStyle(21);
  ee_hEvReco1->SetMarkerSize(0.5);
  ee_hEvReco1->SetMarkerColor(kRed);
  em_hEvReco1->SetStats(0);
  em_hEvReco1->SetMarkerStyle(21);
  em_hEvReco1->SetMarkerSize(0.5);
  em_hEvReco1->SetMarkerColor(kBlue);
  mm_hEvReco1->SetStats(0);
  mm_hEvReco1->SetMarkerStyle(21);
  mm_hEvReco1->SetMarkerSize(0.5);
  mm_hEvReco1->SetMarkerColor(kRed);

  me_hEvReco2->SetStats(0);
  me_hEvReco2->SetMarkerStyle(21);
  me_hEvReco2->SetMarkerSize(0.5);
  me_hEvReco2->SetMarkerColor(kBlue);
  ee_hEvReco2->SetStats(0);
  ee_hEvReco2->SetMarkerStyle(21);
  ee_hEvReco2->SetMarkerSize(0.5);
  ee_hEvReco2->SetMarkerColor(kRed);
  em_hEvReco2->SetStats(0);
  em_hEvReco2->SetMarkerStyle(21);
  em_hEvReco2->SetMarkerSize(0.5);
  em_hEvReco2->SetMarkerColor(kBlue);
  mm_hEvReco2->SetStats(0);
  mm_hEvReco2->SetMarkerStyle(21);
  mm_hEvReco2->SetMarkerSize(0.5);
  mm_hEvReco2->SetMarkerColor(kRed);

  THStack *e_hsEv0 = new THStack("e_hsEv0","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 3.0MeV); Ev_reco (GeV); count");
  e_hsEv0->Add(me_hEvReco0);
  e_hsEv0->Add(ee_hEvReco0);
  THStack *m_hsEv0 = new THStack("m_hsEv0","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 3.0MeV); Ev_reco (GeV); count");
  m_hsEv0->Add(em_hEvReco0);
  m_hsEv0->Add(mm_hEvReco0);
  THStack *e_hsEv1 = new THStack("e_hsEv1","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 0.8MeV); Ev_reco (GeV); count");
  e_hsEv1->Add(me_hEvReco1);
  e_hsEv1->Add(ee_hEvReco1);
  THStack *m_hsEv1 = new THStack("m_hsEv1","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 0.8MeV); Ev_reco (GeV); count");
  m_hsEv1->Add(em_hEvReco1);
  m_hsEv1->Add(mm_hEvReco1);
  THStack *e_hsEv2 = new THStack("e_hsEv2","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 0.5MeV); Ev_reco (GeV); count");
  e_hsEv2->Add(me_hEvReco2);
  e_hsEv2->Add(ee_hEvReco2);
  THStack *m_hsEv2 = new THStack("m_hsEv2","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 0.5MeV); Ev_reco (GeV); count");
  m_hsEv2->Add(em_hEvReco2);
  m_hsEv2->Add(mm_hEvReco2);

  TCanvas *cTarget_Ev = new TCanvas("cTarget_Ev","",2700,1400);
  cTarget_Ev->Divide(3,2);
  cTarget_Ev->cd(1);
  e_hEvReco0->Draw();
  e_hsEv0->Draw("same");
  TLegend *leEv0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv0->AddEntry(e_hEvReco0,"true target");
  leEv0->AddEntry(me_hEvReco0,"true oscillated mu->e");
  leEv0->AddEntry(ee_hEvReco0,"true oscillated + unoscillated");
  leEv0->Draw();
  cTarget_Ev->cd(2);
  e_hEvReco1->Draw();
  e_hsEv1->Draw("same");
  TLegend *leEv1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv1->AddEntry(e_hEvReco1,"true target");
  leEv1->AddEntry(me_hEvReco1,"true oscillated mu->e");
  leEv1->AddEntry(ee_hEvReco1,"true oscillated + unoscillated");
  leEv1->Draw();
  cTarget_Ev->cd(3);
  e_hEvReco2->Draw();
  e_hsEv2->Draw("same");
  TLegend *leEv2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  leEv2->AddEntry(e_hEvReco2,"true target");
  leEv2->AddEntry(me_hEvReco2,"true oscillated mu->e");
  leEv2->AddEntry(ee_hEvReco2,"true oscillated + unoscillated");
  leEv2->Draw();
  cTarget_Ev->cd(4);
  m_hEvReco0->Draw();
  m_hsEv0->Draw("same");
  TLegend *lmEv0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv0->AddEntry(m_hEvReco0,"true target");
  lmEv0->AddEntry(em_hEvReco0,"true oscillated mu->e");
  lmEv0->AddEntry(mm_hEvReco0,"true oscillated + unoscillated");
  lmEv0->Draw();
  cTarget_Ev->cd(5);
  m_hEvReco1->Draw();
  m_hsEv1->Draw("same");
  TLegend *lmEv1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv1->AddEntry(m_hEvReco1,"true target");
  lmEv1->AddEntry(em_hEvReco1,"true oscillated mu->e");
  lmEv1->AddEntry(mm_hEvReco1,"true oscillated + unoscillated");
  lmEv1->Draw();
  cTarget_Ev->cd(6);
  m_hEvReco2->Draw();
  m_hsEv2->Draw("same");
  TLegend *lmEv2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lmEv2->AddEntry(m_hEvReco2,"true target");
  lmEv2->AddEntry(em_hEvReco2,"true oscillated mu->e");
  lmEv2->AddEntry(mm_hEvReco2,"true oscillated + unoscillated");
  lmEv2->Draw();
  cTarget_Ev->SaveAs(Form("true_nue_Ev_target_%d.png",para));

  os_hEvReco0->SetStats(0);
  os_hEvReco0->SetMarkerStyle(21);
  os_hEvReco0->SetMarkerSize(0.5);
  os_hEvReco0->SetMarkerColor(kBlue);
  unos_hEvReco0->SetStats(0);
  unos_hEvReco0->SetMarkerStyle(21);
  unos_hEvReco0->SetMarkerSize(0.5);
  unos_hEvReco0->SetMarkerColor(kRed);
  hEvReco0->SetStats(0);
  hEvReco0->SetMarkerStyle(21);
  hEvReco0->SetMarkerSize(0.5);
  hEvReco0->SetMarkerColor(kRed);

  os_hEvReco1->SetStats(0);
  os_hEvReco1->SetMarkerStyle(21);
  os_hEvReco1->SetMarkerSize(0.5);
  os_hEvReco1->SetMarkerColor(kBlue);
  unos_hEvReco1->SetStats(0);
  unos_hEvReco1->SetMarkerStyle(21);
  unos_hEvReco1->SetMarkerSize(0.5);
  unos_hEvReco1->SetMarkerColor(kRed);
  hEvReco1->SetStats(0);
  hEvReco1->SetMarkerStyle(21);
  hEvReco1->SetMarkerSize(0.5);
  hEvReco1->SetMarkerColor(kRed);

  os_hEvReco2->SetStats(0);
  os_hEvReco2->SetMarkerStyle(21);
  os_hEvReco2->SetMarkerSize(0.5);
  os_hEvReco2->SetMarkerColor(kBlue);
  unos_hEvReco2->SetStats(0);
  unos_hEvReco2->SetMarkerStyle(21);
  unos_hEvReco2->SetMarkerSize(0.5);
  unos_hEvReco2->SetMarkerColor(kRed);
  hEvReco2->SetStats(0);
  hEvReco2->SetMarkerStyle(21);
  hEvReco2->SetMarkerSize(0.5);
  hEvReco2->SetMarkerColor(kRed);

  TCanvas *cEv_sum = new TCanvas("cEv_sum","",2700,700);
  cEv_sum->Divide(3,1);
  cEv_sum->cd(1);
  hEvReco0_ft->Draw();
  os_hEvReco0->Draw("same");
  hEvReco0->Draw("same");
  TLegend *lEv0 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lEv0->AddEntry(hEvReco0_ft,"fluctuated target");
  lEv0->AddEntry(os_hEvReco0,"true oscillated mu->e & e->mu");
  lEv0->AddEntry(hEvReco0,"true target");
  lEv0->Draw();
  cEv_sum->cd(2);
  hEvReco1_ft->Draw();
  os_hEvReco1->Draw("same");
  hEvReco1->Draw("same");
  TLegend *lEv1 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lEv1->AddEntry(hEvReco1_ft,"fluctuated target");
  lEv1->AddEntry(os_hEvReco1,"true oscillated mu->e & e->mu");
  lEv1->AddEntry(hEvReco1,"true target");
  lEv1->Draw();
  cEv_sum->cd(3);
  hEvReco2_ft->Draw();
  os_hEvReco2->Draw("same");
  hEvReco2->Draw("same");
  TLegend *lEv2 = new TLegend(0.65, 0.70, 0.9, 0.9);
  lEv2->AddEntry(hEvReco2_ft,"fluctuated target");
  lEv2->AddEntry(os_hEvReco2,"true oscillated mu->e & e->mu");
  lEv2->AddEntry(unos_hEvReco2,"true target");
  lEv2->Draw();
  cEv_sum->SaveAs(Form("true_nue_Ev_target_%d_sum.png",para));
  
  m_hEvRecoVsEv0->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv1->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.8MeV)");
  m_hEvRecoVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv1->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv1->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.8MeV)");
  e_hEvRecoVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv1->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv2->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hEvRecoVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv2->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv2->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hEvRecoVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv2->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv = new TCanvas("cEvVsEv","",2700,1400);
  cEvVsEv->Divide(3,2);
  cEvVsEv->cd(1);
  e_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(2);
  e_hEvRecoVsEv1->Draw("colz");
  cEvVsEv->cd(3);
  e_hEvRecoVsEv2->Draw("colz");
  cEvVsEv->cd(4);
  m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(5);
  m_hEvRecoVsEv1->Draw("colz");
  cEvVsEv->cd(6);
  m_hEvRecoVsEv2->Draw("colz");
  cEvVsEv->SaveAs(Form("true_nue_Ev_templates_%d.png",para));

  m_hEvRecoVsEv0_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv1_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.8MeV)");
  m_hEvRecoVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv1_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv1_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.8MeV)");
  e_hEvRecoVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv1_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv2_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hEvRecoVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv2_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv2_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hEvRecoVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv2_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv_cov = new TCanvas("cEvVsEv_cov","",2700,1400);
  cEvVsEv_cov->Divide(3,2);
  cEvVsEv_cov->cd(1);
  e_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(2);
  e_hEvRecoVsEv1_cov->Draw("colz");
  cEvVsEv_cov->cd(3);
  e_hEvRecoVsEv2_cov->Draw("colz");
  cEvVsEv_cov->cd(4);
  m_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(5);
  m_hEvRecoVsEv1_cov->Draw("colz");
  cEvVsEv_cov->cd(6);
  m_hEvRecoVsEv2_cov->Draw("colz");
  cEvVsEv_cov->SaveAs(Form("cov_nue_Ev_%d.png",para));

  TF1 *cut0 = new TF1("cut0","sqrt(1.0*1E3/x)",0.5,60);
  TF1 *cut1 = new TF1("cut1","sqrt(0.8*1E3/x)",0.5,60);
  TF1 *cut2 = new TF1("cut2","sqrt(0.5*1E3/x)",0.5,60);
  BthetaVsEe->Scale(scalePOT);
  BthetaVsEe->SetTitle("Btheta Vs Ee");
  BthetaVsEe->GetXaxis()->SetTitle("Ee (GeV)");
  BthetaVsEe->GetYaxis()->SetTitle("Btheta (mrad)");
  TCanvas *cBthetaVsEe = new TCanvas("cBthetaVsEe","",2700,700);
  cBthetaVsEe->Divide(3,1);
  cBthetaVsEe->cd(1);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut0->Draw("same");
  cBthetaVsEe->cd(2);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut1->Draw("same");
  cBthetaVsEe->cd(3);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut2->Draw("same");
  cBthetaVsEe->SaveAs(Form("nue_BthetaVsEe_%d.png",para));

  hEvRes0->SetTitle("(Ev_reco - Ev_true)/Ev_true (Etheta2 < 3.0MeV)");
  hEvRes0->GetXaxis()->SetTitle("(Ev_reco - Ev_true)/Ev_true");
  hEvRes1->SetTitle("(Ev_reco - Ev_true)/Ev_true (Etheta2 < 0.8MeV)");
  hEvRes1->GetXaxis()->SetTitle("(Ev_reco - Ev_true)/Ev_true");
  hEvRes2->SetTitle("(Ev_reco - Ev_true)/Ev_true (Etheta2 < 0.5MeV)");
  hEvRes2->GetXaxis()->SetTitle("(Ev_reco - Ev_true)/Ev_true");
  TCanvas *cEvRes = new TCanvas("cEvRes","",1800,400);
  cEvRes->Divide(3,1);
  cEvRes->cd(1);
  hEvRes0->Draw();
  cEvRes->cd(2);
  hEvRes1->Draw();
  cEvRes->cd(3);
  hEvRes2->Draw();
  cEvRes->SaveAs(Form("nue_EvRes_%d.png",para));

  TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/gen_data/nuescattering/nue_output_%d.root",para),"RECREATE");
  e_hElep0_w->Write();
  me_hElep0_w->Write();
  ee_hElep0_w->Write();
  m_hElep0_w->Write();
  em_hElep0_w->Write();
  mm_hElep0_w->Write();
  hElep0_w->Write();
  hElep0_w_ft->Write();
  m_hElepVsEv0->Write();
  m_hElepVsEv0_w->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv0_w->Write();
  m_hElepVsEv0_cov->Write();
  e_hElepVsEv0_cov->Write();

  e_hElep1_w->Write();
  me_hElep1_w->Write();
  ee_hElep1_w->Write();
  m_hElep1_w->Write();
  em_hElep1_w->Write();
  mm_hElep1_w->Write();
  hElep1_w->Write();
  hElep1_w_ft->Write();
  m_hElepVsEv1->Write();
  m_hElepVsEv1_w->Write();
  e_hElepVsEv1->Write();
  e_hElepVsEv1_w->Write();
  m_hElepVsEv1_cov->Write();
  e_hElepVsEv1_cov->Write();

  e_hElep2_w->Write();
  me_hElep2_w->Write();
  ee_hElep2_w->Write();
  m_hElep2_w->Write();
  em_hElep2_w->Write();
  mm_hElep2_w->Write();
  hElep2_w->Write();
  hElep2_w_ft->Write();
  m_hElepVsEv2->Write();
  m_hElepVsEv2_w->Write();
  e_hElepVsEv2->Write();
  e_hElepVsEv2_w->Write();
  m_hElepVsEv2_cov->Write();
  e_hElepVsEv2_cov->Write();

  e_hEvReco0->Write();
  me_hEvReco0->Write();
  ee_hEvReco0->Write();
  m_hEvReco0->Write();
  em_hEvReco0->Write();
  mm_hEvReco0->Write();
  hEvReco0->Write();
  hEvReco0_ft->Write();
  m_hEvRecoVsEv0->Write();
  m_hEvRecoVsEv0_w->Write();
  e_hEvRecoVsEv0->Write();
  e_hEvRecoVsEv0_w->Write();
  m_hEvRecoVsEv0_cov->Write();
  e_hEvRecoVsEv0_cov->Write();
  hEvRes0->Write();

  e_hEvReco1->Write();
  me_hEvReco1->Write();
  ee_hEvReco1->Write();
  m_hEvReco1->Write();
  em_hEvReco1->Write();
  mm_hEvReco1->Write();
  hEvReco1->Write();
  hEvReco1_ft->Write();
  m_hEvRecoVsEv1->Write();
  m_hEvRecoVsEv1_w->Write();
  e_hEvRecoVsEv1->Write();
  e_hEvRecoVsEv1_w->Write();
  m_hEvRecoVsEv1_cov->Write();
  e_hEvRecoVsEv1_cov->Write();
  hEvRes1->Write();

  e_hEvReco2->Write();
  me_hEvReco2->Write();
  ee_hEvReco2->Write();
  m_hEvReco2->Write();
  em_hEvReco2->Write();
  mm_hEvReco2->Write();
  hEvReco2->Write();
  hEvReco2_ft->Write();
  m_hEvRecoVsEv2->Write();
  m_hEvRecoVsEv2_w->Write();
  e_hEvRecoVsEv2->Write();
  e_hEvRecoVsEv2_w->Write();
  m_hEvRecoVsEv2_cov->Write();
  e_hEvRecoVsEv2_cov->Write();
  hEvRes2->Write();

  out->Close();

  return(0);
}

