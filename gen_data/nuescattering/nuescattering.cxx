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

int main()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "tree", "tree" );
  TChain * meta = new TChain( "meta", "meta" );

  for(int i = 10; i<12; i++){
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

  TH2D *m_hEvRecoVsEv0 = new TH2D("m_hEvRecoVsEv0","",nbinsX,xEdges,100,0,120);  //Etheta2 < 3MeV
  TH2D *m_hEvRecoVsEv1 = new TH2D("m_hEvRecoVsEv1","",nbinsX,xEdges,100,0,120);  //Etheta2 < 1MeV
  TH2D *m_hEvRecoVsEv2 = new TH2D("m_hEvRecoVsEv2","",nbinsX,xEdges,100,0,120);  //Etheta2 < 0.5MeV
  TH2D *e_hEvRecoVsEv0 = new TH2D("e_hEvRecoVsEv0","",nbinsX,xEdges,100,0,120);
  TH2D *e_hEvRecoVsEv1 = new TH2D("e_hEvRecoVsEv1","",nbinsX,xEdges,100,0,120);
  TH2D *e_hEvRecoVsEv2 = new TH2D("e_hEvRecoVsEv2","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv0   = new TH2D("hEvRecoVsEv0","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv1   = new TH2D("hEvRecoVsEv1","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv2   = new TH2D("hEvRecoVsEv2","",nbinsX,xEdges,100,0,120);
  TH2D *m_hEvRecoVsEv0_w = new TH2D("m_hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,120);
  TH2D *m_hEvRecoVsEv1_w = new TH2D("m_hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,120);
  TH2D *m_hEvRecoVsEv2_w = new TH2D("m_hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,120);
  TH2D *e_hEvRecoVsEv0_w = new TH2D("e_hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,120);
  TH2D *e_hEvRecoVsEv1_w = new TH2D("e_hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,120);
  TH2D *e_hEvRecoVsEv2_w = new TH2D("e_hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv0_w   = new TH2D("hEvRecoVsEv0_w","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv1_w   = new TH2D("hEvRecoVsEv1_w","",nbinsX,xEdges,100,0,120);
  TH2D *hEvRecoVsEv2_w   = new TH2D("hEvRecoVsEv2_w","",nbinsX,xEdges,100,0,120);

  TH1D *hEvRes0 = new TH1D("hEvRes0","",100,-1,1);
  TH1D *hEvRes1 = new TH1D("hEvRes1","",100,-1,1);
  TH1D *hEvRes2 = new TH1D("hEvRes2","",100,-1,1);

  TH1D *m_hEvReco0 = new TH1D("m_hEvReco0","",100,0,120);
  TH1D *m_hEvReco1 = new TH1D("m_hEvReco1","",100,0,120);
  TH1D *m_hEvReco2 = new TH1D("m_hEvReco2","",100,0,120);
  TH1D *mm_hEvReco0 = new TH1D("mm_hEvReco0","",100,0,120);
  TH1D *mm_hEvReco1 = new TH1D("mm_hEvReco1","",100,0,120);
  TH1D *mm_hEvReco2 = new TH1D("mm_hEvReco2","",100,0,120);
  TH1D *me_hEvReco0 = new TH1D("me_hEvReco0","",100,0,120);
  TH1D *me_hEvReco1 = new TH1D("me_hEvReco1","",100,0,120);
  TH1D *me_hEvReco2 = new TH1D("me_hEvReco2","",100,0,120);
  TH1D *e_hEvReco0 = new TH1D("e_hEvReco0","",100,0,120);
  TH1D *e_hEvReco1 = new TH1D("e_hEvReco1","",100,0,120);
  TH1D *e_hEvReco2 = new TH1D("e_hEvReco2","",100,0,120);
  TH1D *ee_hEvReco0 = new TH1D("ee_hEvReco0","",100,0,120);
  TH1D *ee_hEvReco1 = new TH1D("ee_hEvReco1","",100,0,120);
  TH1D *ee_hEvReco2 = new TH1D("ee_hEvReco2","",100,0,120);
  TH1D *em_hEvReco0 = new TH1D("em_hEvReco0","",100,0,120);
  TH1D *em_hEvReco1 = new TH1D("em_hEvReco1","",100,0,120);
  TH1D *em_hEvReco2 = new TH1D("em_hEvReco2","",100,0,120);
  TH1D *hEvReco0 = new TH1D("hEvReco0","",100,0,120);
  TH1D *hEvReco1 = new TH1D("hEvReco1","",100,0,120);
  TH1D *hEvReco2 = new TH1D("hEvReco2","",100,0,120);
  TH1D *os_hEvReco0 = new TH1D("os_hEvReco0","",100,0,120);
  TH1D *os_hEvReco1 = new TH1D("os_hEvReco1","",100,0,120);
  TH1D *os_hEvReco2 = new TH1D("os_hEvReco2","",100,0,120);
  TH1D *unos_hEvReco0 = new TH1D("unos_hEvReco0","",100,0,120);
  TH1D *unos_hEvReco1 = new TH1D("unos_hEvReco1","",100,0,120);
  TH1D *unos_hEvReco2 = new TH1D("unos_hEvReco2","",100,0,120);
  TH1D *hEvReco0_ft = new TH1D("hEvReco0_ft","",100,0,120);
  TH1D *hEvReco1_ft = new TH1D("hEvReco1_ft","",100,0,120);
  TH1D *hEvReco2_ft = new TH1D("hEvReco2_ft","",100,0,120);

  TH2D *m_hEvRecoVsEv0_cov = new TH2D("m_hEvRecoVsEv0_cov","",19,mubins,100,0,120);
  TH2D *m_hEvRecoVsEv1_cov = new TH2D("m_hEvRecoVsEv1_cov","",19,mubins,100,0,120);
  TH2D *m_hEvRecoVsEv2_cov = new TH2D("m_hEvRecoVsEv2_cov","",19,mubins,100,0,120);
  TH2D *e_hEvRecoVsEv0_cov = new TH2D("e_hEvRecoVsEv0_cov","",7,ebins,100,0,120);
  TH2D *e_hEvRecoVsEv1_cov = new TH2D("e_hEvRecoVsEv1_cov","",7,ebins,100,0,120);
  TH2D *e_hEvRecoVsEv2_cov = new TH2D("e_hEvRecoVsEv2_cov","",7,ebins,100,0,120);

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

  double Uee2 = 0.04;
  double Umm2 = 0.01;
  double s2ee2 = 4*Uee2*(1 - Uee2); 
  double s2mm2 = 4*Umm2*(1 - Umm2); 
  double s2me2 = 4*Uee2*Umm2;
  double dm2 = 6.0; //eV2
  double L = 500; //m
  double del, pmumu, pmue, pee;

  double me = 510;  //keV
  double Ev_reco;

  const int N = tree->GetEntries();
  //const int N = 10000;
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
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
  
    if(E[0]*Btheta_sm*Btheta_sm/1E3<3){
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

    if(E[0]*Btheta_sm*Btheta_sm/1E3<1){
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

  e_hEvReco0->Add(me_hEvReco0);    e_hEvReco0->Add(ee_hEvReco0);
  m_hEvReco0->Add(em_hEvReco0);    m_hEvReco0->Add(mm_hEvReco0);
  hEvReco0->Add(e_hEvReco0);       hEvReco0->Add(m_hEvReco0);
  os_hEvReco0->Add(me_hEvReco0);   os_hEvReco0->Add(em_hEvReco0);
  unos_hEvReco0->Add(ee_hEvReco0); unos_hEvReco0->Add(mm_hEvReco0);

  e_hElep1_w->Add(me_hElep1_w);    e_hElep1_w->Add(ee_hElep1_w);
  m_hElep1_w->Add(em_hElep1_w);    m_hElep1_w->Add(mm_hElep1_w);
  hElep1_w->Add(e_hElep1_w);       hElep1_w->Add(m_hElep1_w);  
  os_hElep1_w->Add(me_hElep1_w);   os_hElep1_w->Add(em_hElep1_w);
  unos_hElep1_w->Add(ee_hElep1_w); unos_hElep1_w->Add(mm_hElep1_w);

  e_hEvReco1->Add(me_hEvReco1);    e_hEvReco1->Add(ee_hEvReco1);
  m_hEvReco1->Add(em_hEvReco1);    m_hEvReco1->Add(mm_hEvReco1);
  hEvReco1->Add(e_hEvReco1);       hEvReco1->Add(m_hEvReco1);
  os_hEvReco1->Add(me_hEvReco1);   os_hEvReco1->Add(em_hEvReco1);
  unos_hEvReco1->Add(ee_hEvReco1); unos_hEvReco1->Add(mm_hEvReco1);

  e_hElep2_w->Add(me_hElep2_w);    e_hElep2_w->Add(ee_hElep2_w);
  m_hElep2_w->Add(em_hElep2_w);    m_hElep2_w->Add(mm_hElep2_w);
  hElep2_w->Add(e_hElep2_w);       hElep2_w->Add(m_hElep2_w);  
  os_hElep2_w->Add(me_hElep2_w);   os_hElep2_w->Add(em_hElep2_w);
  unos_hElep2_w->Add(ee_hElep2_w); unos_hElep2_w->Add(mm_hElep2_w);

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

  hElep0_w_ft->SetStats(0);
  hElep0_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 3MeV)");
  hElep0_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hElep1_w_ft->SetStats(0);
  hElep1_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 1MeV)");
  hElep1_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hElep2_w_ft->SetStats(0);
  hElep2_w_ft->SetTitle("fluctuated nu+e Elep target (Etheta2 < 0.5MeV)");
  hElep2_w_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cTarget_Elep = new TCanvas("cTarget_Elep","",1800,400);
  cTarget_Elep->Divide(3,1);
  cTarget_Elep->cd(1);
  hElep0_w_ft->Draw();
  cTarget_Elep->cd(2);
  hElep1_w_ft->Draw();
  cTarget_Elep->cd(3);
  hElep2_w_ft->Draw();
  cTarget_Elep->SaveAs("true_nue_ElepReco_target_ft_2.png");

  hEvReco0_ft->SetStats(0);
  hEvReco0_ft->SetTitle("fluctuated nu+e Ev_reco target  (Etheta2 < 3MeV)");
  hEvReco0_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hEvReco1_ft->SetStats(0);
  hEvReco1_ft->SetTitle("fluctuated nu+e Ev_reco target (Etheta2 < 1MeV)");
  hEvReco1_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hEvReco2_ft->SetStats(0);
  hEvReco2_ft->SetTitle("fluctuated nu+e Ev_reco target (Etheta2 < 0.5MeV)");
  hEvReco2_ft->GetXaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cTarget_Ev = new TCanvas("cTarget_Ev","",1800,400);
  cTarget_Ev->Divide(3,1);
  cTarget_Ev->cd(1);
  gPad->SetLogx();
  hEvReco0_ft->Draw();
  cTarget_Ev->cd(2);
  gPad->SetLogx();
  hEvReco1_ft->Draw();
  cTarget_Ev->cd(3);
  gPad->SetLogx();
  hEvReco2_ft->Draw();
  cTarget_Ev->SaveAs("true_nue_EvReco_target_ft_2.png");

  THStack *e_hsElep0_w = new THStack("e_hsElep0_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 3MeV); Elep_reco (GeV); count");
  me_hElep0_w->SetMarkerStyle(21);
  me_hElep0_w->SetMarkerSize(0.5);
  me_hElep0_w->SetMarkerColor(kRed);
  me_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep0_w->Add(me_hElep0_w);
  ee_hElep0_w->SetMarkerStyle(21);
  ee_hElep0_w->SetMarkerSize(0.5);
  ee_hElep0_w->SetMarkerColor(kBlue);
  ee_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep0_w->Add(ee_hElep0_w);
  THStack *m_hsElep0_w = new THStack("m_hsElep0_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 3MeV); Elep_reco (GeV); count");
  em_hElep0_w->SetMarkerStyle(21);
  em_hElep0_w->SetMarkerSize(0.5);
  em_hElep0_w->SetMarkerColor(kRed);
  em_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep0_w->Add(em_hElep0_w);
  mm_hElep0_w->SetMarkerStyle(21);
  mm_hElep0_w->SetMarkerSize(0.5);
  mm_hElep0_w->SetMarkerColor(kBlue);
  mm_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep0_w->Add(mm_hElep0_w);
  THStack *e_hsElep1_w = new THStack("e_hsElep1_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 1MeV); Elep_reco (GeV); count");
  me_hElep1_w->SetMarkerStyle(21);
  me_hElep1_w->SetMarkerSize(0.5);
  me_hElep1_w->SetMarkerColor(kRed);
  me_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep1_w->Add(me_hElep1_w);
  ee_hElep1_w->SetMarkerStyle(21);
  ee_hElep1_w->SetMarkerSize(0.5);
  ee_hElep1_w->SetMarkerColor(kBlue);
  ee_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep1_w->Add(ee_hElep1_w);
  THStack *m_hsElep1_w = new THStack("m_hsElep1_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 1MeV); Elep_reco (GeV); count");
  em_hElep1_w->SetMarkerStyle(21);
  em_hElep1_w->SetMarkerSize(0.5);
  em_hElep1_w->SetMarkerColor(kRed);
  em_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep1_w->Add(em_hElep1_w);
  mm_hElep1_w->SetMarkerStyle(21);
  mm_hElep1_w->SetMarkerSize(0.5);
  mm_hElep1_w->SetMarkerColor(kBlue);
  mm_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep1_w->Add(mm_hElep1_w);
  THStack *e_hsElep2_w = new THStack("e_hsElep2_w","nu+e oscillated + unoscillated nue weighted Elep_reco (Etheta2 < 0.5MeV); Elep_reco (GeV); count");
  me_hElep2_w->SetMarkerStyle(21);
  me_hElep2_w->SetMarkerSize(0.5);
  me_hElep2_w->SetMarkerColor(kRed);
  me_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep2_w->Add(me_hElep2_w);
  ee_hElep2_w->SetMarkerStyle(21);
  ee_hElep2_w->SetMarkerSize(0.5);
  ee_hElep2_w->SetMarkerColor(kBlue);
  ee_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hsElep2_w->Add(ee_hElep2_w);
  THStack *m_hsElep2_w = new THStack("m_hsElep2_w","nu+e oscillated + unoscillated numu weighted Elep_reco (Etheta2 < 0.5MeV); Elep_reco (GeV); count");
  em_hElep2_w->SetMarkerStyle(21);
  em_hElep2_w->SetMarkerSize(0.5);
  em_hElep2_w->SetMarkerColor(kRed);
  em_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep2_w->Add(em_hElep2_w);
  mm_hElep2_w->SetMarkerStyle(21);
  mm_hElep2_w->SetMarkerSize(0.5);
  mm_hElep2_w->SetMarkerColor(kBlue);
  mm_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hsElep2_w->Add(mm_hElep2_w);
  TCanvas *cElep_w = new TCanvas("cElep_w","",1200,1200);
  cElep_w->Divide(2,3);
  cElep_w->cd(1);
  m_hsElep0_w->Draw();
  cElep_w->cd(2);
  e_hsElep0_w->Draw();
  cElep_w->cd(3);
  m_hsElep1_w->Draw();
  cElep_w->cd(4);
  e_hsElep1_w->Draw();
  cElep_w->cd(5);
  m_hsElep2_w->Draw();
  cElep_w->cd(6);
  e_hsElep2_w->Draw();
  cElep_w->SaveAs("true_nue_ElepReco_w_2.png");

  THStack *hsElep0_w = new THStack("hsElep0_w","oscillated + unoscillated weighted Elep_reco (Etheta2 < 3MeV); Elep_reco (GeV); count");
  os_hElep0_w->SetMarkerStyle(21);
  os_hElep0_w->SetMarkerSize(0.5);
  os_hElep0_w->SetMarkerColor(kRed);
  os_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep0_w->Add(os_hElep0_w);
  unos_hElep0_w->SetMarkerStyle(21);
  unos_hElep0_w->SetMarkerSize(0.5);
  unos_hElep0_w->SetMarkerColor(kBlue);
  unos_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep0_w->Add(unos_hElep0_w);
  THStack *hsElep1_w = new THStack("hsElep1_w","oscillated + unoscillated weighted Elep_reco (Etheta2 < 1MeV); Elep_reco (GeV); count");
  os_hElep1_w->SetMarkerStyle(21);
  os_hElep1_w->SetMarkerSize(0.5);
  os_hElep1_w->SetMarkerColor(kRed);
  os_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep1_w->Add(os_hElep1_w);
  unos_hElep1_w->SetMarkerStyle(21);
  unos_hElep1_w->SetMarkerSize(0.5);
  unos_hElep1_w->SetMarkerColor(kBlue);
  unos_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep1_w->Add(unos_hElep1_w);
  THStack *hsElep2_w = new THStack("hsElep2_w","oscillated + unoscillated weighted Elep_reco (Etheta2 < 0.5MeV); Elep_reco (GeV); count");
  os_hElep2_w->SetMarkerStyle(21);
  os_hElep2_w->SetMarkerSize(0.5);
  os_hElep2_w->SetMarkerColor(kRed);
  os_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep2_w->Add(os_hElep2_w);
  unos_hElep2_w->SetMarkerStyle(21);
  unos_hElep2_w->SetMarkerSize(0.5);
  unos_hElep2_w->SetMarkerColor(kBlue);
  unos_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  hsElep2_w->Add(unos_hElep2_w);
  TCanvas *cElep_w_sum = new TCanvas("cElep_w_sum","",1500,400);
  cElep_w_sum->Divide(3,1);
  cElep_w_sum->cd(1);
  hsElep0_w->Draw();
  cElep_w_sum->cd(2);
  hsElep1_w->Draw();
  cElep_w_sum->cd(3);
  hsElep2_w->Draw();
  cElep_w_sum->SaveAs("true_nue_ElepReco_w_2_sum.png");

  THStack *e_hsEv0 = new THStack("e_hsEv0","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 3MeV); Ev_reco (GeV); count");
  me_hEvReco0->SetMarkerStyle(21);
  me_hEvReco0->SetMarkerSize(0.5);
  me_hEvReco0->SetMarkerColor(kRed);
  me_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv0->Add(me_hEvReco0);
  ee_hEvReco0->SetMarkerStyle(21);
  ee_hEvReco0->SetMarkerSize(0.5);
  ee_hEvReco0->SetMarkerColor(kBlue);
  ee_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv0->Add(ee_hEvReco0);
  THStack *m_hsEv0 = new THStack("m_hsEv0","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 3MeV); Ev_reco (GeV); count");
  em_hEvReco0->SetMarkerStyle(21);
  em_hEvReco0->SetMarkerSize(0.5);
  em_hEvReco0->SetMarkerColor(kRed);
  em_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv0->Add(em_hEvReco0);
  mm_hEvReco0->SetMarkerStyle(21);
  mm_hEvReco0->SetMarkerSize(0.5);
  mm_hEvReco0->SetMarkerColor(kBlue);
  mm_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv0->Add(mm_hEvReco0);
  THStack *e_hsEv1 = new THStack("e_hsEv1","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 1MeV); Ev_reco (GeV); count");
  me_hEvReco1->SetMarkerStyle(21);
  me_hEvReco1->SetMarkerSize(0.5);
  me_hEvReco1->SetMarkerColor(kRed);
  me_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv1->Add(me_hEvReco1);
  ee_hEvReco1->SetMarkerStyle(21);
  ee_hEvReco1->SetMarkerSize(0.5);
  ee_hEvReco1->SetMarkerColor(kBlue);
  ee_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv1->Add(ee_hEvReco1);
  THStack *m_hsEv1 = new THStack("m_hsEv1","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 1MeV); Ev_reco (GeV); count");
  em_hEvReco1->SetMarkerStyle(21);
  em_hEvReco1->SetMarkerSize(0.5);
  em_hEvReco1->SetMarkerColor(kRed);
  em_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv1->Add(em_hEvReco1);
  mm_hEvReco1->SetMarkerStyle(21);
  mm_hEvReco1->SetMarkerSize(0.5);
  mm_hEvReco1->SetMarkerColor(kBlue);
  mm_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv1->Add(mm_hEvReco1);
  THStack *e_hsEv2 = new THStack("e_hsEv2","nu+e oscillated + unoscillated nue weighted Ev_reco (Etheta2 < 0.5MeV); Ev_reco (GeV); count");
  me_hEvReco2->SetMarkerStyle(21);
  me_hEvReco2->SetMarkerSize(0.5);
  me_hEvReco2->SetMarkerColor(kRed);
  me_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv2->Add(me_hEvReco2);
  ee_hEvReco2->SetMarkerStyle(21);
  ee_hEvReco2->SetMarkerSize(0.5);
  ee_hEvReco2->SetMarkerColor(kBlue);
  ee_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hsEv2->Add(ee_hEvReco2);
  THStack *m_hsEv2 = new THStack("m_hsEv2","nu+e oscillated + unoscillated numu weighted Ev_reco (Etheta2 < 0.5MeV); Ev_reco (GeV); count");
  em_hEvReco2->SetMarkerStyle(21);
  em_hEvReco2->SetMarkerSize(0.5);
  em_hEvReco2->SetMarkerColor(kRed);
  em_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv2->Add(em_hEvReco2);
  mm_hEvReco2->SetMarkerStyle(21);
  mm_hEvReco2->SetMarkerSize(0.5);
  mm_hEvReco2->SetMarkerColor(kBlue);
  mm_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hsEv2->Add(mm_hEvReco2);
  TCanvas *cEv = new TCanvas("cEv","",1200,1200);
  cEv->Divide(2,3);
  cEv->cd(1);
  gPad->SetLogx();
  m_hsEv0->Draw();
  cEv->cd(2);
  gPad->SetLogx();
  e_hsEv0->Draw();
  cEv->cd(3);
  gPad->SetLogx();
  m_hsEv1->Draw();
  cEv->cd(4);
  gPad->SetLogx();
  e_hsEv1->Draw();
  cEv->cd(5);
  gPad->SetLogx();
  m_hsEv2->Draw();
  cEv->cd(6);
  gPad->SetLogx();
  e_hsEv2->Draw();
  cEv->SaveAs("true_nue_EvReco_w_2.png");

  THStack *hsEv0 = new THStack("hsEv0","oscillated + unoscillated weighted Ev_reco (Etheta2 < 3MeV); Ev_reco (GeV); count");
  os_hEvReco0->SetMarkerStyle(21);
  os_hEvReco0->SetMarkerSize(0.5);
  os_hEvReco0->SetMarkerColor(kRed);
  os_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv0->Add(os_hEvReco0);
  unos_hEvReco0->SetMarkerStyle(21);
  unos_hEvReco0->SetMarkerSize(0.5);
  unos_hEvReco0->SetMarkerColor(kBlue);
  unos_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv0->Add(unos_hEvReco0);
  THStack *hsEv1 = new THStack("hsEv1","oscillated + unoscillated weighted Ev_reco (Etheta2 < 1MeV); Ev_reco (GeV); count");
  os_hEvReco1->SetMarkerStyle(21);
  os_hEvReco1->SetMarkerSize(0.5);
  os_hEvReco1->SetMarkerColor(kRed);
  os_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv1->Add(os_hEvReco1);
  unos_hEvReco1->SetMarkerStyle(21);
  unos_hEvReco1->SetMarkerSize(0.5);
  unos_hEvReco1->SetMarkerColor(kBlue);
  unos_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv1->Add(unos_hEvReco1);
  THStack *hsEv2 = new THStack("hsEv2","oscillated + unoscillated weighted Ev_reco (Etheta2 < 0.5MeV); Ev_reco (GeV); count");
  os_hEvReco2->SetMarkerStyle(21);
  os_hEvReco2->SetMarkerSize(0.5);
  os_hEvReco2->SetMarkerColor(kRed);
  os_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv2->Add(os_hEvReco2);
  unos_hEvReco2->SetMarkerStyle(21);
  unos_hEvReco2->SetMarkerSize(0.5);
  unos_hEvReco2->SetMarkerColor(kBlue);
  unos_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  hsEv2->Add(unos_hEvReco2);
  TCanvas *cEv_sum = new TCanvas("cEv_sum","",1800,400);
  cEv_sum->Divide(3,1);
  cEv_sum->cd(1);
  gPad->SetLogx();
  hsEv0->Draw();
  cEv_sum->cd(2);
  gPad->SetLogx();
  hsEv1->Draw();
  cEv_sum->cd(3);
  gPad->SetLogx();
  hsEv2->Draw();
  cEv_sum->SaveAs("true_nue_EvReco_w_2_sum.png");
  
  e_hElep0_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 3MeV)");
  e_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 3MeV)");
  m_hElep0_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep1_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 1MeV)");
  e_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep1_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 1MeV)");
  m_hElep1_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep2_w->SetTitle("nu+e nue weighted Elep_reco (Etheta2 < 0.5MeV)");
  e_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep2_w->SetTitle("nu+e numu weighted Elep_reco (Etheta2 < 0.5MeV)");
  m_hElep2_w->GetXaxis()->SetTitle("Elep_reco (GeV)");

  e_hEvReco0->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 3MeV)");
  e_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco0->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 3MeV)");
  m_hEvReco0->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco1->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 1MeV)");
  e_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco1->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 1MeV)");
  m_hEvReco1->GetXaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvReco2->SetTitle("nu+e nue weighted Ev_reco (Etheta2 < 0.5MeV)");
  e_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvReco2->SetTitle("nu+e numu weighted Ev_reco (Etheta2 < 0.5MeV)");
  m_hEvReco2->GetXaxis()->SetTitle("Ev_reco (GeV)");

  m_hElepVsEv0->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3MeV)");
  m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3MeV)");
  e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv1->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 1MeV)");
  m_hElepVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv1->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv1->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 1MeV)");
  e_hElepVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv1->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv2->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hElepVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv2->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hElepVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",1200,1200);
  cElepVsEv->Divide(2,3);
  cElepVsEv->cd(1);
  m_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  e_hElepVsEv0->Draw("colz");
  cElepVsEv->cd(3);
  m_hElepVsEv1->Draw("colz");
  cElepVsEv->cd(4);
  e_hElepVsEv1->Draw("colz");
  cElepVsEv->cd(5);
  m_hElepVsEv2->Draw("colz");
  cElepVsEv->cd(6);
  e_hElepVsEv2->Draw("colz");
  cElepVsEv->SaveAs("true_nue_ElepRecoVsEv_2.png");

  m_hEvRecoVsEv0->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 3MeV)");
  m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 3MeV)");
  e_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv1->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 1MeV)");
  m_hEvRecoVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv1->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv1->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 1MeV)");
  e_hEvRecoVsEv1->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv1->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv2->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hEvRecoVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv2->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv2->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hEvRecoVsEv2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv2->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv = new TCanvas("cEvVsEv","",1200,1200);
  cEvVsEv->Divide(2,3);
  cEvVsEv->cd(1);
  m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(2);
  e_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(3);
  m_hEvRecoVsEv1->Draw("colz");
  cEvVsEv->cd(4);
  e_hEvRecoVsEv1->Draw("colz");
  cEvVsEv->cd(5);
  m_hEvRecoVsEv2->Draw("colz");
  cEvVsEv->cd(6);
  e_hEvRecoVsEv2->Draw("colz");
  cEvVsEv->SaveAs("true_nue_EvRecoVsEv_2.png");

  m_hElepVsEv0_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3MeV)");
  m_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv0_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3MeV)");
  e_hElepVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv1_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 1MeV)");
  m_hElepVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv1_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv1_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 1MeV)");
  e_hElepVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv1_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepVsEv2_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hElepVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepVsEv2_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepVsEv2_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hElepVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepVsEv2_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov = new TCanvas("cElepVsEv_cov","",1200,1200);
  cElepVsEv_cov->Divide(2,3);
  cElepVsEv_cov->cd(1);
  m_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(2);
  e_hElepVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(3);
  m_hElepVsEv1_cov->Draw("colz");
  cElepVsEv_cov->cd(4);
  e_hElepVsEv1_cov->Draw("colz");
  cElepVsEv_cov->cd(5);
  m_hElepVsEv2_cov->Draw("colz");
  cElepVsEv_cov->cd(6);
  e_hElepVsEv2_cov->Draw("colz");
  cElepVsEv_cov->SaveAs("cov_nue_ElepRecoVsEv_2.png");

  m_hEvRecoVsEv0_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 3MeV)");
  m_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 3MeV)");
  e_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv1_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 1MeV)");
  m_hEvRecoVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv1_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv1_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 1MeV)");
  e_hEvRecoVsEv1_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv1_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv2_cov->SetTitle("nu+e numu Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  m_hEvRecoVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv2_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv2_cov->SetTitle("nu+e nue Ev_reco Vs true Ev (Etheta2 < 0.5MeV)");
  e_hEvRecoVsEv2_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv2_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv_cov = new TCanvas("cEvVsEv_cov","",1200,1200);
  cEvVsEv_cov->Divide(2,3);
  cEvVsEv_cov->cd(1);
  m_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(2);
  e_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(3);
  m_hEvRecoVsEv1_cov->Draw("colz");
  cEvVsEv_cov->cd(4);
  e_hEvRecoVsEv1_cov->Draw("colz");
  cEvVsEv_cov->cd(5);
  m_hEvRecoVsEv2_cov->Draw("colz");
  cEvVsEv_cov->cd(6);
  e_hEvRecoVsEv2_cov->Draw("colz");
  cEvVsEv_cov->SaveAs("cov_nue_EvRecoVsEv_2.png");

  TF1 *cut0 = new TF1("cut0","sqrt(3*1E3/x)",0.5,60);
  TF1 *cut1 = new TF1("cut1","sqrt(1*1E3/x)",0.5,60);
  TF1 *cut2 = new TF1("cut2","sqrt(0.5*1E3/x)",0.5,60);
  BthetaVsEe->Scale(scalePOT);
  BthetaVsEe->SetTitle("Btheta Vs Ee");
  BthetaVsEe->GetXaxis()->SetTitle("Ee (GeV)");
  BthetaVsEe->GetYaxis()->SetTitle("Btheta (mrad)");
  TCanvas *cBthetaVsEe = new TCanvas("cBthetaVsEe","",1200,800);
  cBthetaVsEe->Divide(2,2);
  cBthetaVsEe->cd(1);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cBthetaVsEe->cd(2);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut0->Draw("same");
  cBthetaVsEe->cd(3);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut1->Draw("same");
  cBthetaVsEe->cd(4);
  cBthetaVsEe->SetLogx();
  BthetaVsEe->Draw("colz");
  cut2->Draw("same");
  cBthetaVsEe->SaveAs("nue_BthetaVsEe_2.png");

  hEvRes0->SetTitle("(Ev_reco - Ev_true)/Ev_true (Etheta2 < 3MeV)");
  hEvRes0->GetXaxis()->SetTitle("(Ev_reco - Ev_true)/Ev_true");
  hEvRes1->SetTitle("(Ev_reco - Ev_true)/Ev_true (Etheta2 < 1MeV)");
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
  cEvRes->SaveAs("EvRes_2.png");

  TFile *out = new TFile("nue_output_2.root","RECREATE");
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

