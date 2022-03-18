#include <string>

void combine()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree_CC = new TChain( "cafTree", "cafTree" );
  TChain * tree_nue = new TChain( "tree", "tree" );

  for(int i = 110; i<111; i++){
  tree_CC->Add( Form("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  std::cout << "File number:" << i << "\n";
  } 

  for(int i = 10; i<11; i++){
  tree_nue->Add( Form("/pnfs/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_0%d.root",i) );
  std::cout << "Nue File number:" << i << "\n";
  } 
  // replace this with the location of your file
  // Note: this step is going to print a bunch of Error failure loading library. These are not actually errors, just ROOT being stupid


  // Get the number of "protons on target", which is essentially the number of neutrinos in this file
/*  TChain * meta = new TChain( "meta", "meta" );
  meta->Add( "/Users/qhuong/Documents/UofR/2021_Summer/data/ND_FHC_FV_00.root" ); // make certain this is the exact same file(s)
*/
  TH1D *mm_hElep_w0 = new TH1D("mm_hElep_w0","",100,0,40);
  //TH1D *mm_hElep_w1 = new TH1D("mm_hElep_w1","",100,0,40);
  //TH1D *mm_hElep_w2 = new TH1D("mm_hElep_w2","",100,0,40);
  TH1D *mm_hElep_w3 = new TH1D("mm_hElep_w3","",100,0,16);
  TH1D *ee_hElep_w0 = new TH1D("ee_hElep_w0","",100,0,40);
  //TH1D *ee_hElep_w1 = new TH1D("ee_hElep_w1","",100,0,40);
  //TH1D *ee_hElep_w2 = new TH1D("ee_hElep_w2","",100,0,40);
  TH1D *ee_hElep_w3 = new TH1D("ee_hElep_w3","",100,0,16);
  TH1D *me_hElep_w0 = new TH1D("me_hElep_w0","",100,0,40);
  //TH1D *me_hElep_w1 = new TH1D("me_hElep_w1","",100,0,40);
  //TH1D *me_hElep_w2 = new TH1D("me_hElep_w2","",100,0,40);
  TH1D *me_hElep_w3 = new TH1D("me_hElep_w3","",100,0,16);
  TH1D *em_hElep_w0 = new TH1D("em_hElep_w0","",100,0,40);
  //TH1D *em_hElep_w1 = new TH1D("em_hElep_w1","",100,0,40);
  //TH1D *em_hElep_w2 = new TH1D("em_hElep_w2","",100,0,40);
  TH1D *em_hElep_w3 = new TH1D("em_hElep_w3","",100,0,16);
  TH1D *e_hElep_w0 = new TH1D("e_hElep_w0","",100,0,40);
  //TH1D *e_hElep_w1 = new TH1D("e_hElep_w1","",100,0,40);
  //TH1D *e_hElep_w2 = new TH1D("e_hElep_w2","",100,0,40);
  TH1D *e_hElep_w3 = new TH1D("e_hElep_w3","",100,0,16);
  TH1D *m_hElep_w0 = new TH1D("m_hElep_w0","",100,0,40);
  //TH1D *m_hElep_w1 = new TH1D("m_hElep_w1","",100,0,40);
  //TH1D *m_hElep_w2 = new TH1D("m_hElep_w2","",100,0,40);
  TH1D *m_hElep_w3 = new TH1D("m_hElep_w3","",100,0,16);

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
  TH2D *nue_m_hElepVsEv0 = new TH2D("nue_m_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *nue_m_hElepVsEv3 = new TH2D("nue_m_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nue_e_hElepVsEv0 = new TH2D("nue_e_hElepVsEv0","",nbinsX0,xEdges0,nbinsY0,yEdges0);
  TH2D *nue_e_hElepVsEv3 = new TH2D("nue_e_hElepVsEv3","",nbinsX,xEdges,nbinsY,yEdges);

  double Uee2 = 0.04;
  double Umm2 = 0.01;
  double s2ee2 = 4*Uee2*(1 - Uee2);
  double s2mm2 = 4*Umm2*(1 - Umm2);
  double s2me2 = 4*Uee2*Umm2;
  double dm2 = 6.0; //eV2
  double L = 500; //m

  double nue_Enu, nue_E[2]; // the energy of the neutrino, in GeV
  int nue_pdg[2];

  tree_nue->SetBranchAddress( "Enu", &nue_Enu );
  tree_nue->SetBranchAddress( "E", &nue_E );
  tree_nue->SetBranchAddress( "pdg", &nue_pdg );

  const int N_nue = tree_nue->GetEntries();
  for( int ii = 0; ii < N_nue; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N_nue, N_nue );
    tree_nue->GetEntry(ii);

    double nue_del = 1.27/1E3*dm2*L/nue_Enu;
    double nue_pmue = s2me2 * pow(sin(nue_del),2);
    double nue_pmumu = 1 - s2mm2 * pow(sin(nue_del),2);
    double nue_pee = 1 - s2ee2 * pow(sin(nue_del),2);

    if(nue_pdg[1] == 14){
      //muon to muon
      mm_hElep_w0->Fill(nue_E[0],nue_pmumu);
      mm_hElep_w3->Fill(nue_E[0],nue_pmumu);

      //muon to electron
      me_hElep_w0->Fill(nue_E[0],nue_pmue);
      me_hElep_w3->Fill(nue_E[0],nue_pmue);

      //nue muon templates
      nue_m_hElepVsEv0->Fill(nue_Enu,nue_E[0]);
      nue_m_hElepVsEv3->Fill(nue_Enu,nue_E[0]);
    }

    if(nue_pdg[1] == 12){
      //electron to electron
      ee_hElep_w0->Fill(nue_E[0],nue_pee);
      ee_hElep_w3->Fill(nue_E[0],nue_pee);

      //electron to muon
      em_hElep_w0->Fill(nue_E[0],nue_pmue);
      em_hElep_w3->Fill(nue_E[0],nue_pmue);

      //nue electron templates
      nue_e_hElepVsEv0->Fill(nue_Enu,nue_E[0]);
      nue_e_hElepVsEv3->Fill(nue_Enu,nue_E[0]);
    }
  }

  me_hElep_w0->Scale(6.);
  me_hElep_w3->Scale(6.);
  em_hElep_w0->Scale(1/6.);
  em_hElep_w3->Scale(1/6.);

  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);


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

  
  tree_CC->SetBranchAddress( "vtx_x", &vtx_x );
  tree_CC->SetBranchAddress( "vtx_y", &vtx_y );
  tree_CC->SetBranchAddress( "vtx_z", &vtx_z );
  tree_CC->SetBranchAddress( "Ev", &Ev );
  tree_CC->SetBranchAddress( "nuPDG", &nuPDG );
  tree_CC->SetBranchAddress( "LepPDG", &LepPDG );
  tree_CC->SetBranchAddress( "LepE", &LepE );
  tree_CC->SetBranchAddress( "LepNuAngle", &LepNuAngle );
  tree_CC->SetBranchAddress( "Ev_reco", &Ev_reco );
  tree_CC->SetBranchAddress( "Elep_reco", &Elep_reco );
  tree_CC->SetBranchAddress( "reco_numu", &reco_numu );
  tree_CC->SetBranchAddress( "reco_nue", &reco_nue );
  tree_CC->SetBranchAddress( "muon_contained", &muon_contained );
  tree_CC->SetBranchAddress( "muon_tracker", &muon_tracker );
  tree_CC->SetBranchAddress( "muon_ecal", &muon_ecal );
  tree_CC->SetBranchAddress( "Ehad_veto", &Ehad_veto );
  tree_CC->SetBranchAddress( "eP", &eP );
  tree_CC->SetBranchAddress( "eN", &eN );
  tree_CC->SetBranchAddress( "ePip", &ePip );
  tree_CC->SetBranchAddress( "ePim", &ePim );
  tree_CC->SetBranchAddress( "ePi0", &ePi0 );
  tree_CC->SetBranchAddress( "eOther", &eOther );
  tree_CC->SetBranchAddress( "LepNuAngle", &LepNuAngle );

  const int N_CC = tree_CC->GetEntries();
  for( int ii = 0; ii < N_CC; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N_CC, N_CC );
    tree_CC->GetEntry(ii);

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
    //if( Ev < 0.2 ) continue;

    double nu = eP + eN + ePip + ePim + ePi0 + eOther;
    double reco_nu = Ev_reco - Elep_reco;
      
    /*double Uee2 = 0.04;
    double Umm2 = 0.01;
    double s2ee2 = 4*Uee2*(1 - Uee2); 
    double s2mm2 = 4*Umm2*(1 - Umm2); 
    double s2me2 = 4*Uee2*Umm2;
    double dm2 = 10.0; //eV2
    double L = 500; //m*/
    double del = 1.27/1E3*dm2*L/Ev;
    double pmue = s2me2 * pow(sin(del),2);
    double pmumu = 1 - s2mm2 * pow(sin(del),2);   
    double pee = 1 - s2ee2 * pow(sin(del),2);
    TRandom *r1 = new TRandom();
    double LepE_sm = r1->Gaus(LepE,0.05);

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
  
  //target: electron = muon-to-electron + electron-to-electron
  e_hElep_w0->Add(me_hElep_w0); e_hElep_w0->Add(ee_hElep_w0);
  e_hElep_w3->Add(me_hElep_w3); e_hElep_w3->Add(ee_hElep_w3);

  //target: muon = electron-to-muon + muon-to-muon
  m_hElep_w0->Add(em_hElep_w0); m_hElep_w0->Add(mm_hElep_w0);
  m_hElep_w3->Add(em_hElep_w3); m_hElep_w3->Add(mm_hElep_w3);
 
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
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
  TCanvas *cElep_w = new TCanvas("cElep_w","",1600,1200);
  cElep_w->Divide(2,2);
  cElep_w->cd(1);
  e_hsElep_w0->Draw();
  cElep_w->cd(2);
  m_hsElep_w0->Draw();
  cElep_w->cd(3);
  e_hsElep_w3->Draw();
  cElep_w->cd(4);
  m_hsElep_w3->Draw();
  cElep_w->SaveAs("Elep_reco_w_4.pdf");

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
  cElepVsEv->SaveAs("CC_ElepVsEv_4.pdf");
  
  nue_m_hElepVsEv0->SetTitle("nu+e numu Elep_reco Vs true Ev");
  nue_m_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nue_m_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  nue_m_hElepVsEv3->SetTitle("nu+e numu Elep_reco Vs true Ev (zoom)");
  nue_m_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nue_m_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  nue_e_hElepVsEv0->SetTitle("nu+e nue Elep_reco Vs true Ev");
  nue_e_hElepVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nue_e_hElepVsEv0->GetYaxis()->SetTitle("Elep_reco");
  nue_e_hElepVsEv3->SetTitle("nu+e nue Elep_reco Vs true Ev (zoom)");
  nue_e_hElepVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nue_e_hElepVsEv3->GetYaxis()->SetTitle("Elep_reco");
  TCanvas *nue_cElepVsEv = new TCanvas("nue_cElepVsEv","",1600,1200);
  nue_cElepVsEv->Divide(2,2);
  nue_cElepVsEv->cd(1);
  nue_m_hElepVsEv0->Draw("colz");
  nue_cElepVsEv->cd(2);
  nue_e_hElepVsEv0->Draw("colz");
  nue_cElepVsEv->cd(3);
  nue_m_hElepVsEv3->Draw("colz");
  nue_cElepVsEv->cd(4);
  nue_e_hElepVsEv3->Draw("colz");
  nue_cElepVsEv->SaveAs("nue_ElepVsEv_4.pdf");

  TFile *out = new TFile("/nashome/q/qvuong/gen_data/output_4.root","RECREATE");
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
  nc_m_hElepVsEv0->Write();
  nc_m_hElepVsEv3->Write();
  e_hElepVsEv0->Write();
  e_hElepVsEv3->Write();
  nue_m_hElepVsEv0->Write();
  nue_m_hElepVsEv3->Write();
  nue_e_hElepVsEv0->Write();
  nue_e_hElepVsEv3->Write();
  out->Close();

}

