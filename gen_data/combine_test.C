#include <string>

void lownu()
{
  TFile f("output_1.root","RECREATE");
  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "tree", "tree" );

  tree->Add( "nueFHC_000.root" );
  
  TH1D *nue_mm_hElep_w0 = new TH1D("nue_mm_hElep_w0","",100,0,40);
  //TH1D *mm_hElep_w1 = new TH1D("mm_hElep_w1","",100,0,40);
  //TH1D *mm_hElep_w2 = new TH1D("mm_hElep_w2","",100,0,40);
  TH1D *nue_mm_hElep_w3 = new TH1D("nue_mm_hElep_w3","",100,0,16);
  TH1D *nue_ee_hElep_w0 = new TH1D("nue_ee_hElep_w0","",100,0,40);
  //TH1D *ee_hElep_w1 = new TH1D("ee_hElep_w1","",100,0,40);
  //TH1D *ee_hElep_w2 = new TH1D("ee_hElep_w2","",100,0,40);
  TH1D *nue_ee_hElep_w3 = new TH1D("nue_ee_hElep_w3","",100,0,16);
  TH1D *nue_me_hElep_w0 = new TH1D("nue_me_hElep_w0","",100,0,40);
  //TH1D *me_hElep_w1 = new TH1D("me_hElep_w1","",100,0,40);
  //TH1D *me_hElep_w2 = new TH1D("me_hElep_w2","",100,0,40);
  TH1D *nue_me_hElep_w3 = new TH1D("nue_me_hElep_w3","",100,0,16);
  TH1D *nue_em_hElep_w0 = new TH1D("nue_em_hElep_w0","",100,0,40);
  //TH1D *me_hElep_w1 = new TH1D("me_hElep_w1","",100,0,40);
  //TH1D *me_hElep_w2 = new TH1D("me_hElep_w2","",100,0,40);
  TH1D *nue_em_hElep_w3 = new TH1D("nue_em_hElep_w3","",100,0,16);
  TH1D *nue_e_hElep_w0 = new TH1D("nue_e_hElep_w0","",100,0,40);
  //TH1D *e_hElep_w1 = new TH1D("e_hElep_w1","",100,0,40);
  //TH1D *e_hElep_w2 = new TH1D("e_hElep_w2","",100,0,40);
  TH1D *nue_e_hElep_w3 = new TH1D("nue_e_hElep_w3","",100,0,16);
  TH1D *nue_m_hElep_w0 = new TH1D("nue_m_hElep_w0","",100,0,40);
  //TH1D *m_hElep_w1 = new TH1D("m_hElep_w1","",100,0,40);
  //TH1D *m_hElep_w2 = new TH1D("m_hElep_w2","",100,0,40);
  TH1D *nue_m_hElep_w3 = new TH1D("nue_m_hElep_w3","",100,0,16);
/*
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
*/
  // Most of them are weights related to systematic uncertainties
  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  double Enu, E[2]; // the energy of the neutrino, in GeV
  double perf_px[2], perf_py[2], perf_pz[2];
  double best_px[2], best_py[2], best_pz[2];
  int pdg[2];

  tree->SetBranchAddress( "Enu", &Enu );
  tree->SetBranchAddress( "E", &E );
  tree->SetBranchAddress( "pdg", &pdg );

  const int N = tree->GetEntries();
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    double Uee2 = 0.04;
    double Umm2 = 0.01;
    double s2ee2 = 4*Uee2*(1 - Uee2);
    double s2mm2 = 4*Umm2*(1 - Umm2);
    double s2me2 = 4*Uee2*Umm2;
    double dm2 = 6.0; //eV2
    double L = 500; //m
    double del = 1.27/1E3*dm2*L/Enu;
    double pmue = s2me2 * pow(sin(del),2);
    double pmumu = 1 - s2mm2 * pow(sin(del),2);
    double pee = 1 - s2ee2 * pow(sin(del),2);

    if(pdg[1] == 14){
      nue_mm_hElep_w0->Fill(E[0],pmumu);
      nue_mm_hElep_w3->Fill(E[0],pmumu);
      
      nue_me_hElep_w0->Fill(E[0],pmue);
      nue_me_hElep_w3->Fill(E[0],pmue);
      
      nue_m_hElepVsEv0->Fill(Enu,E[0]);
      nue_m_hElepVsEv3->Fill(Enu,E[0]);
    }

    if(pdg[1] == 12){
      nue_ee_hElep_w0->Fill(E[0],pee);
      nue_ee_hElep_w3->Fill(E[0],pee);
      
      nue_em_hElep_w0->Fill(E[0],pmue);
      nue_em_hElep_w3->Fill(E[0],pmue);
      
      nue_e_hElepVsEv0->Fill(Enu,E[0]);
      nue_e_hElepVsEv3->Fill(Enu,E[0]);
    }
  }
  nue_me_hElep_w0->Scale(6.);
  nue_me_hElep_w3->Scale(6.);
  nue_em_hElep_w0->Scale(1/6.);
  nue_em_hElep_w3->Scale(1/6.);

  nue_e_hElep_w0->Add(nue_me_hElep_w0); nue_e_hElep_w0->Add(nue_ee_hElep_w0);
  nue_e_hElep_w3->Add(nue_me_hElep_w3); nue_e_hElep_w3->Add(nue_ee_hElep_w3);
  nue_m_hElep_w0->Add(nue_em_hElep_w0); nue_m_hElep_w0->Add(nue_mm_hElep_w0);
  nue_m_hElep_w3->Add(nue_em_hElep_w3); nue_m_hElep_w3->Add(nue_mm_hElep_w3);
  
      
 
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
  TCanvas *e_cElep_w = new TCanvas("e_cElep_w","",1600,600);
  e_cElep_w->Divide(2,1);
  e_cElep_w->cd(1);
  e_hsElep_w0->Draw();
  e_cElep_w->cd(2);
  e_hsElep_w3->Draw();
  e_cElep_w->SaveAs("e_Elep_reco_w_4.pdf");

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
  m_cElep_w->SaveAs("m_Elep_reco_w_4.pdf");
  
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
  cElepVsEv->SaveAs("ElepVsEv_4.pdf");

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
  out->Close();

}

