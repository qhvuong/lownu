void test()
{
/*  TFile *f = new TFile("nue_output_2.root","READ");

  TH1D *hElep_w0_ft = (TH1D*)f->Get("hElep_w0_ft");
  TH1D *em_hElep_w0 = (TH1D*)f->Get("em_hElep_w0");
  TH1D *mm_hElep_w0 = (TH1D*)f->Get("mm_hElep_w0");
  TH1D *ee_hElep_w0 = (TH1D*)f->Get("ee_hElep_w0");
  TH1D *me_hElep_w0 = (TH1D*)f->Get("me_hElep_w0");

  TH1D *os = (TH1D*) em_hElep_w0->Clone();
  TH1D *unos = (TH1D*) ee_hElep_w0->Clone();

  os->Add(me_hElep_w0);
  unos->Add(mm_hElep_w0);

  //e_hElep_w0_ft->SetTitle(	)

  THStack *hs = new THStack("hs","oscillated + unoscillated weighted Elep_reco");
  os->SetMarkerStyle(21);
  os->SetMarkerSize(0.5);
  os->SetMarkerColor(kRed);
  os->GetXaxis()->SetTitle("Elep_reco");
  hs->Add(os);
  unos->SetMarkerStyle(21);
  unos->SetMarkerSize(0.5);
  unos->SetMarkerColor(kBlue);
  unos->GetXaxis()->SetTitle("Elep_reco");
  hs->Add(unos);
  TCanvas *c = new TCanvas("c","",1600,600);
  c->Divide(2,1);
  c->cd(1);
  gPad->SetLogy();
  hs->Draw();
  c->cd(2);
  gPad->SetLogy();
  hElep_w0_ft->SetTitle("fluctuated nu+e target");  
  hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");  
  hElep_w0_ft->Draw();  
  c->SaveAs("true_nue_Elep_reco_w_2_sum_Log.pdf");
*/


  TFile *f  = new TFile("nue_output_1.root","UPDATE");

  f->Delete("hElep_w0_ft;1");
  //f->Delete("hElep_w0_ft;2");

  TH1D *hElep_w0 = (TH1D*)f->Get("hElep_w0");
  TH1D *hElep_w0_ft = new TH1D("hElep_w0_ft","",100,0,16);
  //hElep_w0_flt->Reset();

  TRandom3 *rando = new TRandom3(8888);

  for( int bx = 1; bx <= hElep_w0->GetNbinsX(); bx++ ){
    double mean = hElep_w0->GetBinContent(bx);
    double fluctuated_bin_content = rando->Poisson(mean);
    hElep_w0_ft->AddBinContent(bx, fluctuated_bin_content);
  }

  hElep_w0_ft->SetStats(0);
  hElep_w0_ft->SetTitle("fluctuated nu+e target");
  hElep_w0_ft->GetXaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *c = new TCanvas("c","",800,600);
  hElep_w0_ft->Draw();
  c->SaveAs("true_nue_target_ft_1.pdf");

  hElep_w0_ft->Write();
  f->Close();

}
