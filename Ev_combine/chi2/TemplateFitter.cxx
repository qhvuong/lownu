#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TGraph.h"

TemplateFitter::TemplateFitter(TH1D * CC_templates_m[480], TH1D * CC_templates_m_nc[480], TH1D * CC_templates_e[480], TH1D * nue_templates_m[480], TH1D * nue_templates_m_w[480], TH1D * nue_templates_e[480], TH1D * nue_templates_e_w[480], TH1D * CC_target_e, TH1D * CC_target_m, TH1D * nue_target)
{
  for( int i = 0; i < 480; ++i ) {
  CC_m_templates[i] = CC_templates_m[i];
  CC_nc_m_templates[i] = CC_templates_m_nc[i];
  CC_e_templates[i] = CC_templates_e[i];
  nue_m_templates[i] = nue_templates_m[i];
  nue_w_m_templates[i] = nue_templates_m_w[i];
  nue_e_templates[i] = nue_templates_e[i];
  nue_w_e_templates[i] = nue_templates_e_w[i];
  }
  
  CC_m_target = CC_target_m;
  CC_e_target = CC_target_e;
  target_nue = nue_target;
}

void TemplateFitter::setEnergyBins( double bins[481] )
{
  for( int i = 0; i < 481; ++ i ) {m_energy_bins[i] = bins[i];
  }
}

double TemplateFitter::getPmue( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mue2 = 4 * Uee2 * Umm2;
  double prob = s2mue2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPee( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2ee2 = 4 * Uee2 * (1 - Uee2);
  double prob = 1.0 - s2ee2  * pow(sin(del),2);
  return prob;
}
double TemplateFitter::getPmm( double energy, double Uee2, double Umm2, double dm2 )
{
  double L = 0.5;
  double del = 1.27*L*dm2/energy;
  double s2mm2 = 4 * Umm2 * (1 - Umm2);
  double prob = 1.0 - s2mm2  * pow(sin(del),2);
  return prob;
}

TMatrixD covmtr(300,300);

void TemplateFitter::setCovmtr( double bincontent[301][301] )
{
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      covmtr[i][j] = bincontent[i][j];
    }
  }
}

// function whose return Minuit mimizes, must take const double* and return double
double TemplateFitter::getChi2( double * par )
{ 
  // Create a histogram "temp" from the templates
  // Start with the intrinsic nu_e CC template, which doesn't change with oscillations
  TH1D * CC_tp_e = (TH1D*) CC_m_templates[0]->Clone();
  CC_tp_e->Reset();
  TH1D * CC_tp_me = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_ee = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_m = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_em = (TH1D*) CC_tp_e->Clone();
  TH1D * CC_tp_mm = (TH1D*) CC_tp_e->Clone();

  TH1D * nue_tp = (TH1D*) nue_m_templates[0]->Clone();
  nue_tp->Reset();
  TH1D * nue_tp_em = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_ee = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_me = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_mm = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_os = (TH1D*) nue_tp->Clone();
  TH1D * nue_tp_unos = (TH1D*) nue_tp->Clone();  

  // Add in oscillated neutrinos by taking the nu_mu CC templates and weighting by the oscillation probability
  for( int i = 10; i < 480; ++i ) {
    double mue = 0;
    double ee = 0;
    double mm = 0;
    for(int j = 0; j<1001; j++){
      double e = m_energy_bins[i] + j*(m_energy_bins[i+1] - m_energy_bins[i])/1000.;
      mue = mue + getPmue(e, par[0], par[1], par[2]); // par[0] = Uee2, par[1] = Umm2, par[2] = dm2
      ee  = ee  + getPee(e, par[0], par[1], par[2]);
      mm  = mm  + getPmm(e, par[0], par[1], par[2]);
    } 
    double Pmue = mue/1001.0; 
    double Pee  = ee/1001.0; 
    double Pmm  = mm/1001.0; 

    CC_tp_me->Add(CC_nc_m_templates[i], Pmue);
    CC_tp_ee->Add(CC_e_templates[i], Pee);

    CC_tp_em->Add(CC_e_templates[i], Pmue);
    CC_tp_mm->Add(CC_m_templates[i], Pmm);

    nue_tp_me->Add(nue_w_m_templates[i], Pmue);
    nue_tp_mm->Add(nue_m_templates[i], Pmm);
    nue_tp_em->Add(nue_w_e_templates[i], Pmue);
    nue_tp_ee->Add(nue_e_templates[i], Pee);
  }
  // Now we have nue temp = mu-->e (no reco cut) + e-->e (no reco cut)
  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em); CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);
   
  // calculate the chi2 with the "data" target
  double chi2 = 0.0;
  const int nbins_Elep = 100;

  TMatrixD target(3*nbins_Elep, 1);
  TMatrixD temp(3*nbins_Elep, 1);
  TMatrixD diff(3*nbins_Elep, 1);
  TMatrixD diff_T(1, 3*nbins_Elep);
  TMatrixD unc(3*nbins_Elep, 3*nbins_Elep);
  TMatrixD cov(3*nbins_Elep, 3*nbins_Elep);
  TMatrixD covmtr_tot(3*nbins_Elep, 3*nbins_Elep);

  for( int bx = 0; bx < nbins_Elep; bx++ ) {
    temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    temp[bx+nbins_Elep][0] = CC_tp_e->GetBinContent(bx+1);
    temp[bx+2*nbins_Elep][0] = nue_tp->GetBinContent(bx+1);

    target[bx][0] = CC_m_target->GetBinContent(bx+1);
    target[bx+nbins_Elep][0] = CC_e_target->GetBinContent(bx+1);
    target[bx+2*nbins_Elep][0] = target_nue->GetBinContent(bx+1);

    diff[bx][0] = temp[bx][0] - target[bx][0];
    diff[bx+nbins_Elep][0] = temp[bx+nbins_Elep][0] - target[bx+nbins_Elep][0];
    diff[bx+2*nbins_Elep][0] = temp[bx+2*nbins_Elep][0] - target[bx+2*nbins_Elep][0];

    diff_T[0][bx] = diff[bx][0];
    diff_T[0][bx+nbins_Elep] = diff[bx+nbins_Elep][0];
    diff_T[0][bx+2*nbins_Elep] = diff[bx+2*nbins_Elep][0];

    for( int by = 0; by < nbins_Elep; by++ ) {
      unc[bx][by]              = 0;
      unc[bx][by+nbins_Elep]   = 0;
      unc[bx][by+2*nbins_Elep] = 0;

      unc[bx+nbins_Elep][by]              = 0;
      unc[bx+nbins_Elep][by+nbins_Elep]   = 0;
      unc[bx+nbins_Elep][by+2*nbins_Elep] = 0;

      unc[bx+2*nbins_Elep][by]              = 0;
      unc[bx+2*nbins_Elep][by+nbins_Elep]   = 0;
      unc[bx+2*nbins_Elep][by+2*nbins_Elep] = 0;

      if(bx == by) {
        unc[bx][by]                           = CC_m_target->GetBinContent(bx+1);
        unc[bx+nbins_Elep][by+nbins_Elep]     = CC_e_target->GetBinContent(bx+1);
        unc[bx+2*nbins_Elep][by+2*nbins_Elep] = target_nue->GetBinContent(bx+1);
      }
/*
      cov[bx][by]              = covmtr[bx][by]              * CC_m_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx][by+nbins_Elep]   = covmtr[bx][by+nbins_Elep]   * CC_m_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx][by+2*nbins_Elep] = covmtr[bx][by+2*nbins_Elep] * CC_m_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+nbins_Elep][by]              = covmtr[bx+nbins_Elep][by]              * CC_e_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+nbins_Elep][by+nbins_Elep]   = covmtr[bx+nbins_Elep][by+nbins_Elep]   * CC_e_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+nbins_Elep][by+2*nbins_Elep] = covmtr[bx+nbins_Elep][by+2*nbins_Elep] * CC_e_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+2*nbins_Elep][by]              = covmtr[bx+2*nbins_Elep][by]              * target_nue->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+2*nbins_Elep][by+nbins_Elep]   = covmtr[bx+2*nbins_Elep][by+nbins_Elep]   * target_nue->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+2*nbins_Elep][by+2*nbins_Elep] = covmtr[bx+2*nbins_Elep][by+2*nbins_Elep] * target_nue->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);
*/
    }
  }

  covmtr_tot = unc;
  TDecompSVD svd(covmtr_tot);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, 3*nbins_Elep);
  for( int bx = 0; bx < 3*nbins_Elep; bx++ ) {
    for( int by = 0; by < 3*nbins_Elep; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }

  return chi2;
}


void TemplateFitter::Draw()
{
  int para = 2;
  int cutNo = 0;
  int cutEv = 0;
  int N = 40;
  TH2D *h0 = new TH2D("h0","",N,0,0.1, N,0,12.0);
  TH2D *h1 = new TH2D("h1","",N,0,12.0,N,0,0.1);
  TH2D *h2 = new TH2D("h2","",N,0,0.1, N,0,0.1);

  double seed1[3][3], seed2[3][3];
  seed1[0][0] = 0.01;
  seed1[0][1] = 0.0016;
  seed1[0][2] = 1.3;

  seed2[0][0] = 0.04;
  seed2[0][1] = 0.01;
  seed2[0][2] = 6.0;

  double p0[1], p1[1], p2[1]; //true values
  if(para == 1) {
    p0[0] = seed1[0][0];
    p1[0] = seed1[0][1];
    p2[0] = seed1[0][2];
  }
  else if(para == 2) {
    p0[0] = seed2[0][0];
    p1[0] = seed2[0][1];
    p2[0] = seed2[0][2];
  }

  double par[3], bin[3];
  par[0] = 0.0;
  par[1] = 0.0;
  par[2] = 0.0;
  double chi2 = 0.0;

  par[0] = p0[0];
  for(int j=1; j<=N; j++) {
    par[1] = h0->GetXaxis()->GetBinCenter(j);
    for(int k=1; k<=N; k++) {
      par[2] = h1->GetXaxis()->GetBinCenter(k);
      chi2 = getChi2(par);
      std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
      h0->Fill(par[1], par[2], chi2);
    }
  }

  par[1] = p1[0];
  for(int k=1; k<=N; k++) {
    par[2] = h1->GetXaxis()->GetBinCenter(k);
    for(int i=1; i<=N; i++) {
      par[0] = h2->GetXaxis()->GetBinCenter(i);
      chi2 = getChi2(par);
      std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
      h1->Fill(par[2], par[0], chi2);
    }
  }

  par[2] = p2[0];
  for(int i=1; i<=N; i++) {
    par[0] = h2->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[1] = h0->GetXaxis()->GetBinCenter(j);
      chi2 = getChi2(par);
      std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\n";
      h2->Fill(par[0], par[1], chi2);
    }
  }

  int n = 1;
  TGraph *g0 = new TGraph(n,p1,p2);
  TGraph *g1 = new TGraph(n,p2,p0);
  TGraph *g2 = new TGraph(n,p0,p1);

  g0->SetMarkerColor(kRed);
  g0->SetMarkerSize(5);
  g1->SetMarkerColor(kRed);
  g1->SetMarkerSize(5);
  g2->SetMarkerColor(kRed);
  g2->SetMarkerSize(5);

  h0->SetTitle("Chi2 Surface");
  h0->GetXaxis()->SetTitle("Umm2");
  h0->GetYaxis()->SetTitle("dm2");
  h0->SetStats(0);
  h1->SetTitle("Chi2 Surface");
  h1->GetXaxis()->SetTitle("dm2");
  h1->GetYaxis()->SetTitle("Uee2");
  h1->SetStats(0);
  h2->SetTitle("Chi2 Surface");
  h2->GetXaxis()->SetTitle("Uee2");
  h2->GetYaxis()->SetTitle("Umm2");
  h2->SetStats(0);
  
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
 
  TCanvas *c0 = new TCanvas("c0","",800,600);
  c0->SetLogz();
  h0->SetMaximum(1E6);
  h0->Draw("colz");
  g0->Draw("same C*");
  c0->SaveAs(Form("Ev_chi2Surface_%d%d0_%d.png",para,cutNo,cutEv));
  TCanvas *c1 = new TCanvas("c1","",800,600);
  c1->SetLogz();
  h1->SetMaximum(1E6);
  h1->Draw("colz");
  g1->Draw("same C*");
  c1->SaveAs(Form("Ev_chi2Surface_%d%d1_%d.png",para,cutNo,cutEv));
  TCanvas *c2 = new TCanvas("c2","",800,600);
  c2->SetLogz();
  h2->SetMaximum(1E6);
  h2->Draw("colz");
  g2->Draw("same C*");
  c2->SaveAs(Form("Ev_chi2Surface_%d%d2_%d.png",para,cutNo,cutEv));

  TFile *out = new TFile(Form("Ev_chi2_%d%d_%d.root",para,cutNo,cutEv),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  out->Close();
}
