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

TemplateFitter::TemplateFitter(TH1D * CC_templates_m[480], TH1D * CC_templates_m_nc[480], TH1D * CC_templates_e[480], TH1D * nue_templates_m[480], TH1D * nue_templates_m_w[480], TH1D * nue_templates_e[480], TH1D * nue_templates_e_w[480], TH1D * CC_target_m, TH1D * CC_target_e, TH1D * nue_target)
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

void TemplateFitter::setPara( char var[20], int oscpar, int nuCut, int EvCut, int seed )
{
  name  = var;
  para  = oscpar;
  cutNu = nuCut;
  cutEv = EvCut;
  s     = seed;
}

TMatrixD covmtr(300,300);

void TemplateFitter::setCovmtr( double bincontent[301][301] )
{
  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      covmtr[i][j] = bincontent[i][j];
    }
  }
  std::cout << covmtr[25][25] << "\n";
  std::cout << covmtr[245][123] << "\n";
  std::cout << covmtr[85][2] << "\n";
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
  TH1D * CC_tp_m  = (TH1D*) CC_tp_e->Clone();
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
  CC_tp_m->Add(CC_tp_em);    CC_tp_m->Add(CC_tp_mm);

  nue_tp_os->Add(nue_tp_me);    nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee);  nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os);       nue_tp->Add(nue_tp_unos);
   
  // calculate the chi2 with the "data" target
  const int nbins_E = 100;
  double chi2 = 0.0;

  TMatrixD target(3*nbins_E, 1);
  TMatrixD temp(3*nbins_E, 1);
  TMatrixD diff(3*nbins_E, 1);
  TMatrixD diff_T(1, 3*nbins_E);
  TMatrixD unc(3*nbins_E, 3*nbins_E);
  TMatrixD cov(3*nbins_E, 3*nbins_E);
  TMatrixD covmtr_tot(3*nbins_E, 3*nbins_E);

  for( int bx = 0; bx < nbins_E; bx++ ) {
    temp[bx][0] = CC_tp_m->GetBinContent(bx+1);
    temp[bx+nbins_E][0] = CC_tp_e->GetBinContent(bx+1);
    temp[bx+2*nbins_E][0] = nue_tp->GetBinContent(bx+1);

    target[bx][0] = CC_m_target->GetBinContent(bx+1);
    target[bx+nbins_E][0] = CC_e_target->GetBinContent(bx+1);
    target[bx+2*nbins_E][0] = target_nue->GetBinContent(bx+1);

    diff[bx][0] = temp[bx][0] - target[bx][0];
    diff[bx+nbins_E][0] = temp[bx+nbins_E][0] - target[bx+nbins_E][0];
    diff[bx+2*nbins_E][0] = temp[bx+2*nbins_E][0] - target[bx+2*nbins_E][0];

    diff_T[0][bx] = diff[bx][0];
    diff_T[0][bx+nbins_E] = diff[bx+nbins_E][0];
    diff_T[0][bx+2*nbins_E] = diff[bx+2*nbins_E][0];

    for( int by = 0; by < nbins_E; by++ ) {
      unc[bx][by]              = 0;
      unc[bx][by+nbins_E]   = 0;
      unc[bx][by+2*nbins_E] = 0;

      unc[bx+nbins_E][by]              = 0;
      unc[bx+nbins_E][by+nbins_E]   = 0;
      unc[bx+nbins_E][by+2*nbins_E] = 0;

      unc[bx+2*nbins_E][by]              = 0;
      unc[bx+2*nbins_E][by+nbins_E]   = 0;
      unc[bx+2*nbins_E][by+2*nbins_E] = 0;

      if(bx == by) {
        unc[bx][by]                           = CC_m_target->GetBinContent(bx+1);
        unc[bx+nbins_E][by+nbins_E]     = CC_e_target->GetBinContent(bx+1);
        unc[bx+2*nbins_E][by+2*nbins_E] = target_nue->GetBinContent(bx+1);
      }

      cov[bx][by]              = covmtr[bx][by]              * CC_m_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx][by+nbins_E]   = covmtr[bx][by+nbins_E]   * CC_m_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx][by+2*nbins_E] = covmtr[bx][by+2*nbins_E] * CC_m_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+nbins_E][by]              = covmtr[bx+nbins_E][by]              * CC_e_target->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+nbins_E][by+nbins_E]   = covmtr[bx+nbins_E][by+nbins_E]   * CC_e_target->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+nbins_E][by+2*nbins_E] = covmtr[bx+nbins_E][by+2*nbins_E] * CC_e_target->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);

      cov[bx+2*nbins_E][by]              = covmtr[bx+2*nbins_E][by]              * target_nue->GetBinContent(bx+1) * CC_m_target->GetBinContent(by+1);
      cov[bx+2*nbins_E][by+nbins_E]   = covmtr[bx+2*nbins_E][by+nbins_E]   * target_nue->GetBinContent(bx+1) * CC_e_target->GetBinContent(by+1);
      cov[bx+2*nbins_E][by+2*nbins_E] = covmtr[bx+2*nbins_E][by+2*nbins_E] * target_nue->GetBinContent(bx+1) * target_nue->GetBinContent(by+1);
    }
  }

  covmtr_tot = cov + unc;
  TDecompSVD svd(covmtr_tot);
  TMatrixD inv = svd.Invert();

  TMatrixD diff_cov(1, 3*nbins_E);
  for( int bx = 0; bx < 3*nbins_E; bx++ ) {
    for( int by = 0; by < 3*nbins_E; by++ ) {
      diff_cov[0][bx] += diff_T[0][by] * inv[by][bx];
    }
    chi2 += diff_cov[0][bx] * diff[bx][0];
  }

  return chi2;
}


void TemplateFitter::Draw()
{
  int N = 15;

  double p0[1], p1[1], p2[1];
  if(para == 1){
    p0[0] = 0.01;
    p1[0] = 0.0016;
    p2[0] = 1.3;
  }
  else if(para == 2){
    p0[0] = 0.04;
    p1[0] = 0.01;
    p2[0] = 6.0;
  }

  TH2D *h0 = new TH2D("h0","",N,0,0.1, N,0,12.0);
  TH2D *h1 = new TH2D("h1","",N,0,12.0,N,0,0.1);
  TH2D *h2 = new TH2D("h2","",N,0,0.1, N,0,0.1);
  TH2D *h0L = new TH2D("h0L","",N,0,0.1, N,0,12.0);
  TH2D *h1L = new TH2D("h1L","",N,0,12.0,N,0,0.1);
  TH2D *h2L = new TH2D("h2L","",N,0,0.1, N,0,0.1);
  TH2D *h0z = new TH2D("h0z","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1z = new TH2D("h1z","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2z = new TH2D("h2z","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  TH2D *h0d = new TH2D("h0d","",N,0,0.1, N,0,12.0);
  TH2D *h1d = new TH2D("h1d","",N,0,12.0,N,0,0.1);
  TH2D *h2d = new TH2D("h2d","",N,0,0.1, N,0,0.1);
  TH2D *h0dL = new TH2D("h0dL","",N,0,0.1, N,0,12.0);
  TH2D *h1dL = new TH2D("h1dL","",N,0,12.0,N,0,0.1);
  TH2D *h2dL = new TH2D("h2dL","",N,0,0.1, N,0,0.1);
  TH2D *h0dz = new TH2D("h0dz","",N,0.0      ,2.0*p1[0],N,0.5*p2[0],1.5*p2[0]);
  TH2D *h1dz = new TH2D("h1dz","",N,0.5*p2[0],1.5*p2[0],N,0.0      ,2.0*p0[0]);
  TH2D *h2dz = new TH2D("h2dz","",N,0.0      ,2.0*p0[0],N,0.0      ,2.0*p1[0]);

  double par[3], parz[3], bin[3];
  par[0] = p0[0];
  par[1] = p1[0];
  par[2] = p2[0];
  double chi2 = 0.0;
  double chi2z = 0.0;

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;

  if(cutEv == 0)      Ev = 3.0;
  else if(cutEv == 1) Ev = 0.8;
  else if(cutEv == 2) Ev = 0.5;

  double chi2t = getChi2(par);
  double diff, diffz;

  par[0]  = p0[0];
  parz[0] = p0[0];
  for(int j=1; j<=N; j++) {
    par[1]  = h0->GetXaxis()->GetBinCenter(j);
    parz[1] = h0z->GetXaxis()->GetBinCenter(j);
    for(int k=1; k<=N; k++) {
      par[2]  = h1->GetXaxis()->GetBinCenter(k);
      parz[2] = h1z->GetXaxis()->GetBinCenter(k);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( k==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << j*k*100./(3*N*N) << "%" << "\n";
      h0->Fill(par[1], par[2], chi2);
      h0L->Fill(par[1], par[2], chi2);
      h0z->Fill(parz[1], parz[2], chi2z);
      h0d->Fill(par[1], par[2], diff);
      h0dL->Fill(par[1], par[2], diff);
      h0dz->Fill(parz[1], parz[2], diffz);
    }
  }

  par[1]  = p1[0];
  parz[1] = p1[0];
  for(int k=1; k<=N; k++) {
    par[2]  = h1->GetXaxis()->GetBinCenter(k);
    parz[2] = h1z->GetXaxis()->GetBinCenter(k);
    for(int i=1; i<=N; i++) {
      par[0]  = h2->GetXaxis()->GetBinCenter(i);
      parz[0] = h2z->GetXaxis()->GetBinCenter(i);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff  = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( i==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (N*N+k*i)*100./(3*N*N) << "%" << "\n";
      h1->Fill(par[2], par[0], chi2);
      h1L->Fill(par[2], par[0], chi2);
      h1z->Fill(parz[2], parz[0], chi2z);
      h1d->Fill(par[2], par[0], diff);
      h1dL->Fill(par[2], par[0], diff);
      h1dz->Fill(parz[2], parz[0], diffz);
    }
  }

  par[2]  = p2[0];
  parz[2] = p2[0];
  for(int i=1; i<=N; i++) {
    par[0]  = h2->GetXaxis()->GetBinCenter(i);
    parz[0] = h2z->GetXaxis()->GetBinCenter(i);
    for(int j=1; j<=N; j++) {
      par[1]  = h0->GetXaxis()->GetBinCenter(j);
      parz[1] = h0z->GetXaxis()->GetBinCenter(j);
      chi2  = getChi2(par);
      chi2z = getChi2(parz);
      diff = sqrt(std::fabs(chi2 - chi2t));
      diffz = sqrt(std::fabs(chi2z - chi2t));
      if( j==N ) std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << chi2 << "\t" << diff << "\t" << (2*N*N+i*j)*100./(3*N*N) << "%" << "\n";
      h2->Fill(par[0], par[1], chi2);
      h2L->Fill(par[0], par[1], chi2);
      h2z->Fill(parz[0], parz[1], chi2z);
      h2d->Fill(par[0], par[1], diff);
      h2dL->Fill(par[0], par[1], diff);
      h2dz->Fill(parz[0], parz[1], diffz);
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

  h0->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0->GetXaxis()->SetTitle("Umm2");
  h0->GetYaxis()->SetTitle("dm2");
  h0->SetStats(0);
  h1->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1->GetXaxis()->SetTitle("dm2");
  h1->GetYaxis()->SetTitle("Uee2");
  h1->SetStats(0);
  h2->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2->GetXaxis()->SetTitle("Uee2");
  h2->GetYaxis()->SetTitle("Umm2");
  h2->SetStats(0);
  h0z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0z->GetXaxis()->SetTitle("Umm2");
  h0z->GetYaxis()->SetTitle("dm2");
  h0z->SetStats(0);
  h1z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1z->GetXaxis()->SetTitle("dm2");
  h1z->GetYaxis()->SetTitle("Uee2");
  h1z->SetStats(0);
  h2z->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2z->GetXaxis()->SetTitle("Uee2");
  h2z->GetYaxis()->SetTitle("Umm2");
  h2z->SetStats(0);
  h0L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0L->GetXaxis()->SetTitle("Umm2");
  h0L->GetYaxis()->SetTitle("dm2");
  h0L->SetStats(0);
  h1L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1L->GetXaxis()->SetTitle("dm2");
  h1L->GetYaxis()->SetTitle("Uee2");
  h1L->SetStats(0);
  h2L->SetTitle(Form("%s Chi2 Surface (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2L->GetXaxis()->SetTitle("Uee2");
  h2L->GetYaxis()->SetTitle("Umm2");
  h2L->SetStats(0);

  h0d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0d->GetXaxis()->SetTitle("Umm2");
  h0d->GetYaxis()->SetTitle("dm2");
  h0d->SetStats(0);
  h1d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1d->GetXaxis()->SetTitle("dm2");
  h1d->GetYaxis()->SetTitle("Uee2");
  h1d->SetStats(0);
  h2d->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2d->GetXaxis()->SetTitle("Uee2");
  h2d->GetYaxis()->SetTitle("Umm2");
  h2d->SetStats(0);
  h0dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0dz->GetXaxis()->SetTitle("Umm2");
  h0dz->GetYaxis()->SetTitle("dm2");
  h0dz->SetStats(0);
  h1dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1dz->GetXaxis()->SetTitle("dm2");
  h1dz->GetYaxis()->SetTitle("Uee2");
  h1dz->SetStats(0);
  h2dz->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2dz->GetXaxis()->SetTitle("Uee2");
  h2dz->GetYaxis()->SetTitle("Umm2");
  h2dz->SetStats(0);
  h0dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h0dL->GetXaxis()->SetTitle("Umm2");
  h0dL->GetYaxis()->SetTitle("dm2");
  h0dL->SetStats(0);
  h1dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h1dL->GetXaxis()->SetTitle("dm2");
  h1dL->GetYaxis()->SetTitle("Uee2");
  h1dL->SetStats(0);
  h2dL->SetTitle(Form("%s sqrt(Chi2-Chi2_true) (nu<%.1fGeV & Etheta2<%.1fMeV)",name,nu,Ev));
  h2dL->GetXaxis()->SetTitle("Uee2");
  h2dL->GetYaxis()->SetTitle("Umm2");
  h2dL->SetStats(0);

  double chi2L_max = 1E11;
  double chi2_max = 2E10;
  double chi2z_max = 6E7;
  double diffL_max = 1E6;
  double diff_max = 3E5;
  double diffz_max = 8E3;

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *cchi2_0L = new TCanvas("cchi2_0L","",800,600);
  cchi2_0L->SetLogz(1);
  h0L->SetMaximum(chi2L_max);
  h0L->Draw("colz");
  g0->Draw("same C*");
  cchi2_0L->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_0_Log.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_0 = new TCanvas("cchi2_0","",800,600);
  cchi2_0->SetLogz(0);
  h0->SetMaximum(chi2_max);
  h0->Draw("colz");
  g0->Draw("same C*");
  cchi2_0->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_0.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_0z = new TCanvas("cchi2_0z","",800,600);
  cchi2_0z->SetLogz(0);
  h0z->SetMaximum(chi2z_max);
  h0z->Draw("colz");
  g0->Draw("same C*");
  cchi2_0z->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_0_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cchi2_1L = new TCanvas("cchi2_1L","",800,600);
  cchi2_1L->SetLogz(1);
  h1L->SetMaximum(chi2L_max);
  h1L->Draw("colz");
  g1->Draw("same C*");
  cchi2_1L->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_1_Log.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_1 = new TCanvas("cchi2_1","",800,600);
  cchi2_1->SetLogz(0);
  h1->SetMaximum(chi2_max);
  h1->Draw("colz");
  g1->Draw("same C*");
  cchi2_1->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_1.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_1z = new TCanvas("cchi2_1z","",800,600);
  cchi2_1z->SetLogz(0);
  h1z->SetMaximum(chi2z_max);
  h1z->Draw("colz");
  g1->Draw("same C*");
  cchi2_1z->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_1_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cchi2_2L = new TCanvas("cchi2_2L","",800,600);
  cchi2_2L->SetLogz(1);
  h2L->SetMaximum(chi2L_max);
  h2L->Draw("colz");
  g2->Draw("same C*");
  cchi2_2L->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_2_Log.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_2 = new TCanvas("cchi2_2","",800,600);
  cchi2_2->SetLogz(0);
  h2->SetMaximum(chi2_max);
  h2->Draw("colz");
  g2->Draw("same C*");
  cchi2_2->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_2.png",name,para,cutNu,cutEv));
  TCanvas *cchi2_2z = new TCanvas("cchi2_2z","",800,600);
  cchi2_2z->SetLogz(0);
  h2z->SetMaximum(chi2z_max);
  h2z->Draw("colz");
  g2->Draw("same C*");
  cchi2_2z->SaveAs(Form("%s_chi2Surface_wCov_%d%d%d_2_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_0L = new TCanvas("cdiff_0L","",800,600);
  cdiff_0L->SetLogz(1);
  h0dL->SetMaximum(diffL_max);
  h0dL->Draw("colz");
  g0->Draw("same C*");
  cdiff_0L->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_0_Log.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_0 = new TCanvas("cdiff_0","",800,600);
  cdiff_0->SetLogz(0);
  h0d->SetMaximum(diff_max);
  h0d->Draw("colz");
  g0->Draw("same C*");
  cdiff_0->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_0.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_0z = new TCanvas("cdiff_0z","",800,600);
  cdiff_0z->SetLogz(0);
  h0dz->SetMaximum(diffz_max);
  h0dz->Draw("colz");
  g0->Draw("same C*");
  cdiff_0z->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_0_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_1L = new TCanvas("cdiff_1L","",800,600);
  cdiff_1L->SetLogz(1);
  h1dL->SetMaximum(diffL_max);
  h1dL->Draw("colz");
  g1->Draw("same C*");
  cdiff_1L->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_1_Log.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_1 = new TCanvas("cdiff_1","",800,600);
  cdiff_1->SetLogz(0);
  h1d->SetMaximum(diff_max);
  h1d->Draw("colz");
  g1->Draw("same C*");
  cdiff_1->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_1.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_1z = new TCanvas("cdiff_1z","",800,600);
  cdiff_1z->SetLogz(0);
  h1dz->SetMaximum(diffz_max);
  h1dz->Draw("colz");
  g1->Draw("same C*");
  cdiff_1z->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_1_zoom.png",name,para,cutNu,cutEv));

  TCanvas *cdiff_2L = new TCanvas("cdiff_2L","",800,600);
  cdiff_2L->SetLogz(1);
  h2dL->SetMaximum(diffL_max);
  h2dL->Draw("colz");
  g2->Draw("same C*");
  cdiff_2L->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_2_Log.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_2 = new TCanvas("cdiff_2","",800,600);
  cdiff_2->SetLogz(0);
  h2d->SetMaximum(diff_max);
  h2d->Draw("colz");
  g2->Draw("same C*");
  cdiff_2->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_2.png",name,para,cutNu,cutEv));
  TCanvas *cdiff_2z = new TCanvas("cdiff_2z","",800,600);
  cdiff_2z->SetLogz(0);
  h2dz->SetMaximum(diffz_max);
  h2dz->Draw("colz");
  g2->Draw("same C*");
  cdiff_2z->SaveAs(Form("%s_chi2Diff_wCov_%d%d%d_2_zoom.png",name,para,cutNu,cutEv));

  TFile *out = new TFile(Form("%s_chi2_wCov_%d%d%d.root",name,para,cutNu,cutEv),"RECREATE");
  h0->Write();
  h1->Write();
  h2->Write();
  h0L->Write();
  h1L->Write();
  h2L->Write();
  h0z->Write();
  h1z->Write();
  h2z->Write();
  h0d->Write();
  h1d->Write();
  h2d->Write();
  h0dL->Write();
  h1dL->Write();
  h2dL->Write();
  h0dz->Write();
  h1dz->Write();
  h2dz->Write();
  out->Close();

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

}
